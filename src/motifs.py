"""k-mer enrichment ("motif discovery") on cleaned, cluster-weighted sequences.

The model is deliberately explicit so every number in the output TSV can be
traced back to a formula:

    obs_w = sum_i w_i * count_i(motif)          weighted observed occurrences
    slots = sum_i w_i * (len_i - k + 1)         weighted number of windows
    exp_w = slots * P_background(motif)         null expectation
    fold  = (obs_w + 0.5) / (exp_w + 0.5)       smoothed enrichment
    z     = (obs_w - exp_w) / sqrt(exp_w)
    p     = P(X >= obs_w), X ~ Poisson(exp_w)   right tail (regularised gamma)
    q     = Benjamini-Hochberg over all 4^k k-mers of that k

Everything else (cluster support, coverage, taxon breadth, family de-duplication) is
built on top of those five quantities.
"""

from __future__ import annotations

import math
import os
from itertools import product

import stats
from common import (
    ensure_dir,
    parse_header,
    read_fasta,
    read_tsv,
    to_float,
    write_tsv,
)

MOTIF_FIELDS = [
    "rank", "motif", "k", "gc_pct", "obs_weighted", "exp_weighted", "fold",
    "log2_fold", "z", "p_value", "q_value", "n_clusters_with", "coverage_pct",
    "n_taxa_with", "mean_copies_per_cluster", "family_of", "note",
]

ALL_FIELDS = [field for field in MOTIF_FIELDS if field not in ("rank", "family_of")]


def load_dataset(fasta_path: str, seqinfo_path: str):
    """Load cleaned sequences together with weights / groups / cluster ids."""
    info = {}
    if os.path.exists(seqinfo_path):
        info = {row["record_id"]: row for row in read_tsv(seqinfo_path)}
    records = []
    for header, sequence in read_fasta(fasta_path):
        name, fields = parse_header(header)
        meta = info.get(name, {})
        positional = fields.get("_pos", [])
        group = meta.get("taxon_group") or (positional[0] if positional else "unknown")
        records.append(
            {
                "id": name,
                "seq": sequence.upper(),
                "weight": to_float(meta.get("weight"), 1.0),
                "group": group,
                "cluster": meta.get("cluster_id") or name,
            }
        )
    if not records:
        raise ValueError(f"no sequences found in {fasta_path}")
    return records


def build_background(records, kind: str = "markov1"):
    """Estimate the background model from the (weighted) foreground sequences."""
    sequences = [record["seq"] for record in records]
    weights = [record["weight"] for record in records]
    model = stats.MarkovModel.from_sequences(sequences, weights)
    if kind == "composition":
        # Zero-order: drop the conditional part, keep only base frequencies.
        model.cond = {
            (left, right): model.mono[right]
            for left in stats.BASES
            for right in stats.BASES
        }
    elif kind != "markov1":
        raise ValueError(f"unknown background model: {kind}")
    return model


def _analyze_k(records, model, k: int):
    """Compute enrichment statistics for every k-mer of length ``k``."""
    observed = {}
    support = {}
    groups = {}
    for record in records:
        sequence = record["seq"]
        weight = record["weight"]
        if len(sequence) < k:
            continue
        counter = {}
        for start in range(len(sequence) - k + 1):
            mer = sequence[start : start + k]
            counter[mer] = counter.get(mer, 0) + 1
        for mer, count in counter.items():
            observed[mer] = observed.get(mer, 0.0) + weight * count
            support[mer] = support.get(mer, 0) + 1
            groups.setdefault(mer, set()).add(record["group"])

    slots = sum(
        record["weight"] * max(0, len(record["seq"]) - k + 1) for record in records
    )
    # Hypothesis family = every possible k-mer of this length, not just the ones
    # that happen to be observed (otherwise m is underestimated and FDR is off).
    candidates = ["".join(chars) for chars in product(stats.BASES, repeat=k)]

    expected_values = [model.expected_windows(mer, slots) for mer in candidates]
    p_values = []
    for mer, expected in zip(candidates, expected_values):
        seen = observed.get(mer, 0.0)
        p_values.append(stats.poisson_sf(seen, expected) if seen > 0.0 else 1.0)
    q_values = stats.bh_fdr(p_values)

    total_records = max(1, len(records))
    rows = []
    for mer, expected, p_value, q_value in zip(candidates, expected_values, p_values, q_values):
        seen = observed.get(mer)
        if not seen:
            continue
        safe_expected = max(expected, 1e-9)
        fold = (seen + 0.5) / (expected + 0.5)
        row = {
            "motif": mer,
            "k": k,
            "gc_pct": round(100.0 * sum(1 for ch in mer if ch in "GC") / k, 1),
            "obs_weighted": seen,
            "exp_weighted": expected,
            "fold": fold,
            "log2_fold": math.log2(fold) if fold > 0 else 0.0,
            "z": (seen - expected) / math.sqrt(safe_expected),
            "p_value": p_value,
            "q_value": q_value,
            "n_clusters_with": support.get(mer, 0),
            "coverage_pct": 100.0 * support.get(mer, 0) / total_records,
            "n_taxa_with": len(groups.get(mer, ())),
            "mean_copies_per_cluster": seen / support[mer] if support.get(mer) else 0.0,
            "note": "tiny_expected" if expected < 1e-6 else "",
        }
        rows.append(row)
    rows.sort(key=lambda item: (item["q_value"], -item["fold"]))
    return rows


# --------------------------------------------------------------------------- #
# Pipeline entry point
# --------------------------------------------------------------------------- #

def _evidence_score(row) -> float:
    """Heuristic browsing score (documented as NOT part of the FDR)."""
    fold = max(row["fold"], 1e-9)
    q_value = row["q_value"]
    q_strength = 0.0 if q_value <= 0 else min(-math.log10(q_value), 20.0)
    return math.log2(fold + 1.0) * q_strength * math.log2(1.0 + row["n_clusters_with"])


def run(
    root: str,
    ks=(4, 5, 6, 7, 8),
    bg: str = "markov1",
    min_clusters: int = 3,
    top: int = 50,
    fasta: str = "",
    seqinfo: str = "",
    results_dir: str = "",
    quiet: bool = False,
):
    """Run the k-mer enrichment analysis and write the result tables."""
    clean_dir = os.path.join(root, "data", "clean")
    fasta = fasta or os.path.join(clean_dir, "sequences.fasta")
    seqinfo = seqinfo or os.path.join(clean_dir, "seqinfo.tsv")
    results_dir = ensure_dir(results_dir or os.path.join(root, "results"))

    records = load_dataset(fasta, seqinfo)
    model = build_background(records, bg)
    all_rows = []
    for k in sorted(set(ks)):
        k_rows = _analyze_k(records, model, k)
        all_rows.extend(k_rows)
        if not quiet:
            best = k_rows[0] if k_rows else None
            message = f"[motifs] k={k}: {len(k_rows)} observed k-mer(s)"
            if best:
                message += (
                    f", best {best['motif']} (q={best['q_value']:.3g}, "
                    f"fold={best['fold']:.2f}x)"
                )
            print(message)

    selected = select_and_rank(all_rows, min_clusters=min_clusters, top=top)
    for row in selected:
        row["evidence_score"] = _evidence_score(row)
    for index, row in enumerate(selected, start=1):
        row["rank"] = index

    write_outputs(results_dir, selected, all_rows, records)
    if not quiet:
        print_top(selected)
    return {
        "selected": selected,
        "all_rows": all_rows,
        "records": records,
        "path": os.path.join(results_dir, "motifs.tsv"),
    }


def print_top(rows, limit: int = 15):
    if not rows:
        print("[motifs] no candidate passed the filters")
        return
    print(f"[motifs] top {min(limit, len(rows))} candidate motifs")
    header = f"  {'rank':>4} {'motif':<10} {'k':>2} {'fold':>7} {'q':>10} {'clusters':>8} {'taxa':>4}"
    print(header)
    for row in rows[:limit]:
        print(
            f"  {row['rank']:>4} {row['motif']:<10} {row['k']:>2} {row['fold']:>7.2f} "
            f"{row['q_value']:>10.2e} {row['n_clusters_with']:>8} {row['n_taxa_with']:>4}"
        )


def write_outputs(results_dir, selected, all_rows, records):
    """Write ``motifs.tsv`` (non-redundant top list) and ``motifs_all.tsv``."""
    write_tsv(os.path.join(results_dir, "motifs.tsv"), selected, MOTIF_FIELDS)
    ordered = sorted(all_rows, key=lambda row: (row["q_value"], -row["fold"]))
    write_tsv(os.path.join(results_dir, "motifs_all.tsv"), ordered, ALL_FIELDS)
    write_summary(os.path.join(results_dir, "summary.md"), selected, all_rows, records)


def write_summary(path, selected, all_rows, records):
    """Auto-generated result sheet (a short companion to results/motifs.tsv)."""
    total_bp = sum(len(record["seq"]) for record in records)
    groups = sorted({record["group"] for record in records})
    ks = sorted({row["k"] for row in all_rows})
    hypotheses = sum(4 ** k for k in ks)
    significant = sum(1 for row in all_rows if row["q_value"] < 0.05)
    lines = [
        "# 结果速览（自动生成）",
        "",
        f"- 输入序列（同源簇代表）：{len(records)} 条，共 {total_bp} bp",
        f"- 类群：{', '.join(groups)}",
        f"- 假设族：每个 k 单独做 BH 校正，族大小为 4^k（k = {ks[0]}..{ks[-1]}，合计 {hypotheses}）",
        f"- 观测到并输出的 k-mer：{len(all_rows)}，其中 FDR<0.05：{significant}",
        "",
        "## Top 候选 motif",
        "",
        "| rank | motif | k | fold | q (FDR) | clusters | taxa | coverage% | note |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in selected:
        lines.append(
            "| {rank} | `{motif}` | {k} | {fold:.2f} | {q:.2e} | {ncl} | {ntax} | "
            "{cov:.0f} | {note} |".format(
                rank=row["rank"],
                motif=row["motif"],
                k=row["k"],
                fold=row["fold"],
                q=row["q_value"],
                ncl=row["n_clusters_with"],
                ntax=row["n_taxa_with"],
                cov=row["coverage_pct"],
                note=row.get("note", ""),
            )
        )
    lines += [
        "",
        "> q_value 为 Benjamini-Hochberg FDR（假设族 = 同 k 的全部 4^k 个 k-mer），",
        "> fold 为相对一阶 Markov 自背景的富集倍数。榜上条目是统计候选，不是已验证的调控元件。",
        "",
    ]
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))


def select_and_rank(rows, min_clusters: int = 3, top: int = 50, family_pool: int = 3000):
    """Filter, sort, de-duplicate redundant families and return the top rows.

    ``ATGG`` / ``CTCCA``-style families would otherwise fill the whole top list
    because a short k-mer is always contained in its longer variants. Rows are
    grouped by substring containment and only the best-scoring member is kept.
    """
    filtered = [row for row in rows if row["n_clusters_with"] >= min_clusters]
    filtered.sort(key=lambda row: (row["q_value"], -row["fold"], -row["n_clusters_with"]))

    pool = filtered[: max(family_pool, top * 20)]
    representatives = []
    family_of = {}
    for row in pool:
        motif = row["motif"]
        parent = next((rep for rep in representatives if motif in rep or rep in motif), "")
        family_of[motif] = parent
        if not parent:
            representatives.append(motif)

    family_size = {}
    for motif, parent in family_of.items():
        if parent:
            family_size[parent] = family_size.get(parent, 0) + 1

    selected = []
    for row in filtered:
        motif = row["motif"]
        if motif not in family_of:
            continue
        if family_of[motif]:
            continue
        row = dict(row)
        members = family_size.get(motif, 0)
        if members:
            row["note"] = "; ".join(
                part for part in (row.get("note", ""), f"family_size={members}") if part
            )
        selected.append(row)
        if len(selected) >= top:
            break
    return selected
