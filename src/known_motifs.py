"""Sanity table: where do well-characterised RNA elements land in our ranking?

The pipeline ranks composition outliers. Validated regulatory elements are often
simply common in a viral genome (enrichment ~1) or selected against (fold < 1),
so looking them up explicitly makes that gap visible instead of hiding it.

Only elements whose function can be approximated by a short linear motif are
listed. Structure-dependent elements (HIV-1 TAR, HCV IRES, ...) cannot be reduced
to a fixed k-mer: a hexamer taken from the loop alone neither identifies the
element nor carries its activity. They are therefore left out of this table and
discussed as a limitation instead.

Usage: python src/cli.py known   (after `motifs`)
"""

from __future__ import annotations

import os

from common import ensure_dir, read_tsv, write_tsv

# element, core DNA (5'->3'), core RNA, function, scope, detection, structure context
KNOWN_ELEMENTS = [
    (
        "polyA signal", "AATAAA", "AAUAAA",
        "mRNA 3' 端加工",
        "真核转录本与多种病毒 mRNA",
        "正向 k-mer 计数", "无结构要求",
    ),
    (
        "ARE class I", "ATTTA", "AUUUA",
        "mRNA 稳定性（AU-rich element）",
        "真核转录本，常在 3'UTR",
        "正向 k-mer 计数", "无结构要求",
    ),
    (
        "Kozak context", "GCCACC", "GCCACC",
        "翻译起始上下文",
        "宿主 mRNA 与部分病毒 mRNA",
        "正向 k-mer 计数", "无结构要求",
    ),
    (
        "SARS-CoV-2 TRS-B core", "ACGAAC", "ACGAAC",
        "亚基因组转录（TRS 核心）",
        "SARS-CoV-2 及同属冠状病毒",
        "正向 k-mer 计数", "无结构要求",
    ),
    (
        "influenza poly(U) tract", "TTTTTT", "UUUUUU",
        "聚合酶打滑加 poly(A)，节段末端摆动区",
        "流感 A/B 病毒各节段末端",
        "正向 k-mer 计数", "无结构要求",
    ),
]

FIELDS = [
    "element", "core_dna", "core_rna", "function", "scope", "detection",
    "structure_context", "rank_within_k", "rank_global", "obs_weighted",
    "exp_weighted", "fold", "z", "q_value", "verdict",
]


def run(root: str):
    results_dir = os.path.join(root, "results")
    source = os.path.join(results_dir, "motifs_all.tsv")
    if not os.path.exists(source):
        raise FileNotFoundError(f"{source} not found - run `python src/cli.py motifs` first")
    rows = read_tsv(source)
    for row in rows:
        row["q_value"] = float(row["q_value"]) if row["q_value"] not in ("", "NA") else 1.0
        row["fold"] = float(row["fold"]) if row["fold"] not in ("", "NA") else 0.0

    ordered = sorted(rows, key=lambda item: (item["q_value"], -item["fold"]))
    global_rank = {}
    rank_within_k = {}
    counter = {}
    for position, row in enumerate(ordered, start=1):
        key = (row["motif"], row["k"])
        global_rank.setdefault(key, position)
        counter[row["k"]] = counter.get(row["k"], 0) + 1
        rank_within_k.setdefault(key, counter[row["k"]])
    index = {(row["motif"], row["k"]): row for row in rows}

    out = []
    for element, core_dna, core_rna, function, scope, detection, structure in KNOWN_ELEMENTS:
        key = (core_dna, str(len(core_dna)))
        base = {
            "element": element,
            "core_dna": core_dna,
            "core_rna": core_rna,
            "function": function,
            "scope": scope,
            "detection": detection,
            "structure_context": structure,
        }
        row = index.get(key)
        if row is None:
            out.append(
                {
                    **base,
                    "rank_within_k": "absent",
                    "rank_global": "absent",
                    "obs_weighted": 0,
                    "exp_weighted": 0,
                    "fold": 0,
                    "z": 0,
                    "q_value": "NA",
                    "verdict": "not observed in this dataset",
                }
            )
            continue
        fold = float(row["fold"])
        q_value = float(row["q_value"])
        z_value = float(row["z"])
        if fold >= 1.5 and q_value < 0.05:
            verdict = "enriched"
        elif z_value <= -3.0:
            # The pipeline reports right-tail p only, so depletion is read from z.
            verdict = f"depleted (z={z_value:.1f})"
        else:
            verdict = "not distinguishable"
        out.append(
            {
                **base,
                "rank_within_k": rank_within_k[key],
                "rank_global": global_rank[key],
                "obs_weighted": row["obs_weighted"],
                "exp_weighted": row["exp_weighted"],
                "fold": round(fold, 3),
                "z": round(z_value, 2),
                "q_value": row["q_value"],
                "verdict": verdict,
            }
        )

    ensure_dir(results_dir)
    output = os.path.join(results_dir, "known_motifs.tsv")
    write_tsv(output, out, FIELDS)
    print("[known] known RNA elements in the ranking")
    print(f"  {'element':<26} {'core':<8} {'fold':>6} {'q':>10} {'rank/k':>7}  verdict")
    for row in out:
        print(
            f"  {row['element']:<26} {row['core_dna']:<8} {row['fold']:>6} "
            f"{str(row['q_value']):>10} {str(row['rank_within_k']):>7}  {row['verdict']}"
        )
    print(f"[known] table -> {os.path.relpath(output, root)}")
    return out

