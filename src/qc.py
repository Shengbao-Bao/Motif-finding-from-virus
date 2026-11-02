"""Quality control, de-duplication and homology clustering.

The point of this module is that the *statistical unit* of the downstream
analysis is an independent sequence cluster, not a database record. Public
databases are massively redundant (one SARS-CoV-2 reference plus thousands of
near-identical genomes), and without correction every lineage-specific k-mer
would look "significantly enriched".

Three mechanisms neutralise that bias:

1. normalisation + hard filters (length, ambiguity, missing data);
2. exact and near-duplicate removal via 16-mer Jaccard single-linkage clustering;
3. a cluster weight ``1/cluster_size`` carried all the way into the statistics.
"""

from __future__ import annotations

import os

from common import (
    ensure_dir,
    md5_text,
    parse_header,
    read_fasta,
    read_tsv,
    to_float,
    to_int,
    write_fasta,
    write_tsv,
)

# IUPAC ambiguity codes mapped to a representative base (documented choice: it
# keeps the sequence usable while flagging how many positions were altered).
IUPAC = {
    "R": "A", "Y": "C", "S": "G", "W": "A", "K": "G",
    "M": "C", "B": "C", "D": "A", "H": "A", "V": "A",
}

QC_FIELDS = [
    "record_id", "accession", "organism", "taxon_group", "segment", "mol_type",
    "length_raw", "length_clean", "gc_pct", "n_count", "n_frac", "n_u_to_t",
    "n_iupac", "n_nonstd", "filter_status", "drop_reason", "duplicate_of",
    "cluster_id", "cluster_size", "weight", "is_representative",
]

SEQINFO_FIELDS = [
    "record_id", "accession", "organism", "taxid", "taxon_group", "segment",
    "mol_type", "source_db", "length", "gc_pct", "weight", "cluster_id",
    "cluster_size", "utr5_start", "utr5_end", "utr3_start", "utr3_end",
    "cds_start", "cds_end",
]

CLUSTER_FIELDS = [
    "cluster_id", "record_id", "accession", "taxon_group", "length", "weight",
    "representative", "cluster_size", "gc_pct",
]


def normalise(sequence: str, iupac: str = "map"):
    """Uppercase, drop junk, map U->T and (optionally) IUPAC codes.

    Returns ``(clean_sequence, counters)``; the counters make the QC report
    state *why* a sequence was altered instead of silently changing it.
    """
    counters = {"u_to_t": 0, "iupac": 0, "nonstd": 0}
    out = []
    for char in sequence.upper():
        if char in "ACGT":
            out.append(char)
        elif char == "U":
            out.append("T")
            counters["u_to_t"] += 1
        elif char == "N":
            out.append("N")
        elif char in IUPAC:
            out.append(IUPAC[char] if iupac == "map" else "N")
            counters["iupac"] += 1
        else:
            counters["nonstd"] += 1
    return "".join(out), counters


def n_profile(sequence: str):
    """Return ``(n_count, n_fraction, longest_n_run)``."""
    n_count = sequence.count("N")
    longest = current = 0
    for char in sequence:
        if char == "N":
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    fraction = n_count / len(sequence) if sequence else 1.0
    return n_count, fraction, longest


def gc_pct(sequence: str) -> float:
    acgt = sum(sequence.count(base) for base in "ACGT")
    if not acgt:
        return 0.0
    return 100.0 * (sequence.count("G") + sequence.count("C")) / acgt


def _quantile(sorted_values, q: float) -> float:
    """Linear-interpolated quantile of an already sorted list."""
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = (len(sorted_values) - 1) * q
    low = int(position)
    high = min(low + 1, len(sorted_values) - 1)
    weight = position - low
    return sorted_values[low] * (1 - weight) + sorted_values[high] * weight


def iqr_bounds(lengths):
    """Tukey fences ``(q1 - 1.5*IQR, q3 + 1.5*IQR)``.

    Returns ``None`` when the spread is degenerate (all lengths equal), since a
    zero-width fence would reject everything. Lengths are compared as floats
    directly instead of via ``10**log10(bound)``, which is where rounding used
    to flag a 12000 bp genome as an outlier.
    """
    values = sorted(float(x) for x in lengths)
    if len(values) < 4:
        return None
    q1 = _quantile(values, 0.25)
    q3 = _quantile(values, 0.75)
    spread = q3 - q1
    if spread <= 0:
        return None
    return q1 - 1.5 * spread, q3 + 1.5 * spread


def kmer_set(sequence: str, k: int = 16):
    """Set of hashed k-mers - hashing keeps the memory footprint small."""
    if len(sequence) < k:
        return set()
    return {hash(sequence[i : i + k]) for i in range(len(sequence) - k + 1)}


def cluster_sequences(records, k: int = 16, threshold: float = 0.9):
    """Single-linkage clustering on 16-mer Jaccard similarity.

    Two shortcuts keep this fast enough for a demo:

    * ``Jaccard <= min(|A|,|B|) / max(|A|,|B|)`` bounds the achievable score, so
      differently sized sequences (e.g. 30 kb vs 2 kb) are skipped cheaply;
    * the exact intersection is only computed for the surviving pairs.
    """
    sets = [kmer_set(record["seq"], k) for record in records]
    parent = list(range(len(records)))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    comparisons = 0
    for i in range(len(records)):
        set_i = sets[i]
        if not set_i:
            continue
        for j in range(i + 1, len(records)):
            set_j = sets[j]
            if not set_j:
                continue
            upper = min(len(set_i), len(set_j)) / max(len(set_i), len(set_j))
            if upper < threshold:
                continue
            comparisons += 1
            shared = len(set_i & set_j)
            jaccard = shared / (len(set_i) + len(set_j) - shared)
            if jaccard >= threshold:
                union(i, j)

    groups = {}
    for index in range(len(records)):
        groups.setdefault(find(index), []).append(index)
    return list(groups.values()), comparisons


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #

def run(
    root: str,
    min_len: int = 100,
    max_len: int = 50000,
    max_n_frac: float = 0.05,
    max_n_run: int = 10,
    iupac: str = "map",
    use_iqr: bool = False,
    cluster_k: int = 16,
    cluster_threshold: float = 0.9,
):
    """Run the QC pipeline and write ``data/clean/*``; return a summary dict."""
    raw_dir = os.path.join(root, "data", "raw")
    clean_dir = ensure_dir(os.path.join(root, "data", "clean"))
    manifest_path = os.path.join(raw_dir, "manifest.tsv")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"{manifest_path} not found - run `python src/cli.py fetch` first"
        )
    meta_by_id = {row["record_id"]: row for row in read_tsv(manifest_path)}

    records = []
    for header, sequence in read_fasta(os.path.join(raw_dir, "raw.fasta")):
        name, fields = parse_header(header)
        meta = meta_by_id.get(name, {})
        clean, counters = normalise(sequence, iupac=iupac)
        n_count, n_fraction, longest_run = n_profile(clean)
        record = {
            "record_id": name,
            "meta": meta,
            "seq": clean,
            "length_raw": len(sequence),
            "length_clean": len(clean),
            "gc_pct": round(gc_pct(clean), 2),
            "n_count": n_count,
            "n_frac": n_fraction,
            "longest_n_run": longest_run,
            "n_u_to_t": counters["u_to_t"],
            "n_iupac": counters["iupac"],
            "n_nonstd": counters["nonstd"],
            "status": "kept",
            "reason": "",
            "duplicate_of": "",
            "cluster_id": "",
            "cluster_size": 0,
            "weight": 0.0,
            "is_representative": False,
        }
        for condition, reason in (
            (record["length_clean"] < min_len, "too_short"),
            (record["length_clean"] > max_len, "too_long"),
            (n_fraction > max_n_frac, "n_fraction>5%"),
            (longest_run >= max_n_run, f"n_run>={max_n_run}"),
            (record["length_clean"] == 0, "empty"),
        ):
            if condition:
                record["status"] = "dropped"
                record["reason"] = reason
                break
        records.append(record)

    # Optional Tukey length filter. Off by default: this dataset deliberately
    # mixes 30 kb coronavirus genomes with 0.9-2.3 kb influenza segments, and a
    # global IQR assumption would silently discard the segmented virus.
    if use_iqr:
        survivors = [rec for rec in records if rec["status"] == "kept"]
        bounds = iqr_bounds([rec["length_clean"] for rec in survivors])
        if bounds:
            low, high = bounds
            for record in survivors:
                if record["length_clean"] < low or record["length_clean"] > high:
                    record["status"] = "dropped"
                    record["reason"] = "iqr_length_outlier"

    survivors = [rec for rec in records if rec["status"] == "kept"]

    # Exact duplicates: keep the RefSeq representative (then the longest).
    survivors.sort(
        key=lambda rec: (
            0 if rec["meta"].get("source_db") == "refseq" else 1,
            -rec["length_clean"],
        )
    )
    seen = {}
    for record in survivors:
        digest = md5_text(record["seq"])
        record["md5"] = digest
        if digest in seen:
            record["status"] = "dropped"
            record["reason"] = "exact_duplicate"
            record["duplicate_of"] = seen[digest]
        else:
            seen[digest] = record["record_id"]
    survivors = [rec for rec in records if rec["status"] == "kept"]

    clusters, comparisons = cluster_sequences(survivors, cluster_k, cluster_threshold)
    clusters.sort(key=lambda members: -max(survivors[i]["length_clean"] for i in members))
    representatives = []
    for number, members in enumerate(clusters, start=1):
        cluster_id = f"C{number:02d}"
        weight = 1.0 / len(members)
        best = min(
            members,
            key=lambda i: (
                0 if survivors[i]["meta"].get("source_db") == "refseq" else 1,
                -survivors[i]["length_clean"],
            ),
        )
        for index in members:
            record = survivors[index]
            record["cluster_id"] = cluster_id
            record["cluster_size"] = len(members)
            record["weight"] = weight
            record["is_representative"] = index == best
        representatives.append(survivors[best])

    write_qc_outputs(root, clean_dir, records, survivors, representatives)
    summary = {
        "records_total": len(records),
        "records_kept": len(survivors),
        "records_dropped": len(records) - len(survivors),
        "clusters": len(clusters),
        "jaccard_comparisons": comparisons,
        "bp_kept": sum(rec["length_clean"] for rec in survivors),
        "bp_effective": sum(rec["length_clean"] * rec["weight"] for rec in survivors),
    }
    print(
        f"[qc] {summary['records_total']} records -> {summary['records_kept']} kept "
        f"in {summary['clusters']} homology cluster(s); "
        f"{summary['bp_kept']} bp before clustering -> {summary['bp_effective']:.0f} bp "
        f"after collapsing redundant lineages"
    )
    reasons = {}
    for record in records:
        if record["status"] != "kept":
            reasons[record["reason"]] = reasons.get(record["reason"], 0) + 1
    for reason, count in sorted(reasons.items(), key=lambda item: -item[1]):
        print(f"       dropped: {reason} x{count}")
    return summary


def write_qc_outputs(root, clean_dir, records, survivors, representatives):
    """Write ``qc_report.tsv``, ``clusters.tsv``, ``seqinfo.tsv`` and the FASTA."""
    qc_rows = []
    for record in records:
        meta = record["meta"]
        qc_rows.append(
            {
                "record_id": record["record_id"],
                "accession": meta.get("accession_version", record["record_id"]),
                "organism": meta.get("organism", ""),
                "taxon_group": meta.get("taxon_group", ""),
                "segment": meta.get("segment", ""),
                "mol_type": meta.get("mol_type", ""),
                "length_raw": record["length_raw"],
                "length_clean": record["length_clean"],
                "gc_pct": record["gc_pct"],
                "n_count": record["n_count"],
                "n_frac": round(record["n_frac"], 4),
                "n_u_to_t": record["n_u_to_t"],
                "n_iupac": record["n_iupac"],
                "n_nonstd": record["n_nonstd"],
                "filter_status": record["status"],
                "drop_reason": record["reason"],
                "duplicate_of": record["duplicate_of"],
                "cluster_id": record["cluster_id"],
                "cluster_size": record["cluster_size"],
                "weight": round(record["weight"], 4),
                "is_representative": "yes" if record["is_representative"] else "",
            }
        )
    write_tsv(os.path.join(clean_dir, "qc_report.tsv"), qc_rows, QC_FIELDS)

    cluster_rows = []
    for record in sorted(survivors, key=lambda rec: (rec["cluster_id"], rec["record_id"])):
        meta = record["meta"]
        cluster_rows.append(
            {
                "cluster_id": record["cluster_id"],
                "record_id": record["record_id"],
                "accession": meta.get("accession_version", record["record_id"]),
                "taxon_group": meta.get("taxon_group", ""),
                "length": record["length_clean"],
                "weight": round(record["weight"], 4),
                "representative": "yes" if record["is_representative"] else "",
                "cluster_size": record["cluster_size"],
                "gc_pct": record["gc_pct"],
            }
        )
    write_tsv(os.path.join(clean_dir, "clusters.tsv"), cluster_rows, CLUSTER_FIELDS)

    seqinfo_rows = []
    fasta = []
    for record in sorted(representatives, key=lambda rec: rec["cluster_id"]):
        meta = record["meta"]
        # Downstream statistical unit = one independent lineage, so every cluster
        # representative enters the analysis with weight 1.0 ("one lineage, one
        # vote").  The within-cluster 1/cluster_size weight stays in
        # qc_report.tsv / clusters.tsv purely as redundancy documentation.
        unit_weight = 1.0
        seqinfo_rows.append(
            {
                "record_id": record["record_id"],
                "accession": meta.get("accession_version", record["record_id"]),
                "organism": meta.get("organism", ""),
                "taxid": meta.get("taxid", ""),
                "taxon_group": meta.get("taxon_group", ""),
                "segment": meta.get("segment", ""),
                "mol_type": meta.get("mol_type", ""),
                "source_db": meta.get("source_db", ""),
                "length": record["length_clean"],
                "gc_pct": record["gc_pct"],
                "weight": unit_weight,
                "cluster_id": record["cluster_id"],
                "cluster_size": record["cluster_size"],
                "utr5_start": meta.get("utr5_start", ""),
                "utr5_end": meta.get("utr5_end", ""),
                "utr3_start": meta.get("utr3_start", ""),
                "utr3_end": meta.get("utr3_end", ""),
                "cds_start": meta.get("cds_start", ""),
                "cds_end": meta.get("cds_end", ""),
            }
        )
        header = (
            f"{record['record_id']}|{meta.get('taxon_group','')}"
            f"|len={record['length_clean']}|seg={meta.get('segment','')}"
            f"|cluster={record['cluster_id']}|w={unit_weight:.1f}"
        )
        fasta.append((header, record["seq"]))

    write_tsv(os.path.join(clean_dir, "seqinfo.tsv"), seqinfo_rows, SEQINFO_FIELDS)
    write_fasta(os.path.join(clean_dir, "sequences.fasta"), fasta)



