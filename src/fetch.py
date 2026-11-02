"""Download viral RNA records from NCBI and build a provenance manifest.

Standard library only. Two modes:

* online  -- ``esearch`` + ``efetch`` (GenBank flat file, ``gbwithparts``),
             cached under ``data/raw/<group>_<accession>.gb``;
* offline -- ``--offline`` simply re-parses whatever ``.gb`` files are cached,
             so every downstream step stays reproducible without the network.

The manifest is the single source of truth for provenance: accession/version,
organism, taxid, molecule type, topology, segment, length, GC/N, annotated
5'UTR / 3'UTR / CDS spans, source database, release date, fetch date and an
md5 of the sequence.
"""

from __future__ import annotations

import glob
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

from common import ensure_dir, md5_text, to_int, write_fasta, write_tsv

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
SLEEP_SECONDS = 0.11 if API_KEY else 0.4
USER_AGENT = "viral-motif-demo/1.0 (educational prototype)"

# Query panel: five taxonomically distant groups so that "conserved across
# taxa" and "lineage-specific" signals can be told apart later.  Both RefSeq
# representatives and a few GenBank genomes per group are pulled on purpose:
# the duplicates that this creates are exactly what the QC/clustering step has
# to neutralise.
GROUPS = {
    "sars2": {
        "label": "Betacoronavirus (SARS-CoV-2)",
        "refseq": "txid2697049[Organism] AND srcdb_refseq[PROP]",
        # Curated complete genomes instead of a generic search: the generic
        # "complete genome" hits are N-heavy consensus sequences that QC
        # (correctly) rejects, which would leave no near-duplicates to demo the
        # homology-clustering step with.
        "genbank": (
            "txid2697049[Organism] AND (MN908947.3[ACCN] OR MN985325.1[ACCN] OR "
            "MT020781.1[ACCN] OR MT027064.1[ACCN] OR LC528232.1[ACCN] OR MT039890.1[ACCN])"
        ),
        "refseq_max": 1,
        "genbank_max": 6,
    },
    "hcv": {
        "label": "Hepacivirus C (HCV)",
        "refseq": "txid11103[Organism] AND srcdb_refseq[PROP]",
        "genbank": "txid11103[Organism] AND srcdb_genbank[PROP] AND 9000:9800[SLEN]",
        "refseq_max": 6,
        "genbank_max": 2,
    },
    "hiv1": {
        "label": "Human immunodeficiency virus 1",
        "refseq": "txid11676[Organism] AND srcdb_refseq[PROP]",
        "genbank": "txid11676[Organism] AND srcdb_genbank[PROP] AND 9000:9800[SLEN]",
        "refseq_max": 1,
        "genbank_max": 3,
    },
    "zika": {
        "label": "Zika virus",
        "refseq": "txid64320[Organism] AND srcdb_refseq[PROP]",
        "genbank": "txid64320[Organism] AND srcdb_genbank[PROP] AND 10000:11000[SLEN]",
        "refseq_max": 2,
        "genbank_max": 1,
    },
    "flua": {
        "label": "Influenza A virus (segmented)",
        "refseq": "txid11320[Organism] AND srcdb_refseq[PROP]",
        "genbank": "",
        "refseq_max": 16,
        "genbank_max": 0,
    },
}

MANIFEST_FIELDS = [
    "record_id", "accession", "accession_version", "organism", "taxid",
    "taxon_group", "mol_type", "topology", "segment", "length", "gc_pct",
    "n_count", "utr5_start", "utr5_end", "utr3_start", "utr3_end",
    "cds_start", "cds_end", "source_db", "release_date", "fetch_date",
    "md5", "notes",
]


# --------------------------------------------------------------------------- #
# HTTP helpers (throttled, retried, cached on disk)
# --------------------------------------------------------------------------- #

def _http_get(url: str, params: dict, attempts: int = 3) -> str:
    """GET with NCBI throttling and exponential backoff."""
    if API_KEY:
        params = dict(params, api_key=API_KEY)
    full = url + "?" + urllib.parse.urlencode(params)
    last_error = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(full, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=90) as response:
                return response.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last_error = exc
            time.sleep(SLEEP_SECONDS * (2 ** attempt) + 1.0)
    raise RuntimeError(f"NCBI request failed after {attempts} attempts: {last_error}")


def esearch(term: str, retmax: int):
    """Return up to ``retmax`` nuccore UIDs matching ``term``."""
    if not term or retmax <= 0:
        return []
    time.sleep(SLEEP_SECONDS)
    text = _http_get(
        f"{EUTILS}/esearch.fcgi",
        {"db": "nuccore", "term": term, "retmax": retmax, "retmode": "json"},
    )
    try:
        payload = json.loads(text)
        uids = payload.get("esearchresult", {}).get("idlist", [])
    except ValueError:
        uids = re.findall(r'"(?:id|uid)":"(\d+)"', text)
    return list(dict.fromkeys(str(uid) for uid in uids))[:retmax]


def efetch(uids):
    """Fetch GenBank flat-file records for ``uids``."""
    if not uids:
        return ""
    time.sleep(SLEEP_SECONDS)
    return _http_get(
        f"{EUTILS}/efetch.fcgi",
        {"db": "nuccore", "id": ",".join(uids), "rettype": "gbwithparts", "retmode": "text"},
    )


def _split_records(text: str):
    """Split a GenBank flat file into per-record blocks."""
    blocks = []
    current = []
    for line in text.splitlines():
        current.append(line)
        if line.strip() == "//":
            blocks.append("\n".join(current))
            current = []
    if current and any(chunk.strip() for chunk in current):
        blocks.append("\n".join(current))
    return [block for block in blocks if block.lstrip().startswith("LOCUS")]


def download(root: str, groups=None, offline: bool = False):
    """Fetch records for the requested groups into ``data/raw`` and return the dir."""
    raw_dir = ensure_dir(os.path.join(root, "data", "raw"))
    selected = [g.strip() for g in (groups or list(GROUPS)) if g.strip() in GROUPS]
    if offline:
        print("[fetch] offline mode -> reusing cached GenBank files")
        return raw_dir
    for key in selected:
        spec = GROUPS[key]
        plan = (
            ("refseq", spec.get("refseq", ""), spec.get("refseq_max", 0)),
            ("genbank", spec.get("genbank", ""), spec.get("genbank_max", 0)),
        )
        for source, term, limit in plan:
            if not term or limit <= 0:
                continue
            try:
                uids = esearch(term, limit)
            except RuntimeError as exc:
                print(f"[fetch] {key}/{source}: search failed ({exc})")
                continue
            print(f"[fetch] {key}/{source}: {len(uids)} record(s) selected")
            # NCBI's nuccore efetch returns only the first record when several
            # comma-separated UIDs are passed via GET, so fetch them one by one.
            for uid in uids:
                try:
                    text = efetch([uid])
                except RuntimeError as exc:
                    print(f"[fetch] {key}/{source}: fetch failed ({exc})")
                    continue
                for block in _split_records(text):
                    record = parse_record(block)
                    if not record:
                        continue
                    record["taxon_group"] = key
                    path = os.path.join(raw_dir, f"{key}_{record['accession_version']}.gb")
                    with open(path, "w", encoding="utf-8", newline="\n") as fh:
                        fh.write(block.rstrip() + "\n")
                    print(f"        saved {os.path.basename(path)} ({record['length']} bp, "
                          f"{record.get('mol_type') or '?'})")
    return raw_dir


# --------------------------------------------------------------------------- #
# GenBank flat-file parsing
# --------------------------------------------------------------------------- #

LOCUS_RE = re.compile(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp\s+(.*)$")
FEATURE_RE = re.compile(r"^ {5}(\S+)\s+(\S.*)?$")
QUALIFIER_RE = re.compile(r"^ {21}/([A-Za-z0-9_']+)=(?:\"([^\"]*)\"|(\S+))$")
CONTINUATION_RE = re.compile(r"^ {21}([0-9<>(),a-zA-Z]+)$")


def parse_location(loc: str):
    """Reduce a GenBank location string to ``(start, end, strand)``.

    Handles ``complement(...)``, ``join(a..b,c..d)``, partial markers (``<1..90``)
    and bare positions. For joins we keep the outer envelope, which is all the
    manifest needs (exon structure is not used by the k-mer analysis).
    """
    if not loc:
        return None
    strand = "-" if "complement" in loc else "+"
    spans = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc)]
    if spans:
        return min(s for s, _ in spans), max(e for _, e in spans), strand
    singles = [int(x) for x in re.findall(r"(?<![.\d])(\d+)(?![.\d])", loc)]
    if singles:
        return min(singles), max(singles), strand
    return None


def _parse_features(lines):
    """Parse the FEATURES table into ``{key: [{'loc':..., 'quals':{...}}, ...]}``."""
    features = {}
    current = None
    in_features = False
    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        stripped = line.strip()
        if stripped.startswith(("ORIGIN", "BASE COUNT", "CONTIG")):
            break
        if not stripped:
            continue
        match = FEATURE_RE.match(line)
        if match:
            current = {"loc": (match.group(2) or "").strip(), "quals": {}}
            features.setdefault(match.group(1), []).append(current)
            continue
        if current is None:
            continue
        qual = QUALIFIER_RE.match(line)
        if qual:
            value = qual.group(2) if qual.group(2) is not None else qual.group(3)
            current["quals"][qual.group(1)] = value
            continue
        if CONTINUATION_RE.match(line):
            current["loc"] += stripped
    return features


def _extract_taxid(db_xref: str) -> str:
    match = re.search(r"taxon:(\d+)", db_xref or "")
    return match.group(1) if match else ""


def parse_record(block: str):
    """Parse one GenBank flat-file record into a dict (or None if unparsable)."""
    lines = block.splitlines()
    if not lines:
        return None
    match = LOCUS_RE.match(lines[0])
    if not match:
        return None
    name = match.group(1)
    length = int(match.group(2))
    rest = match.group(3).split()
    release_date = rest[-1] if rest else ""
    division = rest[-2] if len(rest) >= 2 else ""
    middle = rest[:-2] if len(rest) >= 2 else []
    topology = next((token for token in middle if token in ("linear", "circular")), "")
    mol_type = " ".join(t for t in middle if t not in ("linear", "circular"))

    accession = version = definition = organism = ""
    for index, line in enumerate(lines):
        if line.startswith("ACCESSION") and not accession:
            accession = line.split()[1] if len(line.split()) > 1 else ""
        elif line.startswith("VERSION") and not version:
            version = line.split()[1] if len(line.split()) > 1 else ""
        elif line.startswith("DEFINITION") and not definition:
            definition = line[12:].strip()
            cursor = index + 1
            while cursor < len(lines) and lines[cursor].startswith(" " * 12):
                definition += " " + lines[cursor].strip()
                cursor += 1
        elif line.startswith("  ORGANISM") and not organism:
            organism = line[12:].strip()

    features = _parse_features(lines)
    source_quals = features.get("source", [{}])[0].get("quals", {})
    organism = organism or source_quals.get("organism", "")
    mol_type = source_quals.get("mol_type", mol_type)
    taxid = _extract_taxid(source_quals.get("db_xref", ""))
    segment = source_quals.get("segment", "") or source_quals.get("isolate", "")

    seq_chunks = []
    in_origin = False
    for line in lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            seq_chunks.append(re.sub(r"[^A-Za-z]", "", line))
    sequence = "".join(seq_chunks).upper()

    def span(keys):
        spans = []
        for key in keys:
            for feature in features.get(key, []):
                parsed = parse_location(feature["loc"])
                if parsed:
                    spans.append(parsed)
        if not spans:
            return None
        return min(s for s, _, _ in spans), max(e for _, e, _ in spans)

    utr5 = span(["5'UTR", "five_prime_UTR"])
    utr3 = span(["3'UTR", "three_prime_UTR"])
    cds = span(["CDS"])

    return {
        "accession": accession or name,
        "accession_version": version or name,
        "definition": definition,
        "organism": organism,
        "taxid": taxid,
        "mol_type": mol_type,
        "topology": topology,
        "segment": segment,
        "length": length,
        "release_date": release_date,
        "division": division,
        "taxon_group": "",
        "source_db": "refseq" if re.match(r"^N[CWZR]_", name) else "genbank",
        "utr5_start": utr5[0] if utr5 else "",
        "utr5_end": utr5[1] if utr5 else "",
        "utr3_start": utr3[0] if utr3 else "",
        "utr3_end": utr3[1] if utr3 else "",
        "cds_start": cds[0] if cds else "",
        "cds_end": cds[1] if cds else "",
        "n_features": sum(len(v) for v in features.values()),
        "seq": sequence,
    }


# --------------------------------------------------------------------------- #
# Manifest assembly
# --------------------------------------------------------------------------- #

def _gc_pct(seq: str) -> float:
    acgt = sum(seq.count(base) for base in "ACGT")
    if not acgt:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / acgt


def build_manifest(root: str):
    """Parse every cached ``.gb`` file into ``manifest.tsv`` + ``raw.fasta``."""
    raw_dir = os.path.join(root, "data", "raw")
    paths = sorted(glob.glob(os.path.join(raw_dir, "*.gb")))
    if not paths:
        # The GenBank caches are not committed; when they are absent keep the
        # existing manifest/raw.fasta rather than overwriting them with nothing.
        if os.path.exists(os.path.join(raw_dir, "manifest.tsv")):
            print("[fetch] no cached .gb files; keeping existing manifest.tsv / raw.fasta")
        else:
            print("[fetch] no cached .gb files; run without --offline to download")
        return []
    records = []
    for path in paths:
        group = os.path.basename(path).split("_", 1)[0]
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                record = parse_record(fh.read())
        except OSError as exc:
            print(f"[fetch] cannot read {os.path.basename(path)}: {exc}")
            continue
        if not record:
            print(f"[fetch] could not parse {os.path.basename(path)}")
            continue
        record["taxon_group"] = group
        record["fetch_date"] = time.strftime(
            "%Y-%m-%d", time.localtime(os.path.getmtime(path))
        )
        records.append(record)
    records.sort(key=lambda rec: (rec["taxon_group"], rec["accession_version"]))

    rows = []
    fasta = []
    for record in records:
        sequence = record["seq"]
        row = {field: record.get(field, "") for field in MANIFEST_FIELDS}
        row.update(
            {
                "record_id": record["accession_version"],
                "length": record["length"],
                "gc_pct": round(_gc_pct(sequence), 2),
                "n_count": sequence.count("N"),
                "md5": md5_text(sequence),
                "notes": f"features={record['n_features']};division={record['division']}",
            }
        )
        rows.append(row)
        fasta.append(
            (
                f"{record['accession_version']}|{record['taxon_group']}"
                f"|len={record['length']}|seg={record['segment']}",
                sequence,
            )
        )

    write_tsv(os.path.join(raw_dir, "manifest.tsv"), rows, MANIFEST_FIELDS)
    write_fasta(os.path.join(raw_dir, "raw.fasta"), fasta)
    total_bp = sum(to_int(row["length"]) for row in rows)
    print(f"[fetch] manifest: {len(rows)} record(s), {total_bp} bp "
          f"from {len(set(row['taxon_group'] for row in rows))} group(s)")
    return records


if __name__ == "__main__":  # pragma: no cover - thin convenience wrapper
    import argparse

    parser = argparse.ArgumentParser(description="Fetch viral RNA records from NCBI")
    parser.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    parser.add_argument("--groups", default=",".join(GROUPS))
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    download(args.root, args.groups.split(","), offline=args.offline)
    build_manifest(args.root)


