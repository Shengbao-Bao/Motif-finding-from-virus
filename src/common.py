"""Shared IO helpers for the viral RNA motif pipeline (standard library only).

Everything here is deliberately simple so that the whole pipeline can be read
and explained line by line.
"""

from __future__ import annotations

import csv
import hashlib
import os

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def ensure_dir(path: str) -> str:
    """Create ``path`` (and parents) if missing and return it."""
    if path:
        os.makedirs(path, exist_ok=True)
    return path


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path: str):
    """Yield ``(header, sequence)`` tuples. Header excludes the leading '>'."""
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(path: str, records, width: int = 60) -> None:
    """Write ``[(header, sequence), ...]`` wrapped at ``width`` columns."""
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for header, seq in records:
            fh.write(">" + header + "\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def read_tsv(path: str):
    """Read a TSV into a list of dicts (values stay strings)."""
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: str, rows, fields) -> None:
    """Write dicts to TSV, keeping the column order given by ``fields``."""
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=list(fields), delimiter="\t", extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _fmt(row.get(k)) for k in fields})


def _fmt(value):
    """Stable, compact number formatting for TSV output."""
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:  # NaN
            return "NA"
        return f"{value:.6g}"
    return value


def parse_header(header: str):
    """Split a FASTA header into ``(name, fields_dict)``.

    Header grammar used by this project::

        NC_045512.2|sars2|len=29903|mol=genomic RNA|cluster=C1

    The first token is the record name; tokens containing '=' become
    key/value pairs; remaining tokens are collected under ``_pos``.
    """
    parts = header.split("|")
    name = parts[0].strip()
    fields = {"_pos": []}
    for token in parts[1:]:
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, _, value = token.partition("=")
            fields[key.strip()] = value.strip()
        else:
            fields["_pos"].append(token)
    return name, fields


def md5_text(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def md5_file(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def to_float(value, default=0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int(value, default=0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default
