"""Single entry point for the viral RNA motif pipeline.

Usage
-----
    python src/cli.py fetch    [--groups sars2,hcv,hiv1,zika,flua] [--offline]
    python src/cli.py qc       [--iqr] [--iupac map|drop] [--min-len 100]
    python src/cli.py motifs   [--k 4,5,6,7,8] [--bg markov1] [--top 50] [--fasta FILE]
    python src/cli.py known
    python src/cli.py all      [--refresh]

Standard library only. Every run drops ``results/run_metadata.json`` with Python/OS
versions, argument values, and script/input md5s.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DEFAULT = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import fetch  # noqa: E402
import known_motifs  # noqa: E402
import motifs  # noqa: E402
import qc  # noqa: E402
from common import ensure_dir, md5_file  # noqa: E402

INPUT_FILES = [
    "data/raw/manifest.tsv",
    "data/raw/raw.fasta",
    "data/clean/sequences.fasta",
    "data/clean/seqinfo.tsv",
    "data/clean/clusters.tsv",
]
SCRIPT_FILES = [
    "src/common.py",
    "src/stats.py",
    "src/fetch.py",
    "src/qc.py",
    "src/motifs.py",
    "src/known_motifs.py",
    "src/cli.py",
]


def parse_ints(text: str):
    return tuple(sorted({int(part) for part in str(text).split(",") if part.strip()}))


def write_metadata(root: str, args, timings, extra=None):
    """Record everything needed to reproduce (and audit) this run."""
    results_dir = ensure_dir(os.path.join(root, "results"))
    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "command": args.command,
        "arguments": {
            key: value
            for key, value in vars(args).items()
            if isinstance(value, (str, int, float, bool, type(None)))
        },
        "timings_seconds": {key: round(value, 3) for key, value in timings.items()},
        "inputs_md5": {},
        "scripts_md5": {},
    }
    for relative in INPUT_FILES:
        path = os.path.join(root, relative)
        if os.path.exists(path):
            payload["inputs_md5"][relative] = md5_file(path)
    for relative in SCRIPT_FILES:
        path = os.path.join(root, relative)
        if os.path.exists(path):
            payload["scripts_md5"][relative] = md5_file(path)
    if extra:
        payload.update(extra)
    path = os.path.join(results_dir, "run_metadata.json")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    return path


def cmd_fetch(args, timings):
    start = time.time()
    fetch.download(args.root, args.groups.split(","), offline=args.offline)
    fetch.build_manifest(args.root)
    timings["fetch"] = time.time() - start


def cmd_qc(args, timings):
    start = time.time()
    qc.run(
        args.root,
        min_len=args.min_len,
        max_len=args.max_len,
        max_n_frac=args.max_n_frac,
        iupac=args.iupac,
        use_iqr=args.iqr,
        cluster_threshold=args.cluster_threshold,
    )
    timings["qc"] = time.time() - start


def cmd_motifs(args, timings):
    start = time.time()
    motifs.run(
        args.root,
        ks=parse_ints(args.k),
        bg=args.bg,
        min_clusters=args.min_clusters,
        top=args.top,
        fasta=getattr(args, "fasta", "") or "",
        seqinfo=getattr(args, "seqinfo", "") or "",
        results_dir=getattr(args, "results_dir", "") or "",
    )
    timings["motifs"] = time.time() - start


def cmd_known(args, timings):
    start = time.time()
    known_motifs.run(args.root)
    timings["known"] = time.time() - start


def cmd_all(args, timings):
    manifest = os.path.join(args.root, "data", "raw", "manifest.tsv")
    args.offline = not args.refresh and os.path.exists(manifest)
    cmd_fetch(args, timings)
    cmd_qc(args, timings)
    cmd_motifs(args, timings)
    cmd_known(args, timings)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="cli.py",
        description="Candidate regulatory motif discovery in viral RNA (k-mer enrichment)",
    )
    parser.add_argument(
        "--root", default=ROOT_DEFAULT, help="project root (default: repository root)"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    fetch_parser = sub.add_parser("fetch", help="download GenBank records + build the manifest")
    fetch_parser.add_argument("--groups", default=",".join(fetch.GROUPS))
    fetch_parser.add_argument("--offline", action="store_true")
    fetch_parser.set_defaults(func=cmd_fetch)

    qc_parser = sub.add_parser("qc", help="normalise, filter, de-duplicate, cluster")
    qc_parser.add_argument("--min-len", type=int, default=100)
    qc_parser.add_argument("--max-len", type=int, default=50000)
    qc_parser.add_argument("--max-n-frac", type=float, default=0.05)
    qc_parser.add_argument("--iupac", choices=("map", "drop"), default="map")
    qc_parser.add_argument(
        "--iqr", action="store_true", help="enable the Tukey length filter (off by default)"
    )
    qc_parser.add_argument("--cluster-threshold", type=float, default=0.9)
    qc_parser.set_defaults(func=cmd_qc)

    motifs_parser = sub.add_parser("motifs", help="k-mer enrichment and candidate ranking")
    motifs_parser.add_argument("--k", default="4,5,6,7,8")
    motifs_parser.add_argument("--bg", choices=("markov1", "composition"), default="markov1")
    motifs_parser.add_argument("--min-clusters", type=int, default=3)
    motifs_parser.add_argument("--top", type=int, default=50)
    motifs_parser.add_argument(
        "--fasta", default="",
        help="analyse an external FASTA (e.g. examples/sample.fasta) instead of data/clean/sequences.fasta",
    )
    motifs_parser.add_argument(
        "--seqinfo", default="",
        help="optional TSV with record_id / taxon_group / weight / cluster_id",
    )
    motifs_parser.add_argument("--results-dir", default="", help="output directory (default: results/)")
    motifs_parser.set_defaults(func=cmd_motifs)

    known_parser = sub.add_parser("known", help="look up known RNA elements in the ranking")
    known_parser.set_defaults(func=cmd_known)

    all_parser = sub.add_parser("all", help="fetch -> qc -> motifs -> known")
    all_parser.add_argument("--groups", default=",".join(fetch.GROUPS))
    all_parser.add_argument("--offline", action="store_true")
    all_parser.add_argument("--refresh", action="store_true", help="force re-download")
    all_parser.add_argument("--k", default="4,5,6,7,8")
    all_parser.add_argument("--bg", choices=("markov1", "composition"), default="markov1")
    all_parser.add_argument("--min-clusters", type=int, default=3)
    all_parser.add_argument("--top", type=int, default=50)
    all_parser.add_argument("--min-len", type=int, default=100)
    all_parser.add_argument("--max-len", type=int, default=50000)
    all_parser.add_argument("--max-n-frac", type=float, default=0.05)
    all_parser.add_argument("--iupac", choices=("map", "drop"), default="map")
    all_parser.add_argument("--iqr", action="store_true")
    all_parser.add_argument("--cluster-threshold", type=float, default=0.9)
    all_parser.set_defaults(func=cmd_all)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.root = os.path.abspath(args.root)
    timings = {}
    started = time.time()
    args.func(args, timings)
    timings["total"] = time.time() - started
    path = write_metadata(args.root, args, timings)
    print(f"[cli] {args.command} done in {timings['total']:.1f}s; metadata -> "
          f"{os.path.relpath(path, args.root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


