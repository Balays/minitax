#!/usr/bin/env python3
"""
Build a contig-level custom minitax database from MAG FASTA files and a taxonomy table.

The input genome FASTA files may contain arbitrary contig headers. This script rewrites every
contig header to a unique minitax-safe sequence ID that contains the parent genome ID, while the
output taxonomy table maps every contig back to the same genome-level taxonomy.

Outputs:
  - <prefix>.fna or <prefix>.fna.gz
  - db_data.tsv
  - MAG.db.tsv
  - MAG.db.uni.tsv
  - contig_manifest.tsv
  - build_manifest.tsv
  - unmatched_genomes.tsv
  - missing_files.tsv

Example:
  scripts/build_mag_contig_db_for_minitax.py \
    --genomes-dir /path/to/dereplicated_genomes \
    --taxonomy-tsv mag_taxonomy_gtdbtk_summary.tsv \
    --outdir /path/to/minitax_mag_db \
    --prefix CanMAG_contigs \
    --gzip \
    --build-index \
    --minimap2 mm2-fast \
    --threads 20 \
    --index-batch-size 128G
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
GTDB_PREFIX_TO_RANK = {
    "d": "superkingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
DEFAULT_EXTENSIONS = [
    ".fa", ".fna", ".fasta",
    ".fa.gz", ".fna.gz", ".fasta.gz",
]


def log(message: str) -> None:
    print(message, file=sys.stderr)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def safe_id(value: str) -> str:
    value = str(value).strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "NA"


def strip_fasta_suffix(name: str) -> str:
    for suffix in [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def open_output(path: Path, gzip_output: bool):
    if gzip_output:
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def read_table(path: Path, sep: str = "\t") -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=sep)
        if reader.fieldnames is None:
            fail(f"Taxonomy table has no header: {path}")
        return [dict(row) for row in reader]


def parse_gtdb_classification(classification: str) -> Dict[str, str]:
    parsed = {rank: "" for rank in RANKS}
    if classification is None:
        return parsed
    classification = str(classification).strip()
    if classification.upper() in {"", "NA", "NAN", "NULL"}:
        return parsed
    for part in classification.split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^([dpcofgs])__(.*)$", part)
        if not m:
            continue
        rank = GTDB_PREFIX_TO_RANK[m.group(1)]
        value = m.group(2).strip()
        parsed[rank] = value
    return parsed


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    seq_parts: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(re.sub(r"\s+", "", line))
        if header is not None:
            yield header, "".join(seq_parts)


def write_wrapped(handle, seq: str, width: int = 80) -> None:
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def file_candidates(row: Dict[str, str], genome_id_col: str, filename_col: Optional[str]) -> List[str]:
    candidates: List[str] = []
    if filename_col and filename_col in row and row.get(filename_col):
        candidates.append(row[filename_col])
    if genome_id_col in row and row.get(genome_id_col):
        candidates.append(row[genome_id_col])
    if "genome_path" in row and row.get("genome_path"):
        candidates.append(Path(row["genome_path"]).name)
    expanded: List[str] = []
    for item in candidates:
        base = Path(str(item)).name
        expanded.extend([base, strip_fasta_suffix(base)])
    seen = set()
    out = []
    for x in expanded:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return out


def discover_genomes(genomes_dir: Path, extensions: List[str], recursive: bool) -> Dict[str, Path]:
    paths: Iterable[Path]
    if recursive:
        paths = (p for p in genomes_dir.rglob("*") if p.is_file())
    else:
        paths = (p for p in genomes_dir.iterdir() if p.is_file())

    index: Dict[str, Path] = {}
    for path in paths:
        name = path.name
        if not any(name.endswith(ext) for ext in extensions):
            continue
        keys = {name, strip_fasta_suffix(name)}
        for key in keys:
            if key not in index:
                index[key] = path
    return index


def resolve_genome_file(
    row: Dict[str, str],
    genome_index: Dict[str, Path],
    genome_id_col: str,
    filename_col: Optional[str],
) -> Optional[Path]:
    for candidate in file_candidates(row, genome_id_col, filename_col):
        if candidate in genome_index:
            return genome_index[candidate]
    return None


def make_seq_id(genome_id: str, contig_index: int, original_header: str, mode: str) -> str:
    genome_safe = safe_id(genome_id)
    original_first_token = safe_id(original_header.split()[0] if original_header else f"contig_{contig_index}")
    if mode == "numbered":
        return f"{genome_safe}__ctg{contig_index:06d}"
    if mode == "original":
        return f"{genome_safe}__{original_first_token}"
    digest = hashlib.sha1(original_header.encode("utf-8", errors="replace")).hexdigest()[:10]
    return f"{genome_safe}__ctg{contig_index:06d}__{digest}"


def maybe_build_index(args, fasta_path: Path, index_path: Path) -> None:
    if not args.build_index:
        return
    minimap2 = shutil.which(args.minimap2) if os.sep not in args.minimap2 else args.minimap2
    if not minimap2 or not Path(minimap2).exists():
        fail(f"minimap2 executable not found: {args.minimap2}")
    cmd = [
        minimap2,
        "-d",
        str(index_path),
        "-I",
        args.index_batch_size,
        "-t",
        str(args.threads),
        str(fasta_path),
    ]
    log("Building minimap2 index:")
    log("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a contig-level minitax FASTA/database from MAG FASTAs and a GTDB-Tk-style taxonomy table."
    )
    parser.add_argument("--genomes-dir", required=True, type=Path, help="Folder containing genome/MAG FASTA files.")
    parser.add_argument("--taxonomy-tsv", required=True, type=Path, help="Taxonomy table containing genome IDs and GTDB-Tk classification.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output database directory.")
    parser.add_argument("--prefix", default="MAG_contigs", help="Output FASTA/index prefix. Default: MAG_contigs")
    parser.add_argument("--genome-id-col", default="genome_id", help="Genome ID column. Default: genome_id")
    parser.add_argument("--filename-col", default="file_name", help="Filename column. Use empty string to disable. Default: file_name")
    parser.add_argument("--taxonomy-col", default="classification", help="GTDB-Tk classification column. Default: classification")
    parser.add_argument("--recursive", action="store_true", help="Search genomes-dir recursively.")
    parser.add_argument("--extensions", default=",".join(DEFAULT_EXTENSIONS), help="Comma-separated FASTA extensions to include.")
    parser.add_argument("--min-contig-length", type=int, default=0, help="Discard contigs shorter than this length. Default: 0")
    parser.add_argument("--header-mode", choices=["numbered", "original", "hash"], default="numbered", help="How to name rewritten contigs. Default: numbered")
    parser.add_argument("--gzip", action="store_true", help="Compress output FASTA with gzip.")
    parser.add_argument("--force", action="store_true", help="Overwrite output files if present.")
    parser.add_argument("--build-index", action="store_true", help="Also build a minimap2 index from the output FASTA.")
    parser.add_argument("--minimap2", default="minimap2", help="minimap2/mm2-fast executable. Default: minimap2")
    parser.add_argument("--threads", type=int, default=4, help="Threads for minimap2 index build. Default: 4")
    parser.add_argument("--index-batch-size", default="128G", help="minimap2 -I value. Default: 128G")
    args = parser.parse_args()

    if args.min_contig_length < 0:
        fail("--min-contig-length must be >= 0")
    if args.threads < 1:
        fail("--threads must be >= 1")
    if not args.genomes_dir.is_dir():
        fail(f"Genome folder not found: {args.genomes_dir}")
    if not args.taxonomy_tsv.is_file():
        fail(f"Taxonomy table not found: {args.taxonomy_tsv}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = args.outdir / f"{args.prefix}.fna"
    if args.gzip:
        fasta_path = fasta_path.with_suffix(fasta_path.suffix + ".gz")
    index_path = args.outdir / f"{args.prefix}.idx"

    output_paths = [
        fasta_path,
        args.outdir / "db_data.tsv",
        args.outdir / "MAG.db.tsv",
        args.outdir / "MAG.db.uni.tsv",
        args.outdir / "contig_manifest.tsv",
        args.outdir / "build_manifest.tsv",
        args.outdir / "unmatched_genomes.tsv",
        args.outdir / "missing_files.tsv",
    ]
    if not args.force:
        existing = [str(p) for p in output_paths if p.exists()]
        if existing:
            fail("Output file(s) already exist; use --force to overwrite: " + ", ".join(existing))

    filename_col = args.filename_col if args.filename_col else None
    extensions = [x.strip() for x in args.extensions.split(",") if x.strip()]
    rows = read_table(args.taxonomy_tsv)
    if not rows:
        fail("Taxonomy table has no rows.")
    fieldnames = set(rows[0].keys())
    for col in [args.genome_id_col, args.taxonomy_col]:
        if col not in fieldnames:
            fail(f"Required column not found in taxonomy table: {col}")
    if filename_col and filename_col not in fieldnames:
        log(f"WARNING: filename column '{filename_col}' not found; matching by genome_id/genome_path only.")
        filename_col = None

    genome_index = discover_genomes(args.genomes_dir, extensions, args.recursive)
    log(f"Discovered {len(set(genome_index.values()))} genome FASTA file(s) in {args.genomes_dir}")

    db_rows: List[Dict[str, str]] = []
    manifest_rows: List[Dict[str, str]] = []
    missing_rows: List[Dict[str, str]] = []
    used_files = set()
    seen_seqnames = set()
    total_bases = 0
    total_contigs = 0
    kept_genomes = 0

    with open_output(fasta_path, args.gzip) as fasta_out:
        for row in rows:
            genome_id = str(row.get(args.genome_id_col, "")).strip()
            if not genome_id:
                continue
            genome_path = resolve_genome_file(row, genome_index, args.genome_id_col, filename_col)
            if genome_path is None:
                missing_rows.append({
                    "genome_id": genome_id,
                    "file_name": row.get(filename_col, "") if filename_col else "",
                    "searched_keys": ";".join(file_candidates(row, args.genome_id_col, filename_col)),
                })
                continue

            taxonomy = parse_gtdb_classification(row.get(args.taxonomy_col, ""))
            contig_count = 0
            kept_contig_count = 0
            genome_bases = 0
            used_files.add(genome_path)

            for original_header, seq in iter_fasta(genome_path):
                contig_count += 1
                seq = seq.upper()
                if len(seq) < args.min_contig_length:
                    continue
                kept_contig_count += 1
                seqname = make_seq_id(genome_id, kept_contig_count, original_header, args.header_mode)
                if seqname in seen_seqnames:
                    fail(f"Duplicate generated seqname: {seqname}. Try --header-mode hash.")
                seen_seqnames.add(seqname)

                fasta_out.write(f">{seqname} genome_id={safe_id(genome_id)} original_header={safe_id(original_header.split()[0] if original_header else '')}\n")
                write_wrapped(fasta_out, seq)

                db_row = {
                    "seqnames": seqname,
                    "taxid": safe_id(genome_id),
                    **taxonomy,
                }
                db_rows.append(db_row)
                manifest_rows.append({
                    "seqnames": seqname,
                    "genome_id": genome_id,
                    "genome_file": str(genome_path),
                    "original_header": original_header,
                    "length": str(len(seq)),
                    "taxid": safe_id(genome_id),
                    **taxonomy,
                })
                total_bases += len(seq)
                genome_bases += len(seq)
                total_contigs += 1

            if kept_contig_count > 0:
                kept_genomes += 1
            else:
                log(f"WARNING: no contigs retained for genome {genome_id}: {genome_path}")

    db_fields = ["seqnames", "taxid", *RANKS]
    with open(args.outdir / "db_data.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=db_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(db_rows)
    shutil.copyfile(args.outdir / "db_data.tsv", args.outdir / "MAG.db.tsv")
    shutil.copyfile(args.outdir / "db_data.tsv", args.outdir / "MAG.db.uni.tsv")

    manifest_fields = ["seqnames", "genome_id", "genome_file", "original_header", "length", "taxid", *RANKS]
    with open(args.outdir / "contig_manifest.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(manifest_rows)

    with open(args.outdir / "missing_files.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome_id", "file_name", "searched_keys"], delimiter="\t")
        writer.writeheader()
        writer.writerows(missing_rows)

    unmatched_paths = sorted(set(genome_index.values()) - used_files)
    with open(args.outdir / "unmatched_genomes.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome_file"], delimiter="\t")
        writer.writeheader()
        for p in unmatched_paths:
            writer.writerow({"genome_file": str(p)})

    with open(args.outdir / "build_manifest.tsv", "w", encoding="utf-8") as handle:
        handle.write("key\tvalue\n")
        handle.write(f"taxonomy_tsv\t{args.taxonomy_tsv}\n")
        handle.write(f"genomes_dir\t{args.genomes_dir}\n")
        handle.write(f"genome_id_col\t{args.genome_id_col}\n")
        handle.write(f"filename_col\t{filename_col or ''}\n")
        handle.write(f"taxonomy_col\t{args.taxonomy_col}\n")
        handle.write(f"fasta\t{fasta_path}\n")
        handle.write(f"db_data\t{args.outdir / 'db_data.tsv'}\n")
        handle.write(f"genomes_in_taxonomy\t{len(rows)}\n")
        handle.write(f"genomes_with_files\t{kept_genomes}\n")
        handle.write(f"missing_genome_files\t{len(missing_rows)}\n")
        handle.write(f"contigs_written\t{total_contigs}\n")
        handle.write(f"total_bases\t{total_bases}\n")
        handle.write(f"min_contig_length\t{args.min_contig_length}\n")
        handle.write(f"header_mode\t{args.header_mode}\n")
        handle.write(f"index\t{index_path}\n")
        handle.write(f"index_batch_size\t{args.index_batch_size}\n")

    maybe_build_index(args, fasta_path, index_path)

    log("Done.")
    log(f"  FASTA: {fasta_path}")
    log(f"  db_data.tsv: {args.outdir / 'db_data.tsv'}")
    log(f"  contigs written: {total_contigs}")
    log(f"  total bases: {total_bases}")
    log(f"  genomes with files: {kept_genomes}")
    log(f"  missing genome files: {len(missing_rows)}")
    if missing_rows:
        log(f"  inspect missing files: {args.outdir / 'missing_files.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
