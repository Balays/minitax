#!/usr/bin/env python3
"""
Create per-sample read accounting tables for minitax outputs.

This script compares BAM primary-read counts with minitax final taxonomic summaries.
It makes the denominator explicit: total primary reads, primary mapped reads,
primary unmapped reads, final assigned reads, and mapped-but-unassigned/filtered reads.

It uses primary alignments as the read-level denominator:
  primary_total    = samtools view -c -F 0x900 BAM
  primary_mapped   = samtools view -c -F 0x904 BAM
  primary_unmapped = samtools view -c -f 0x4 -F 0x900 BAM

The total BAM record count is also reported, but should not be used as a read-count denominator
for long-read/minimap2 BAMs because secondary and supplementary alignments inflate it.

Example:
  python3 scripts/minitax_read_accounting.py \
    --outdir minitax_CanMAG \
    --method BestAln \
    --threads 20 \
    --write-normalized
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def eprint(*args) -> None:
    print(*args, file=sys.stderr, flush=True)


def fail(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def run_count(cmd: List[str], label: str = "") -> int:
    if label:
        eprint(f"    running: {label}")
    try:
        out = subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE).strip()
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else ""
        fail(f"Command failed: {' '.join(cmd)}\n{stderr}")
    try:
        value = int(out)
    except ValueError:
        fail(f"Could not parse integer from command output: {' '.join(cmd)} -> {out!r}")
    if label:
        eprint(f"    done: {label} = {value}")
    return value


def sample_from_bam(path: Path) -> str:
    name = path.name
    if name.endswith(".bam"):
        return name[:-4]
    return path.stem


def find_bams(outdir: Path) -> List[Path]:
    bamdir = outdir / "bam"
    if not bamdir.is_dir():
        fail(f"BAM directory not found: {bamdir}")
    return sorted(p for p in bamdir.glob("*.bam") if not p.name.endswith(".bai"))


def find_method_table(outdir: Path, method: str, explicit: Optional[Path] = None) -> Path:
    if explicit is not None:
        if not explicit.is_file():
            fail(f"Method table not found: {explicit}")
        return explicit
    patterns = [
        f"*_{method}_taxa.all.sum.tsv",
        f"{method}_taxa.all.sum.tsv",
    ]
    hits: List[Path] = []
    for pattern in patterns:
        hits.extend(outdir.glob(pattern))
    hits = sorted(set(hits))
    if len(hits) == 0:
        fail(f"Could not find method summary table in {outdir} using patterns: {', '.join(patterns)}")
    if len(hits) > 1:
        eprint("WARNING: multiple method summary tables found; using:", hits[0])
        for hit in hits[1:]:
            eprint("  ignored:", hit)
    return hits[0]


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"TSV has no header: {path}")
        rows = [dict(row) for row in reader]
    return list(reader.fieldnames), rows


def numeric(value: str) -> float:
    if value is None:
        return 0.0
    value = str(value).strip()
    if value == "" or value.upper() in {"NA", "NAN", "NULL"}:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def assigned_counts_by_sample(method_table: Path, method: str) -> Dict[str, float]:
    header, rows = read_tsv(method_table)
    if "sample" not in header:
        fail(f"Method table is missing 'sample' column: {method_table}")
    if method == "SpeciesEstimate" and "norm_count" in header:
        count_col = "norm_count"
    elif "count" in header:
        count_col = "count"
    elif "norm_count" in header:
        count_col = "norm_count"
    else:
        fail(f"Method table has neither 'count' nor 'norm_count': {method_table}")

    counts: Dict[str, float] = {}
    for row in rows:
        sample = row.get("sample", "")
        if not sample:
            continue
        counts[sample] = counts.get(sample, 0.0) + numeric(row.get(count_col, "0"))
    return counts


def fraction(num: float, denom: float) -> str:
    if denom == 0 or math.isnan(denom):
        return "NA"
    return f"{num / denom:.10g}"


def samtools_view_count_args(samtools: str, samtools_threads: int, extra_args: List[str], bam: Path) -> List[str]:
    cmd = [samtools, "view"]
    if samtools_threads > 0:
        cmd.extend(["-@", str(samtools_threads)])
    cmd.extend(["-c"])
    cmd.extend(extra_args)
    cmd.append(str(bam))
    return cmd


def bam_counts(bam: Path, samtools: str, samtools_threads: int = 0, verbose: bool = True) -> Dict[str, int]:
    sample = sample_from_bam(bam)

    def count(metric: str, extra_args: List[str]) -> int:
        label = f"{sample}: {metric}" if verbose else ""
        return run_count(samtools_view_count_args(samtools, samtools_threads, extra_args, bam), label=label)

    return {
        "bam_records_total": count("all BAM records", []),
        "primary_total_reads": count("primary total reads (-F 0x900)", ["-F", "0x900"]),
        "primary_mapped_reads": count("primary mapped reads (-F 0x904)", ["-F", "0x904"]),
        "primary_unmapped_reads": count("primary unmapped reads (-f 0x4 -F 0x900)", ["-f", "0x4", "-F", "0x900"]),
        "secondary_records": count("secondary records (-f 0x100)", ["-f", "0x100"]),
        "supplementary_records": count("supplementary records (-f 0x800)", ["-f", "0x800"]),
    }


def build_accounting_row(
    bam: Path,
    method: str,
    method_table: Path,
    assigned: Dict[str, float],
    samtools: str,
    samtools_threads: int,
    verbose: bool,
) -> Dict[str, object]:
    sample = sample_from_bam(bam)
    started = time.time()
    if verbose:
        eprint(f"START {sample}: {bam}")
    counts = bam_counts(bam, samtools, samtools_threads=samtools_threads, verbose=verbose)
    final_assigned = assigned.get(sample, 0.0)
    primary_mapped = counts["primary_mapped_reads"]
    primary_total = counts["primary_total_reads"]
    mapped_unassigned = max(0.0, float(primary_mapped) - float(final_assigned))
    elapsed = time.time() - started
    if verbose:
        eprint(
            f"FINISH {sample}: primary_total={primary_total}, primary_mapped={primary_mapped}, "
            f"final_assigned={final_assigned:.10g}, elapsed={elapsed:.1f}s"
        )
    return {
        "sample": sample,
        "method": method,
        "bam": str(bam),
        "method_table": str(method_table),
        **counts,
        "final_assigned_reads": f"{final_assigned:.10g}",
        "mapped_unassigned_or_filtered_reads": f"{mapped_unassigned:.10g}",
        "mapped_fraction_of_primary": fraction(float(primary_mapped), float(primary_total)),
        "unmapped_fraction_of_primary": fraction(float(counts["primary_unmapped_reads"]), float(primary_total)),
        "assigned_fraction_of_primary": fraction(float(final_assigned), float(primary_total)),
        "assigned_fraction_of_primary_mapped": fraction(float(final_assigned), float(primary_mapped)),
    }


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def add_normalized_rows(method_table: Path, accounting_rows: List[Dict[str, object]], method: str, out_path: Path) -> None:
    header, rows = read_tsv(method_table)
    if "sample" not in header:
        fail(f"Method table is missing 'sample' column: {method_table}")
    count_col = "norm_count" if method == "SpeciesEstimate" and "norm_count" in header else "count"
    if count_col not in header:
        fail(f"Cannot add normalized rows because {count_col!r} is missing from {method_table}")

    rank_cols = [c for c in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"] if c in header]
    extra_cols = ["tax.identity", "tax.identity.level", "lineage"]
    for col in extra_cols:
        if col not in header:
            header.append(col)

    acc_by_sample = {str(r["sample"]): r for r in accounting_rows}
    out_rows = list(rows)
    for sample, acc in acc_by_sample.items():
        for label, value_key in [
            ("unmapped", "primary_unmapped_reads"),
            ("mapped_unassigned_or_filtered", "mapped_unassigned_or_filtered_reads"),
        ]:
            value = float(acc.get(value_key, 0) or 0)
            if value <= 0:
                continue
            row = {col: "" for col in header}
            row["sample"] = sample
            row[count_col] = f"{value:.10g}"
            for rank in rank_cols:
                row[rank] = label
            row["tax.identity"] = label
            row["tax.identity.level"] = label
            row["lineage"] = label
            out_rows.append(row)

    write_tsv(out_path, out_rows, header)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write read_accounting.tsv for a minitax output directory.")
    parser.add_argument("--outdir", required=True, type=Path, help="minitax output directory containing bam/")
    parser.add_argument("--method", default="BestAln", help="Method to account against, e.g. BestAln, LCA, SpeciesEstimate. Default: BestAln")
    parser.add_argument("--method-table", type=Path, default=None, help="Explicit method taxa.all.sum.tsv table. Default: auto-detect in outdir")
    parser.add_argument("--samtools", default="samtools", help="samtools executable. Default: samtools")
    parser.add_argument("--threads", type=int, default=1, help="Number of BAM files to process in parallel. Default: 1")
    parser.add_argument("--samtools-threads", type=int, default=0, help="Threads passed to each samtools view call using -@. Default: 0. Use carefully together with --threads.")
    parser.add_argument("--quiet", action="store_true", help="Suppress per-command progress messages.")
    parser.add_argument("--output", type=Path, default=None, help="Output accounting TSV. Default: <outdir>/read_accounting.tsv")
    parser.add_argument("--write-normalized", action="store_true", help="Also write a taxa table with synthetic unmapped and mapped_unassigned_or_filtered rows.")
    args = parser.parse_args()

    if args.threads < 1:
        fail("--threads must be >= 1")
    if args.samtools_threads < 0:
        fail("--samtools-threads must be >= 0")
    if not args.outdir.is_dir():
        fail(f"Output directory not found: {args.outdir}")

    method_table = find_method_table(args.outdir, args.method, args.method_table)
    assigned = assigned_counts_by_sample(method_table, args.method)
    bams = find_bams(args.outdir)
    if len(bams) == 0:
        fail(f"No BAM files found under {args.outdir / 'bam'}")

    verbose = not args.quiet
    workers = min(args.threads, len(bams))
    eprint(f"Method table: {method_table}")
    eprint(f"Found {len(assigned)} sample(s) in method table")
    eprint(f"Found {len(bams)} BAM file(s) in {args.outdir / 'bam'}")
    eprint(f"Processing with {workers} BAM worker(s)")
    if args.samtools_threads > 0:
        eprint(f"Each samtools view call will use -@ {args.samtools_threads}")
    eprint("Counting plan per BAM: all records; primary total; primary mapped; primary unmapped; secondary; supplementary")

    rows_by_sample: Dict[str, Dict[str, object]] = {}
    if workers == 1:
        for i, bam in enumerate(bams, start=1):
            eprint(f"[{i}/{len(bams)}] Counting {bam.name}")
            row = build_accounting_row(bam, args.method, method_table, assigned, args.samtools, args.samtools_threads, verbose)
            rows_by_sample[str(row["sample"])] = row
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_bam = {
                executor.submit(build_accounting_row, bam, args.method, method_table, assigned, args.samtools, args.samtools_threads, verbose): bam
                for bam in bams
            }
            completed = 0
            for future in as_completed(future_to_bam):
                bam = future_to_bam[future]
                completed += 1
                try:
                    row = future.result()
                except Exception as exc:
                    fail(f"Failed while processing {bam}: {exc}")
                rows_by_sample[str(row["sample"])] = row
                eprint(f"[{completed}/{len(bams)}] Completed {bam.name}")

    rows = [rows_by_sample[sample_from_bam(bam)] for bam in bams]

    output = args.output or (args.outdir / "read_accounting.tsv")
    fieldnames = [
        "sample", "method", "bam", "method_table",
        "primary_total_reads", "primary_mapped_reads", "primary_unmapped_reads",
        "final_assigned_reads", "mapped_unassigned_or_filtered_reads",
        "mapped_fraction_of_primary", "unmapped_fraction_of_primary",
        "assigned_fraction_of_primary", "assigned_fraction_of_primary_mapped",
        "bam_records_total", "secondary_records", "supplementary_records",
    ]
    write_tsv(output, rows, fieldnames)
    eprint(f"Wrote: {output}")

    if args.write_normalized:
        normalized = args.outdir / f"{args.method}_taxa.all.sum.input_normalized.tsv"
        add_normalized_rows(method_table, rows, args.method, normalized)
        eprint(f"Wrote: {normalized}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
