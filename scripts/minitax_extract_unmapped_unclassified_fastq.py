#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Set, TextIO

PRIMARY_EXCLUDE = 0x100 | 0x800
UNMAPPED = 0x4
NA_VALUES = {'', 'NA', 'NAN', 'NULL', 'UNCLASSIFIED', 'UNASSIGNED'}


def log(*x):
    print(*x, file=sys.stderr, flush=True)


def die(msg):
    raise SystemExit('ERROR: ' + msg)


def sample_from_bam(bam: Path) -> str:
    return bam.name[:-4] if bam.name.endswith('.bam') else bam.stem


def find_bams(outdir: Path) -> List[Path]:
    bamdir = outdir / 'bam'
    if not bamdir.is_dir():
        die(f'BAM directory not found: {bamdir}')
    return sorted(p for p in bamdir.glob('*.bam') if not p.name.endswith('.bai'))


def is_na(v) -> bool:
    return str(v or '').strip().upper() in NA_VALUES


def load_unclassified_qnames(cache_tsv: Path) -> Set[str]:
    if not cache_tsv.is_file():
        log('WARNING: cache not found, exporting only unmapped reads:', cache_tsv)
        return set()
    qnames = set()
    with open(cache_tsv, encoding='utf-8', errors='replace', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames or 'qname' not in reader.fieldnames:
            log('WARNING: cache has no qname column:', cache_tsv)
            return set()
        has_tax = 'tax.identity' in reader.fieldnames
        has_level = 'tax.identity.level' in reader.fieldnames
        rank_cols = [c for c in ['superkingdom','phylum','class','order','family','genus','species'] if c in reader.fieldnames]
        for row in reader:
            q = str(row.get('qname') or '').strip()
            if not q:
                continue
            unclassified = False
            if has_tax and is_na(row.get('tax.identity')):
                unclassified = True
            if has_level and is_na(row.get('tax.identity.level')):
                unclassified = True
            if (not has_tax) and rank_cols and all(is_na(row.get(c)) for c in rank_cols):
                unclassified = True
            if unclassified:
                qnames.add(q)
    return qnames


def open_out(path: Path, gz: bool) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        return gzip.open(path, 'wt', encoding='utf-8', newline='')
    return open(path, 'w', encoding='utf-8', newline='')


def write_fastq(out: TextIO, qname: str, seq: str, qual: str, source: str) -> bool:
    if not seq or seq == '*':
        return False
    if not qual or qual == '*' or len(qual) != len(seq):
        qual = 'I' * len(seq)
    out.write(f'@{qname} source={source}\n{seq}\n+\n{qual}\n')
    return True


def export_one(bam: Path, outdir: Path, samtools: str, samtools_threads: int, gz: bool, force: bool, quiet: bool) -> Dict[str, object]:
    sample = sample_from_bam(bam)
    cache = outdir / 'best_alignments_w_taxa' / f'{sample}_best_alignments_w_taxa.tsv'
    out_fastq = outdir / 'unmapped_unclassified_fastq' / (f'{sample}_unmapped_unclassified.fastq.gz' if gz else f'{sample}_unmapped_unclassified.fastq')
    if out_fastq.exists() and not force:
        die(f'Output exists; use --force to overwrite: {out_fastq}')

    t0 = time.time()
    if not quiet:
        log('START', sample)
    unclassified = load_unclassified_qnames(cache)
    if not quiet:
        log(sample, 'unclassified qnames from cache:', len(unclassified))

    cmd = [samtools, 'view']
    if samtools_threads > 0:
        cmd += ['-@', str(samtools_threads)]
    cmd.append(str(bam))

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='replace')
    assert proc.stdout is not None
    seen = set()
    primary_seen = 0
    unmapped_written = 0
    unclassified_written = 0
    total_written = 0

    with open_out(out_fastq, gz) as out:
        for line in proc.stdout:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                continue
            qname = fields[0]
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            if flag & PRIMARY_EXCLUDE:
                continue
            primary_seen += 1
            is_unmapped = bool(flag & UNMAPPED)
            is_unclassified = qname in unclassified
            if not (is_unmapped or is_unclassified) or qname in seen:
                continue
            source = 'unmapped' if is_unmapped else 'minitax_unclassified'
            if write_fastq(out, qname, fields[9], fields[10], source):
                seen.add(qname)
                total_written += 1
                if is_unmapped:
                    unmapped_written += 1
                else:
                    unclassified_written += 1
    stderr = proc.stderr.read() if proc.stderr else ''
    rc = proc.wait()
    if rc != 0:
        die(f'samtools view failed for {bam}: {stderr}')
    elapsed = time.time() - t0
    if not quiet:
        log('FINISH', sample, 'primary_seen=', primary_seen, 'unmapped=', unmapped_written, 'unclassified=', unclassified_written, 'total=', total_written)
    return {
        'sample': sample,
        'bam': str(bam),
        'cache_tsv': str(cache),
        'output_fastq': str(out_fastq),
        'primary_records_seen': primary_seen,
        'unclassified_qnames_from_cache': len(unclassified),
        'unmapped_reads_written': unmapped_written,
        'unclassified_reads_written': unclassified_written,
        'total_reads_written': total_written,
        'elapsed_sec': f'{elapsed:.3f}',
    }


def write_manifest(path: Path, rows: List[Dict[str, object]]):
    fields = ['sample','bam','cache_tsv','output_fastq','primary_records_seen','unclassified_qnames_from_cache','unmapped_reads_written','unclassified_reads_written','total_reads_written','elapsed_sec']
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--samtools', default='samtools')
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--samtools-threads', type=int, default=0)
    ap.add_argument('--gzip', action='store_true')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--manifest', type=Path, default=None)
    args = ap.parse_args()
    if args.threads < 1:
        die('--threads must be >= 1')
    bams = find_bams(args.outdir)
    if not bams:
        die('No BAM files found')
    workers = min(args.threads, len(bams))
    manifest = args.manifest or args.outdir / 'unmapped_unclassified_fastq' / 'manifest.tsv'
    log('Found', len(bams), 'BAM(s); workers:', workers)
    rows_by_sample = {}
    if workers == 1:
        for bam in bams:
            row = export_one(bam, args.outdir, args.samtools, args.samtools_threads, args.gzip, args.force, args.quiet)
            rows_by_sample[row['sample']] = row
    else:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(export_one, bam, args.outdir, args.samtools, args.samtools_threads, args.gzip, args.force, args.quiet): bam for bam in bams}
            done = 0
            for fut in as_completed(futs):
                row = fut.result()
                rows_by_sample[row['sample']] = row
                done += 1
                log(f'[{done}/{len(bams)}] completed', futs[fut].name)
    rows = [rows_by_sample[sample_from_bam(b)] for b in bams]
    write_manifest(manifest, rows)
    log('Wrote:', manifest)


if __name__ == '__main__':
    main()
