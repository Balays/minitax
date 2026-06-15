#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build an NCBI RefSeq TargetedLoci 16S database for minitax.

This creates a minitax-compatible 16S reference bundle from NCBI RefSeq
TargetedLoci Bacteria/Archaea files:
  ncbi_refseq_16s.fna.gz
  ncbi_refseq_16s.idx
  NCBI.db.tsv
  NCBI.db.uni.tsv
  NCBI.db.uni.spec.tsv
  NCBI_16S_seqlengths.txt
  metadata/ncbi_refseq_16s_accession_taxid.tsv

Usage:
  scripts/build_ncbi_refseq_16s_db.sh --outdir DIR [options]

Options:
  --outdir DIR              Output database directory. Required.
  --kingdoms LIST           TargetedLoci kingdoms to download. Default: Bacteria,Archaea
  --min-length N            Keep only FASTA records >= N bp. Default: 0.
                            For full-length V1-V9, 1200 is often useful.
  --threads N               Threads for pigz/minimap2 indexing. Default: 4
  --index-batch-size SIZE   minimap2 -I batch size. Default: 128G
  --rscript PATH            Rscript path. Default: command -v Rscript
  --minimap2 PATH           minimap2 path. Default: command -v minimap2
  --skip-index              Do not build the minimap2 index.
  --force                   Rebuild combined FASTA and tables from downloaded files.
  -h, --help                Show this help.

Examples:
  scripts/build_ncbi_refseq_16s_db.sh \
    --outdir /mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s \
    --min-length 1200 \
    --threads 20

Use the resulting bundle in minitax_config.txt as:
  db               all_NCBI_genomes
  db.dir           /path/to/ncbi_targeted_loci_16s
  mm2_ref          ncbi_refseq_16s.fna.gz
  mm2_index        ncbi_refseq_16s.idx
EOF
}

fail() { echo "ERROR: $*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

download_file() {
  local url="$1"
  local dest="$2"
  if [[ -s "$dest" ]]; then
    echo "Using existing: $dest"
    return 0
  fi
  echo "Downloading: $url"
  if have curl; then
    curl -L --fail --retry 5 --retry-delay 5 -o "$dest" "$url"
  elif have wget; then
    wget -O "$dest" "$url"
  else
    fail "Need curl or wget for downloads."
  fi
}

csv_to_lines() {
  printf '%s' "$1" | awk -v RS=',' '{gsub(/^[ \t]+|[ \t]+$/, ""); if ($0 != "") print}'
}

OUTDIR=""
KINGDOMS="Bacteria,Archaea"
MIN_LENGTH=0
THREADS=4
INDEX_BATCH_SIZE="128G"
RSCRIPT_BIN="${RSCRIPT:-}"
MINIMAP2_BIN="${MINIMAP2:-}"
SKIP_INDEX=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --kingdoms) KINGDOMS="${2:-}"; shift 2 ;;
    --min-length) MIN_LENGTH="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --index-batch-size) INDEX_BATCH_SIZE="${2:-}"; shift 2 ;;
    --rscript) RSCRIPT_BIN="${2:-}"; shift 2 ;;
    --minimap2) MINIMAP2_BIN="${2:-}"; shift 2 ;;
    --skip-index) SKIP_INDEX=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$OUTDIR" ]] || { usage; fail "--outdir is required."; }
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || fail "--threads must be a positive integer."
[[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] || fail "--min-length must be a non-negative integer."
[[ "$INDEX_BATCH_SIZE" =~ ^[0-9]+[KkMmGgTt]?$ ]] || fail "--index-batch-size must look like 128G, 64000M, or 500000000."

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAKE_DB_R="${SCRIPT_DIR}/make_ncbi_refseq_16s_db.R"
[[ -f "$MAKE_DB_R" ]] || fail "Cannot find helper script: $MAKE_DB_R"

if [[ -z "$RSCRIPT_BIN" ]]; then
  RSCRIPT_BIN="$(command -v Rscript || true)"
fi
[[ -n "$RSCRIPT_BIN" && -x "$RSCRIPT_BIN" ]] || fail "Rscript not found. Activate an R environment or pass --rscript PATH."

if [[ "$SKIP_INDEX" -eq 0 && -z "$MINIMAP2_BIN" ]]; then
  MINIMAP2_BIN="$(command -v minimap2 || true)"
elif [[ "$SKIP_INDEX" -eq 0 && "$MINIMAP2_BIN" != */* ]]; then
  MINIMAP2_BIN="$(command -v "$MINIMAP2_BIN" || true)"
fi
if [[ "$SKIP_INDEX" -eq 0 ]]; then
  [[ -n "$MINIMAP2_BIN" && -x "$MINIMAP2_BIN" ]] || fail "minimap2 not found. Activate an environment or pass --minimap2 PATH."
fi

mkdir -p "$OUTDIR"/{downloads,metadata,logs,taxdump}

RAW_FASTA="${OUTDIR}/ncbi_refseq_16s.raw.fna"
FILTERED_FASTA="${OUTDIR}/ncbi_refseq_16s.fna"
FASTA_GZ="${FILTERED_FASTA}.gz"
SEQLENGTHS="${OUTDIR}/NCBI_16S_seqlengths.txt"
TAXDUMP="${OUTDIR}/taxdump/taxdump.tar.gz"
ACCESSION_TAXID="${OUTDIR}/metadata/ncbi_refseq_16s_accession_taxid.tsv"
GBFF_LIST="${OUTDIR}/metadata/gbff_files.txt"

if [[ "$FORCE" -eq 1 ]]; then
  rm -f "$RAW_FASTA" "$FILTERED_FASTA" "$FASTA_GZ" "$SEQLENGTHS" "$ACCESSION_TAXID" "$GBFF_LIST" \
        "${OUTDIR}/NCBI.db.tsv" "${OUTDIR}/NCBI.db.uni.tsv" "${OUTDIR}/NCBI.db.uni.spec.tsv" \
        "${OUTDIR}/ncbi_refseq_16s.idx"
fi

: > "$RAW_FASTA"
: > "$GBFF_LIST"

while read -r kingdom; do
  [[ -n "$kingdom" ]] || continue
  kingdom_lc="$(printf '%s' "$kingdom" | tr '[:upper:]' '[:lower:]')"
  fna="${kingdom_lc}.16SrRNA.fna.gz"
  gbff="${kingdom_lc}.16SrRNA.gbff.gz"
  base_url="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/${kingdom}"
  download_file "${base_url}/${fna}" "${OUTDIR}/downloads/${fna}"
  download_file "${base_url}/${gbff}" "${OUTDIR}/downloads/${gbff}"
  echo "${OUTDIR}/downloads/${gbff}" >> "$GBFF_LIST"
  gzip -cd "${OUTDIR}/downloads/${fna}" >> "$RAW_FASTA"
done < <(csv_to_lines "$KINGDOMS")

# Normalize FASTA headers to first token only, and optionally length-filter.
awk -v min_len="$MIN_LENGTH" -v seq_lengths="$SEQLENGTHS" '
  function flush_record() {
    if (name != "") {
      if (len >= min_len) {
        print ">" name
        for (i = 1; i <= nseq; i++) print seq[i]
        print name "\t" len >> seq_lengths
      }
    }
    delete seq
    nseq = 0
    len = 0
  }
  BEGIN { OFS = "\t"; print "seqnames\tseqlengths" > seq_lengths }
  /^>/ {
    flush_record()
    header = $0
    sub(/^>/, "", header)
    split(header, fields, /[ \t]/)
    name = fields[1]
    next
  }
  {
    gsub(/\r/, "")
    if ($0 != "") {
      nseq++
      seq[nseq] = $0
      len += length($0)
    }
  }
  END { flush_record() }
' "$RAW_FASTA" > "$FILTERED_FASTA"

if have pigz; then
  pigz -f -p "$THREADS" "$FILTERED_FASTA"
else
  gzip -f "$FILTERED_FASTA"
fi

# Taxdump is always downloaded because it is small enough and ensures current lineage names.
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "$TAXDUMP"
tar -xzf "$TAXDUMP" -C "${OUTDIR}/taxdump" names.dmp nodes.dmp merged.dmp

"$RSCRIPT_BIN" "$MAKE_DB_R" \
  --gbff-list "$GBFF_LIST" \
  --fasta "$FASTA_GZ" \
  --seq-lengths "$SEQLENGTHS" \
  --taxdump-dir "${OUTDIR}/taxdump" \
  --outdir "$OUTDIR"

cat > "${OUTDIR}/build_manifest.tsv" <<EOF
key	value
database	NCBI RefSeq TargetedLoci 16S
kingdoms	$KINGDOMS
min_length	$MIN_LENGTH
fasta	$FASTA_GZ
index	${OUTDIR}/ncbi_refseq_16s.idx
index_batch_size	$INDEX_BATCH_SIZE
EOF

if [[ "$SKIP_INDEX" -eq 0 ]]; then
  echo "Building minimap2 index with -I ${INDEX_BATCH_SIZE}..."
  "${SCRIPT_DIR}/build_minimap2_index.sh" \
    --db-dir "$OUTDIR" \
    --fasta "$FASTA_GZ" \
    --index ncbi_refseq_16s.idx \
    --threads "$THREADS" \
    --index-batch-size "$INDEX_BATCH_SIZE" \
    --minimap2 "$MINIMAP2_BIN" \
    --force
else
  echo "Skipping minimap2 index."
fi

echo "Done. Use this in minitax_config.txt:"
echo -e "  db\tall_NCBI_genomes"
echo -e "  db.dir\t$OUTDIR"
echo -e "  mm2_ref\tncbi_refseq_16s.fna.gz"
echo -e "  mm2_index\tncbi_refseq_16s.idx"
