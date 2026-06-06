#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build an NCBI genome database for minitax.

This creates:
  all_NCBI_genomes.fna.gz
  all_NCBI_genomes.idx
  NCBI.db.tsv
  NCBI.db.uni.tsv
  NCBI.db.uni.spec.tsv
  NCBI.db.genomesize.tsv
  NCBI_genome_collection_seqlengths.txt

Usage:
  scripts/build_ncbi_minitax_db.sh --outdir DIR [options]

Options:
  --outdir DIR                Output database directory. Required.
  --sources LIST              refseq, genbank, or both. Default: refseq
  --tax-groups LIST           NCBI genome groups. Default: bacteria,archaea,viral,fungi,protozoa
  --assembly-levels LIST      Assembly levels. Default: Complete Genome,Chromosome,Scaffold
  --refseq-categories LIST    RefSeq categories. Default: reference genome,representative genome
                              Viral records with category "na" are also kept by default.
                              Use "any" to keep all categories for every group.
  --assembly-accessions FILE  Optional GCF_/GCA_ accession allow-list, one per line.
  --max-genomes N             Keep only the first N selected assemblies. Useful for testing.
  --threads N                 Threads for minimap2 indexing. Default: 4
  --index-batch-size SIZE     minimap2 -I batch size. Default: 128G.
                              Use a value larger than the reference size to avoid
                              multi-part indexes with unreliable MAPQ.
  --rscript PATH              Rscript path. Default: command -v Rscript
  --minimap2 PATH             minimap2 path. Default: command -v minimap2
  --skip-index                Do not build the minimap2 index.
  --dry-run                   Download/filter metadata only; do not download genomes.
                              Prints selected assembly count and total genome size.
  --force                     Rebuild combined FASTA and metadata from downloaded genome files.
  -h, --help                  Show this help.

Examples:
  scripts/build_ncbi_minitax_db.sh --outdir /mnt/d/data/databases/all_NCBI_genomes \
    --threads 32

  scripts/build_ncbi_minitax_db.sh --outdir /tmp/minitax_ncbi_test --max-genomes 5 --dry-run

Recommended environment:
  mamba activate minitax
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

have() {
  command -v "$1" >/dev/null 2>&1
}

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

csv_to_awk_set() {
  printf '%s' "$1" | awk -v RS=',' '{gsub(/^[ \t]+|[ \t]+$/, ""); if ($0 != "") print}'
}

OUTDIR=""
SOURCES="refseq"
TAX_GROUPS="bacteria,archaea,viral,fungi,protozoa"
ASSEMBLY_LEVELS="Complete Genome,Chromosome,Scaffold"
REFSEQ_CATEGORIES="reference genome,representative genome"
ASSEMBLY_ACCESSIONS=""
MAX_GENOMES=""
THREADS=4
INDEX_BATCH_SIZE="128G"
RSCRIPT_BIN="${RSCRIPT:-}"
MINIMAP2_BIN="${MINIMAP2:-}"
SKIP_INDEX=0
DRY_RUN=0
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --sources) SOURCES="${2:-}"; shift 2 ;;
    --tax-groups) TAX_GROUPS="${2:-}"; shift 2 ;;
    --assembly-levels) ASSEMBLY_LEVELS="${2:-}"; shift 2 ;;
    --refseq-categories) REFSEQ_CATEGORIES="${2:-}"; shift 2 ;;
    --assembly-accessions) ASSEMBLY_ACCESSIONS="${2:-}"; shift 2 ;;
    --max-genomes) MAX_GENOMES="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --index-batch-size) INDEX_BATCH_SIZE="${2:-}"; shift 2 ;;
    --rscript) RSCRIPT_BIN="${2:-}"; shift 2 ;;
    --minimap2) MINIMAP2_BIN="${2:-}"; shift 2 ;;
    --skip-index) SKIP_INDEX=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$OUTDIR" ]] || { usage; fail "--outdir is required."; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || fail "--threads must be a positive integer."
[[ "$THREADS" -gt 0 ]] || fail "--threads must be greater than 0."
[[ "$INDEX_BATCH_SIZE" =~ ^[0-9]+[KkMmGgTt]?$ ]] || fail "--index-batch-size must look like 128G, 64000M, or 500000000."
if [[ -n "$MAX_GENOMES" && ! "$MAX_GENOMES" =~ ^[0-9]+$ ]]; then
  fail "--max-genomes must be a positive integer."
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAKE_DB_R="${SCRIPT_DIR}/make_ncbi_minitax_db.R"
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

mkdir -p "$OUTDIR"/{metadata,downloads,logs,taxdump}

ALL_ASSEMBLIES="${OUTDIR}/metadata/assembly_summary.all.tsv"
SELECTED_ASSEMBLIES="${OUTDIR}/metadata/assembly_summary.filtered.tsv"
SEQ_MAP="${OUTDIR}/metadata/seq_to_assembly.tsv"
SEQLENGTHS="${OUTDIR}/NCBI_genome_collection_seqlengths.txt"
FASTA="${OUTDIR}/all_NCBI_genomes.fna"
FASTA_GZ="${FASTA}.gz"
TAXDUMP="${OUTDIR}/taxdump/taxdump.tar.gz"

echo "Output directory: $OUTDIR"
echo "Sources: $SOURCES"
echo "Tax groups: $TAX_GROUPS"
echo "Assembly levels: $ASSEMBLY_LEVELS"
echo "RefSeq categories: $REFSEQ_CATEGORIES"
echo "Index batch size: $INDEX_BATCH_SIZE"

tmp_levels="${OUTDIR}/metadata/.levels.txt"
tmp_categories="${OUTDIR}/metadata/.categories.txt"
csv_to_awk_set "$ASSEMBLY_LEVELS" > "$tmp_levels"
csv_to_awk_set "$REFSEQ_CATEGORIES" > "$tmp_categories"

{
  printf 'assembly_accession\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tassembly_level\tftp_path\tsource\ttax_group\tgenome_size\n'
  IFS=',' read -r -a source_array <<< "$SOURCES"
  IFS=',' read -r -a group_array <<< "$TAX_GROUPS"
  for source in "${source_array[@]}"; do
    source="$(printf '%s' "$source" | xargs)"
    [[ "$source" == "refseq" || "$source" == "genbank" ]] || fail "Unsupported source: $source"
    for group in "${group_array[@]}"; do
      group="$(printf '%s' "$group" | xargs)"
      [[ -n "$group" ]] || continue
      summary="${OUTDIR}/metadata/assembly_summary_${source}_${group}.txt"
      url="https://ftp.ncbi.nlm.nih.gov/genomes/${source}/${group}/assembly_summary.txt"
      download_file "$url" "$summary" >&2
      awk -F '\t' -v OFS='\t' -v src="$source" -v grp="$group" '
        $0 !~ /^#/ && NF >= 20 {
          print $1, $5, $6, $7, $8, $9, $10, $12, $20, src, grp, $26
        }
      ' "$summary"
    done
  done
} > "$ALL_ASSEMBLIES"

FILTERED_RAW="${OUTDIR}/metadata/assembly_summary.filtered.raw.tsv"
awk -F '\t' -v OFS='\t' \
  -v levels_file="$tmp_levels" \
  -v categories_file="$tmp_categories" \
  -v categories="$REFSEQ_CATEGORIES" '
  BEGIN {
    while ((getline line < levels_file) > 0) level_ok[line] = 1
    while ((getline line < categories_file) > 0) category_ok[line] = 1
    keep_any_category = (tolower(categories) == "any")
  }
  NR == 1 { print; next }
  $9 == "" || $9 == "na" { next }
  !($8 in level_ok) { next }
  !keep_any_category && !($2 in category_ok) && !($11 == "viral" && $2 == "na") { next }
  { print }
' "$ALL_ASSEMBLIES" > "$FILTERED_RAW"

{
  head -n 1 "$FILTERED_RAW"
  tail -n +2 "$FILTERED_RAW" | sort -t $'\t' -u -k1,1
} > "${SELECTED_ASSEMBLIES}.tmp"

if [[ -n "$ASSEMBLY_ACCESSIONS" ]]; then
  [[ -f "$ASSEMBLY_ACCESSIONS" ]] || fail "Assembly accession file not found: $ASSEMBLY_ACCESSIONS"
  awk -F '\t' -v OFS='\t' '
    NR == FNR {
      gsub(/\r/, "", $1)
      if ($1 != "" && $1 !~ /^#/) wanted[$1] = 1
      next
    }
    FNR == 1 || ($1 in wanted)
  ' "$ASSEMBLY_ACCESSIONS" "${SELECTED_ASSEMBLIES}.tmp" > "${SELECTED_ASSEMBLIES}.allow.tmp"
  mv "${SELECTED_ASSEMBLIES}.allow.tmp" "${SELECTED_ASSEMBLIES}.tmp"
fi

if [[ -n "$MAX_GENOMES" && "$MAX_GENOMES" -gt 0 ]]; then
  awk -v max="$MAX_GENOMES" 'NR == 1 || NR <= max + 1' "${SELECTED_ASSEMBLIES}.tmp" > "$SELECTED_ASSEMBLIES"
else
  mv "${SELECTED_ASSEMBLIES}.tmp" "$SELECTED_ASSEMBLIES"
fi
rm -f "${SELECTED_ASSEMBLIES}.tmp" "$FILTERED_RAW" "$tmp_levels" "$tmp_categories"

selected_count=$(( $(wc -l < "$SELECTED_ASSEMBLIES") - 1 ))
selected_bases="$(awk -F '\t' 'NR > 1 && $12 ~ /^[0-9]+$/ {sum += $12} END {printf "%.0f", sum + 0}' "$SELECTED_ASSEMBLIES")"
selected_gbp="$(awk -v bases="$selected_bases" 'BEGIN {printf "%.2f", bases / 1000000000}')"
echo "Selected assemblies: $selected_count"
echo "Selected genome size: ${selected_bases} bp (${selected_gbp} Gbp)"
[[ "$selected_count" -gt 0 ]] || fail "No assemblies selected. Relax filters or check tax groups/source."

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry run complete. Selected assembly table:"
  echo "  $SELECTED_ASSEMBLIES"
  exit 0
fi

if [[ "$FORCE" -eq 1 ]]; then
  rm -f "$SEQ_MAP" "$SEQLENGTHS" "$FASTA" "$FASTA_GZ"
fi

printf 'seqnames\tident\ttaxid\n' > "$SEQ_MAP"
: > "$SEQLENGTHS"
: > "$FASTA"

tail -n +2 "$SELECTED_ASSEMBLIES" | while IFS=$'\t' read -r acc refcat taxid species_taxid organism infra isolate level ftp source group genome_size; do
  [[ -n "$acc" ]] || continue
  base="$(basename "$ftp")"
  genome_url="${ftp}/${base}_genomic.fna.gz"
  genome_file="${OUTDIR}/downloads/${acc}_${base}_genomic.fna.gz"
  download_file "$genome_url" "$genome_file"
  echo "Adding assembly: $acc"
  gzip -cd "$genome_file" | awk -v OFS='\t' -v acc="$acc" -v taxid="$taxid" -v seq_map="$SEQ_MAP" -v seq_lengths="$SEQLENGTHS" '
    /^>/ {
      if (name != "") print name, len >> seq_lengths
      header = $0
      sub(/^>/, "", header)
      split(header, fields, /[ \t]/)
      name = fields[1]
      print name, acc, taxid >> seq_map
      len = 0
      print
      next
    }
    {
      gsub(/\r/, "")
      len += length($0)
      print
    }
    END {
      if (name != "") print name, len >> seq_lengths
    }
  ' >> "$FASTA"
done

if have pigz; then
  pigz -f -p "$THREADS" "$FASTA"
else
  gzip -f "$FASTA"
fi

download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "$TAXDUMP"
tar -xzf "$TAXDUMP" -C "${OUTDIR}/taxdump" names.dmp nodes.dmp

"$RSCRIPT_BIN" "$MAKE_DB_R" \
  --assemblies "$SELECTED_ASSEMBLIES" \
  --seq-map "$SEQ_MAP" \
  --seq-lengths "$SEQLENGTHS" \
  --taxdump-dir "${OUTDIR}/taxdump" \
  --outdir "$OUTDIR"

cat > "${OUTDIR}/build_manifest.tsv" <<EOF
key	value
sources	$SOURCES
tax_groups	$TAX_GROUPS
assembly_levels	$ASSEMBLY_LEVELS
refseq_categories	$REFSEQ_CATEGORIES
selected_assemblies	$selected_count
selected_genome_size_bp	$selected_bases
fasta	$FASTA_GZ
index	${OUTDIR}/all_NCBI_genomes.idx
index_batch_size	$INDEX_BATCH_SIZE
EOF

if [[ "$SKIP_INDEX" -eq 0 ]]; then
  echo "Building minimap2 index with -I ${INDEX_BATCH_SIZE}..."
  "${SCRIPT_DIR}/build_minimap2_index.sh" \
    --db-dir "$OUTDIR" \
    --fasta "$FASTA_GZ" \
    --index all_NCBI_genomes.idx \
    --threads "$THREADS" \
    --index-batch-size "$INDEX_BATCH_SIZE" \
    --minimap2 "$MINIMAP2_BIN" \
    --force
else
  echo "Skipping minimap2 index."
fi

echo "Done. Use this in minitax_config.txt:"
echo "  db.dir	$OUTDIR"
echo "  mm2_index	all_NCBI_genomes.idx"
