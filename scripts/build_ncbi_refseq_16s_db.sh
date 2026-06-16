#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build an NCBI RefSeq TargetedLoci 16S database for minitax.

Outputs:
  ncbi_refseq_16s.fna.gz
  ncbi_refseq_16s.idx
  NCBI.db.tsv
  NCBI.db.uni.tsv
  NCBI.db.uni.spec.tsv
  db_data.tsv
  NCBI_16S_seqlengths.txt
  metadata/ncbi_refseq_16s_accession_taxid.tsv

Usage:
  scripts/build_ncbi_refseq_16s_db.sh --outdir DIR [options]

Options:
  --outdir DIR              Output database directory. Required.
  --kingdoms LIST           TargetedLoci kingdoms. Default: Bacteria,Archaea
  --min-length N            Keep FASTA records >= N bp. Default: 0. Use 1200 for near-full-length V1-V9.
  --threads N               Threads for pigz/minimap2 indexing. Default: 4
  --index-batch-size SIZE   minimap2 -I batch size. Default: 128G
  --minimap2 PATH           minimap2 path. Default: command -v minimap2
  --taxonkit PATH           taxonkit path. Default: command -v taxonkit
  --skip-index              Do not build minimap2 index.
  --force                   Rebuild combined FASTA and tables.
  -h, --help                Show this help.

Example:
  scripts/build_ncbi_refseq_16s_db.sh \
    --outdir /mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s \
    --min-length 1200 \
    --threads 20

Use in minitax_config.txt:
  db         ncbi_refseq_16S
  db.dir     /path/to/ncbi_targeted_loci_16s
  mm2_ref    ncbi_refseq_16s.fna.gz
  mm2_index  ncbi_refseq_16s.idx
EOF
}

fail() { echo "ERROR: $*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

download_file() {
  local url="$1" dest="$2"
  if [[ -s "$dest" ]]; then echo "Using existing: $dest"; return 0; fi
  echo "Downloading: $url"
  if have curl; then curl -L --fail --retry 5 --retry-delay 5 -o "$dest" "$url"
  elif have wget; then wget -O "$dest" "$url"
  else fail "Need curl or wget for downloads."
  fi
}

csv_to_lines() { printf '%s' "$1" | awk -v RS=',' '{gsub(/^[ \t]+|[ \t]+$/, ""); if ($0 != "") print}'; }

OUTDIR=""; KINGDOMS="Bacteria,Archaea"; MIN_LENGTH=0; THREADS=4; INDEX_BATCH_SIZE="128G"
MINIMAP2_BIN="${MINIMAP2:-}"; TAXONKIT_BIN="${TAXONKIT:-}"; SKIP_INDEX=0; FORCE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --kingdoms) KINGDOMS="${2:-}"; shift 2 ;;
    --min-length) MIN_LENGTH="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --index-batch-size) INDEX_BATCH_SIZE="${2:-}"; shift 2 ;;
    --minimap2) MINIMAP2_BIN="${2:-}"; shift 2 ;;
    --taxonkit) TAXONKIT_BIN="${2:-}"; shift 2 ;;
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
if [[ -z "$TAXONKIT_BIN" ]]; then TAXONKIT_BIN="$(command -v taxonkit || true)"; elif [[ "$TAXONKIT_BIN" != */* ]]; then TAXONKIT_BIN="$(command -v "$TAXONKIT_BIN" || true)"; fi
[[ -n "$TAXONKIT_BIN" && -x "$TAXONKIT_BIN" ]] || fail "taxonkit not found. Install with: mamba install -c bioconda -c conda-forge taxonkit"
if [[ "$SKIP_INDEX" -eq 0 && -z "$MINIMAP2_BIN" ]]; then MINIMAP2_BIN="$(command -v minimap2 || true)"; elif [[ "$SKIP_INDEX" -eq 0 && "$MINIMAP2_BIN" != */* ]]; then MINIMAP2_BIN="$(command -v "$MINIMAP2_BIN" || true)"; fi
if [[ "$SKIP_INDEX" -eq 0 ]]; then [[ -n "$MINIMAP2_BIN" && -x "$MINIMAP2_BIN" ]] || fail "minimap2 not found. Activate an environment or pass --minimap2 PATH."; fi

mkdir -p "$OUTDIR"/{downloads,metadata,logs,taxdump}
RAW_FASTA="${OUTDIR}/ncbi_refseq_16s.raw.fna"
FILTERED_FASTA="${OUTDIR}/ncbi_refseq_16s.fna"
FASTA_GZ="${FILTERED_FASTA}.gz"
SEQLENGTHS="${OUTDIR}/NCBI_16S_seqlengths.txt"
TAXDUMP="${OUTDIR}/taxdump/taxdump.tar.gz"
ACCESSION_TAXID="${OUTDIR}/metadata/ncbi_refseq_16s_accession_taxid.tsv"
ACCESSION_TAXID_RAW="${OUTDIR}/metadata/ncbi_refseq_16s_accession_taxid.raw.tsv"
TAXIDS="${OUTDIR}/metadata/ncbi_refseq_16s.taxids.txt"
LINEAGES="${OUTDIR}/metadata/ncbi_refseq_16s.taxid_lineages.tsv"

if [[ "$FORCE" -eq 1 ]]; then
  rm -f "$RAW_FASTA" "$FILTERED_FASTA" "$FASTA_GZ" "$SEQLENGTHS" "$ACCESSION_TAXID" "$ACCESSION_TAXID_RAW" "$TAXIDS" "$LINEAGES" \
        "${OUTDIR}/NCBI.db.tsv" "${OUTDIR}/NCBI.db.uni.tsv" "${OUTDIR}/NCBI.db.uni.spec.tsv" "${OUTDIR}/db_data.tsv" "${OUTDIR}/ncbi_refseq_16s.idx"
fi
: > "$RAW_FASTA"
printf 'seqnames\taccession\ttaxid\tsuperkingdom\torganism_name\n' > "$ACCESSION_TAXID_RAW"

while read -r kingdom; do
  [[ -n "$kingdom" ]] || continue
  case "$kingdom" in
    Bacteria|Archaea) ;;
    *) fail "Unsupported TargetedLoci kingdom for 16S builder: $kingdom. Expected Bacteria or Archaea." ;;
  esac
  kingdom_lc="$(printf '%s' "$kingdom" | tr '[:upper:]' '[:lower:]')"
  fna="${kingdom_lc}.16SrRNA.fna.gz"; gbff="${kingdom_lc}.16SrRNA.gbff.gz"
  base_url="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/${kingdom}"
  download_file "${base_url}/${fna}" "${OUTDIR}/downloads/${fna}"
  download_file "${base_url}/${gbff}" "${OUTDIR}/downloads/${gbff}"
  gzip -cd "${OUTDIR}/downloads/${fna}" >> "$RAW_FASTA"
  gzip -cd "${OUTDIR}/downloads/${gbff}" | awk -v OFS='\t' -v superkingdom="$kingdom" '
    function flush_record(){ if(seqname!="" && taxid!="") print seqname, accession, taxid, superkingdom, organism }
    /^LOCUS/ { flush_record(); seqname=""; accession=""; taxid=""; organism=""; in_source=0; next }
    /^ACCESSION[ ]+/ { accession=$0; sub(/^ACCESSION[ ]+/,"",accession); split(accession,a,/[ \t]/); accession=a[1]; next }
    /^VERSION[ ]+/ { version=$0; sub(/^VERSION[ ]+/,"",version); split(version,v,/[ \t]/); seqname=v[1]; next }
    /^  ORGANISM[ ]+/ { organism=$0; sub(/^  ORGANISM[ ]+/,"",organism); gsub(/^[ \t]+|[ \t]+$/,"",organism); next }
    /^     source[ ]+/ { in_source=1; next }
    /^     [A-Za-z_]+[ ]+/ && $0 !~ /^     source[ ]+/ { in_source=0 }
    in_source && /\/db_xref="taxon:[0-9]+"/ { taxid=$0; sub(/^.*taxon:/,"",taxid); sub(/".*$/,"",taxid); next }
    /^\/\// { flush_record(); seqname=""; accession=""; taxid=""; organism=""; in_source=0 }
    END { flush_record() }
  ' >> "$ACCESSION_TAXID_RAW"
done < <(csv_to_lines "$KINGDOMS")

awk -v min_len="$MIN_LENGTH" -v seq_lengths="$SEQLENGTHS" '
  function flush_record(){ if(name!="" && len>=min_len){ print ">" name; for(i=1;i<=nseq;i++) print seq[i]; print name "\t" len >> seq_lengths } delete seq; nseq=0; len=0 }
  BEGIN { print "seqnames\tseqlengths" > seq_lengths }
  /^>/ { flush_record(); header=$0; sub(/^>/,"",header); split(header,fields,/[ \t]/); name=fields[1]; next }
  { gsub(/\r/,""); if($0!=""){ nseq++; seq[nseq]=$0; len+=length($0) } }
  END { flush_record() }
' "$RAW_FASTA" > "$FILTERED_FASTA"

awk -F '\t' -v OFS='\t' 'NR==FNR{if(NR>1) keep[$1]=1; next} FNR==1{print; next} $1 in keep' "$SEQLENGTHS" "$ACCESSION_TAXID_RAW" > "$ACCESSION_TAXID"
cut -f3 "$ACCESSION_TAXID" | tail -n +2 | sort -u > "$TAXIDS"
[[ -s "$TAXIDS" ]] || fail "No TaxIDs could be parsed from the GBFF files."

download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "$TAXDUMP"
tar -xzf "$TAXDUMP" -C "${OUTDIR}/taxdump" names.dmp nodes.dmp merged.dmp
"$TAXONKIT_BIN" lineage --data-dir "${OUTDIR}/taxdump" "$TAXIDS" \
  | "$TAXONKIT_BIN" reformat --data-dir "${OUTDIR}/taxdump" -F -f '{p}\t{c}\t{o}\t{f}\t{g}\t{s}' \
  | awk -F '\t' -v OFS='\t' 'BEGIN{print "taxid","lineage","phylum","class","order","family","genus","species"} {print $1,$2,$3,$4,$5,$6,$7,$8}' > "$LINEAGES"

awk -F '\t' -v OFS='\t' '
  NR==FNR {
    if (FNR > 1) lin[$1] = $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8
    next
  }
  FNR==1 {print "seqnames","ident","taxid","superkingdom","phylum","class","order","family","genus","species"; next}
  {
    seq=$1; tax=$3; superkingdom=$4; org=$5
    if (tax in lin) print seq,seq,tax,superkingdom,lin[tax]
    else print seq,seq,tax,superkingdom,"","","","","",org
  }
' "$LINEAGES" "$ACCESSION_TAXID" > "${OUTDIR}/NCBI.db.tsv"

awk -F '\t' -v OFS='\t' 'NR==1{next} {print $1,$3,$4,$5,$6,$7,$8,$9,$10}' "${OUTDIR}/NCBI.db.tsv" | sort -t $'\t' -u > "${OUTDIR}/NCBI.db.uni.tsv.tmp"
{ printf 'seqnames\ttaxid\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n'; cat "${OUTDIR}/NCBI.db.uni.tsv.tmp"; } > "${OUTDIR}/NCBI.db.uni.tsv"
rm -f "${OUTDIR}/NCBI.db.uni.tsv.tmp"

awk -F '\t' -v OFS='\t' 'NR==1{next} {print $4,$5,$6,$7,$8,$9,$10}' "${OUTDIR}/NCBI.db.tsv" | sort -t $'\t' -u > "${OUTDIR}/NCBI.db.uni.spec.tsv.tmp"
{ printf 'superkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n'; cat "${OUTDIR}/NCBI.db.uni.spec.tsv.tmp"; } > "${OUTDIR}/NCBI.db.uni.spec.tsv"
rm -f "${OUTDIR}/NCBI.db.uni.spec.tsv.tmp"

cp "${OUTDIR}/NCBI.db.uni.tsv" "${OUTDIR}/db_data.tsv"
if have pigz; then pigz -f -p "$THREADS" "$FILTERED_FASTA"; else gzip -f "$FILTERED_FASTA"; fi

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
  "${SCRIPT_DIR}/build_minimap2_index.sh" --db-dir "$OUTDIR" --fasta "$FASTA_GZ" --index ncbi_refseq_16s.idx --threads "$THREADS" --index-batch-size "$INDEX_BATCH_SIZE" --minimap2 "$MINIMAP2_BIN" --force
else
  echo "Skipping minimap2 index."
fi

echo "Done. Use this in minitax_config.txt:"
echo -e "  db\tncbi_refseq_16S"
echo -e "  db.dir\t$OUTDIR"
echo -e "  mm2_ref\tncbi_refseq_16s.fna.gz"
echo -e "  mm2_index\tncbi_refseq_16s.idx"
