#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Package a runnable minitax database bundle for publication.

Usage:
  scripts/package_minitax_db.sh --db-dir DIR [options]

Options:
  --db-dir DIR       Database directory to package. Required.
  --name NAME        Bundle name. Default: basename(db-dir)_YYYYMMDD.minimap2-I128G
  --outdir DIR       Output directory. Default: db-dir/bundles
  --part-size SIZE   Split archive part size. Default: 45G.
                    Use a value below the target host's per-file limit.
  -h, --help         Show this help.

The bundle contains the FASTA, minimap2 index, taxonomy tables, sequence lengths,
build manifest, and index build log. Downloaded raw genome files are not included.
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

have() {
  command -v "$1" >/dev/null 2>&1
}

sha256_file() {
  local path="$1"
  if have sha256sum; then
    sha256sum "$path"
  elif have shasum; then
    shasum -a 256 "$path"
  else
    fail "Need sha256sum or shasum."
  fi
}

DB_DIR=""
NAME=""
OUTDIR=""
PART_SIZE="45G"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-dir) DB_DIR="${2:-}"; shift 2 ;;
    --name) NAME="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --part-size) PART_SIZE="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$DB_DIR" ]] || { usage; fail "--db-dir is required."; }
[[ -d "$DB_DIR" ]] || fail "Database directory not found: $DB_DIR"
[[ "$PART_SIZE" =~ ^[0-9]+[KkMmGgTt]?$ ]] || fail "--part-size must look like 45G or 500M."

DB_DIR="$(cd "$DB_DIR" && pwd -P)"
date_stamp="$(date +%Y%m%d)"
if [[ -z "$NAME" ]]; then
  NAME="$(basename "$DB_DIR")_${date_stamp}.minimap2-I128G"
fi
if [[ -z "$OUTDIR" ]]; then
  OUTDIR="${DB_DIR}/bundles"
fi
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd -P)"

required_files=(
  all_NCBI_genomes.fna.gz
  all_NCBI_genomes.idx
  NCBI.db.tsv
  NCBI.db.uni.tsv
  NCBI.db.genomesize.tsv
  NCBI_genome_collection_seqlengths.txt
  build_manifest.tsv
)
optional_files=(
  NCBI.db.uni.spec.tsv
  logs/build_minimap2_index.log
)

for file in "${required_files[@]}"; do
  [[ -s "${DB_DIR}/${file}" ]] || fail "Required database file missing or empty: ${DB_DIR}/${file}"
done

include_files=("${required_files[@]}")
for file in "${optional_files[@]}"; do
  if [[ -s "${DB_DIR}/${file}" ]]; then
    include_files+=("$file")
  fi
done

manifest="${OUTDIR}/${NAME}.bundle_manifest.tsv"
component_sha="${OUTDIR}/${NAME}.component_sha256.tsv"
part_sha="${OUTDIR}/${NAME}.SHA256SUMS"
part_prefix="${OUTDIR}/${NAME}.tar.gz.part-"

rm -f "${part_prefix}"* "$manifest" "$component_sha" "$part_sha"

{
  printf 'key\tvalue\n'
  printf 'bundle_name\t%s\n' "$NAME"
  printf 'created_at\t%s\n' "$(date -Is)"
  printf 'db_dir\t%s\n' "$DB_DIR"
  printf 'archive_format\ttar.gz split by %s\n' "$PART_SIZE"
  printf 'included_files\t%s\n' "$(IFS=';'; echo "${include_files[*]}")"
} > "$manifest"

{
  printf 'sha256\tfile\n'
  for file in "${include_files[@]}"; do
    sha256_file "${DB_DIR}/${file}" | awk -v file="$file" '{print $1 "\t" file}'
  done
} > "$component_sha"

echo "Creating split archive parts in: $OUTDIR" >&2
tar -C "$DB_DIR" -czf - "${include_files[@]}" | split -b "$PART_SIZE" - "$part_prefix"

{
  for part in "${part_prefix}"*; do
    [[ -f "$part" ]] || continue
    sha256_file "$part" | awk '{name=$2; sub(/^.*\//, "", name); print $1 "  " name}'
  done
  sha256_file "$manifest" | awk '{name=$2; sub(/^.*\//, "", name); print $1 "  " name}'
  sha256_file "$component_sha" | awk '{name=$2; sub(/^.*\//, "", name); print $1 "  " name}'
} > "$part_sha"

echo "Bundle created:"
echo "  manifest: $manifest"
echo "  component checksums: $component_sha"
echo "  part checksums: $part_sha"
echo "  parts:"
ls -lh "${part_prefix}"*
