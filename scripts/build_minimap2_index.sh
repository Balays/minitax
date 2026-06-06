#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Build a minimap2 index safely for an existing minitax database.

Usage:
  scripts/build_minimap2_index.sh --db-dir DIR [options]

Options:
  --db-dir DIR     Database directory. Required.
  --fasta FILE     FASTA/FASTA.gz path or filename inside db-dir.
                   Default: all_NCBI_genomes.fna.gz
  --index FILE     Index path or filename inside db-dir.
                   Default: all_NCBI_genomes.idx
  --threads N      minimap2 threads. Default: 8
  --index-batch-size SIZE
                   minimap2 -I batch size. Default: 128G.
                   Use a value larger than the reference size to avoid
                   multi-part indexes with unreliable MAPQ.
  --minimap2 PATH  minimap2 executable. Default: command -v minimap2
  --force          Rebuild even if the final index exists.
  -h, --help       Show this help.

The script writes:
  logs/build_minimap2_index.log
  logs/build_minimap2_index.pid
  logs/build_minimap2_index.status

The index is first written to INDEX.tmp and moved into place only on success.
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

inside_or_absolute() {
  local db_dir="$1"
  local path="$2"
  if [[ "$path" == /* ]]; then
    printf '%s' "$path"
  else
    printf '%s/%s' "$db_dir" "$path"
  fi
}

DB_DIR=""
FASTA="all_NCBI_genomes.fna.gz"
INDEX="all_NCBI_genomes.idx"
THREADS=8
INDEX_BATCH_SIZE="128G"
MINIMAP2_BIN="${MINIMAP2:-}"
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-dir) DB_DIR="${2:-}"; shift 2 ;;
    --fasta) FASTA="${2:-}"; shift 2 ;;
    --index) INDEX="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --index-batch-size) INDEX_BATCH_SIZE="${2:-}"; shift 2 ;;
    --minimap2) MINIMAP2_BIN="${2:-}"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$DB_DIR" ]] || { usage; fail "--db-dir is required."; }
[[ -d "$DB_DIR" ]] || fail "Database directory not found: $DB_DIR"
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || fail "--threads must be a positive integer."
[[ "$INDEX_BATCH_SIZE" =~ ^[0-9]+[KkMmGgTt]?$ ]] || fail "--index-batch-size must look like 128G, 64000M, or 500000000."

if [[ -z "$MINIMAP2_BIN" ]]; then
  MINIMAP2_BIN="$(command -v minimap2 || true)"
elif [[ "$MINIMAP2_BIN" != */* ]]; then
  MINIMAP2_BIN="$(command -v "$MINIMAP2_BIN" || true)"
fi
[[ -n "$MINIMAP2_BIN" && -x "$MINIMAP2_BIN" ]] || fail "minimap2 not found. Activate an environment or pass --minimap2 PATH."

FASTA_PATH="$(inside_or_absolute "$DB_DIR" "$FASTA")"
INDEX_PATH="$(inside_or_absolute "$DB_DIR" "$INDEX")"
TMP_INDEX="${INDEX_PATH}.tmp"
LOG_DIR="${DB_DIR}/logs"
LOG_FILE="${LOG_DIR}/build_minimap2_index.log"
PID_FILE="${LOG_DIR}/build_minimap2_index.pid"
STATUS_FILE="${LOG_DIR}/build_minimap2_index.status"
MANIFEST="${DB_DIR}/build_manifest.tsv"

[[ -s "$FASTA_PATH" ]] || fail "FASTA not found or empty: $FASTA_PATH"
mkdir -p "$LOG_DIR"

if [[ -s "$INDEX_PATH" && "$FORCE" -eq 0 ]]; then
  echo "Index already exists: $INDEX_PATH"
  echo "0" > "$STATUS_FILE"
  exit 0
fi

echo "$$" > "$PID_FILE"
echo "running" > "$STATUS_FILE"
rm -f "$TMP_INDEX"

sha256_file() {
  local path="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$path" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$path" | awk '{print $1}'
  else
    printf 'NA'
  fi
}

manifest_set() {
  local key="$1"
  local value="$2"
  local tmp="${MANIFEST}.tmp.$$"
  if [[ ! -f "$MANIFEST" ]]; then
    printf 'key\tvalue\n' > "$MANIFEST"
  fi
  awk -F '\t' -v OFS='\t' -v key="$key" -v value="$value" '
    NR == 1 && $1 == "key" { print; next }
    $1 == key { print key, value; seen = 1; next }
    { print }
    END { if (!seen) print key, value }
  ' "$MANIFEST" > "$tmp"
  mv -f "$tmp" "$MANIFEST"
}

minimap2_version="$("$MINIMAP2_BIN" --version 2>/dev/null || true)"
[[ -n "$minimap2_version" ]] || minimap2_version="unknown"
build_command="$MINIMAP2_BIN -t $THREADS -I $INDEX_BATCH_SIZE -d $TMP_INDEX $FASTA_PATH"

{
  echo "started: $(date -Is)"
  echo "minimap2: $MINIMAP2_BIN"
  echo "minimap2_version: $minimap2_version"
  echo "threads: $THREADS"
  echo "index_batch_size: $INDEX_BATCH_SIZE"
  echo "fasta: $FASTA_PATH"
  echo "index: $INDEX_PATH"
  echo "tmp_index: $TMP_INDEX"
  echo "command: $build_command"
  echo
} > "$LOG_FILE"

set +e
"$MINIMAP2_BIN" -t "$THREADS" -I "$INDEX_BATCH_SIZE" -d "$TMP_INDEX" "$FASTA_PATH" >> "$LOG_FILE" 2>&1
status=$?
set -e

index_parts=0
if [[ "$status" -eq 0 ]]; then
  index_parts="$(grep -c 'loaded/built the index' "$LOG_FILE" || true)"
  if [[ "$index_parts" -gt 1 ]]; then
    status=50
    {
      echo
      echo "ERROR: minimap2 built $index_parts index parts."
      echo "ERROR: Multi-part indexes make MAPQ unreliable for minitax; increase --index-batch-size."
    } >> "$LOG_FILE"
  fi
fi

{
  echo
  echo "ended: $(date -Is)"
  echo "index_parts: $index_parts"
  echo "status: $status"
} >> "$LOG_FILE"

echo "$status" > "$STATUS_FILE"

if [[ "$status" -eq 0 ]]; then
  mv -f "$TMP_INDEX" "$INDEX_PATH"
  echo "Index built: $INDEX_PATH" >> "$LOG_FILE"
  build_date="$(date -Is)"
  fasta_sha256="$(sha256_file "$FASTA_PATH")"
  index_sha256="$(sha256_file "$INDEX_PATH")"
  manifest_set "index_batch_size" "$INDEX_BATCH_SIZE"
  manifest_set "minimap2_version" "$minimap2_version"
  manifest_set "minimap2_index_command" "$build_command"
  manifest_set "fasta_sha256" "$fasta_sha256"
  manifest_set "index_sha256" "$index_sha256"
  manifest_set "index_build_date" "$build_date"
  manifest_set "index_parts" "$index_parts"
else
  echo "Index build failed; temporary index kept for inspection: $TMP_INDEX" >> "$LOG_FILE"
fi

exit "$status"
