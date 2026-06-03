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
MINIMAP2_BIN="${MINIMAP2:-}"
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-dir) DB_DIR="${2:-}"; shift 2 ;;
    --fasta) FASTA="${2:-}"; shift 2 ;;
    --index) INDEX="${2:-}"; shift 2 ;;
    --threads) THREADS="${2:-}"; shift 2 ;;
    --minimap2) MINIMAP2_BIN="${2:-}"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$DB_DIR" ]] || { usage; fail "--db-dir is required."; }
[[ -d "$DB_DIR" ]] || fail "Database directory not found: $DB_DIR"
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || fail "--threads must be a positive integer."

if [[ -z "$MINIMAP2_BIN" ]]; then
  MINIMAP2_BIN="$(command -v minimap2 || true)"
fi
[[ -n "$MINIMAP2_BIN" && -x "$MINIMAP2_BIN" ]] || fail "minimap2 not found. Activate an environment or pass --minimap2 PATH."

FASTA_PATH="$(inside_or_absolute "$DB_DIR" "$FASTA")"
INDEX_PATH="$(inside_or_absolute "$DB_DIR" "$INDEX")"
TMP_INDEX="${INDEX_PATH}.tmp"
LOG_DIR="${DB_DIR}/logs"
LOG_FILE="${LOG_DIR}/build_minimap2_index.log"
PID_FILE="${LOG_DIR}/build_minimap2_index.pid"
STATUS_FILE="${LOG_DIR}/build_minimap2_index.status"

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

{
  echo "started: $(date -Is)"
  echo "minimap2: $MINIMAP2_BIN"
  echo "threads: $THREADS"
  echo "fasta: $FASTA_PATH"
  echo "index: $INDEX_PATH"
  echo "tmp_index: $TMP_INDEX"
  echo
} > "$LOG_FILE"

set +e
"$MINIMAP2_BIN" -t "$THREADS" -d "$TMP_INDEX" "$FASTA_PATH" >> "$LOG_FILE" 2>&1
status=$?
set -e

{
  echo
  echo "ended: $(date -Is)"
  echo "status: $status"
} >> "$LOG_FILE"

echo "$status" > "$STATUS_FILE"

if [[ "$status" -eq 0 ]]; then
  mv -f "$TMP_INDEX" "$INDEX_PATH"
  echo "Index built: $INDEX_PATH" >> "$LOG_FILE"
else
  echo "Index build failed; temporary index kept for inspection: $TMP_INDEX" >> "$LOG_FILE"
fi

exit "$status"
