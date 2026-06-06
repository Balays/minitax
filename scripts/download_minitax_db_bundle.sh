#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Download, verify, and extract a published minitax database bundle.

Usage:
  scripts/download_minitax_db_bundle.sh --url-base URL --bundle NAME --outdir DIR

Options:
  --url-base URL   Directory URL containing NAME.SHA256SUMS and archive parts.
                  Works with Zenodo, Hugging Face dataset raw file URLs, or
                  any stable HTTPS file directory.
  --bundle NAME    Bundle basename, for example:
                  minitax_bacteria_refseq_20260605.minimap2-I128G
  --outdir DIR     Directory where the database should be extracted. Required.
  --keep-parts     Keep downloaded archive parts after extraction.
  -h, --help       Show this help.

The script downloads NAME.SHA256SUMS, all listed NAME.tar.gz.part-* files,
verifies SHA256 checksums, concatenates the parts, and extracts the database.
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
  if have curl; then
    curl -L --fail --retry 5 --retry-delay 5 -o "$dest" "$url"
  elif have wget; then
    wget -O "$dest" "$url"
  else
    fail "Need curl or wget."
  fi
}

sha256_check() {
  local sums="$1"
  if have sha256sum; then
    sha256sum -c "$sums"
  elif have shasum; then
    shasum -a 256 -c "$sums"
  else
    fail "Need sha256sum or shasum."
  fi
}

URL_BASE=""
BUNDLE=""
OUTDIR=""
KEEP_PARTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --url-base) URL_BASE="${2:-}"; shift 2 ;;
    --bundle) BUNDLE="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --keep-parts) KEEP_PARTS=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$URL_BASE" ]] || { usage; fail "--url-base is required."; }
[[ -n "$BUNDLE" ]] || { usage; fail "--bundle is required."; }
[[ -n "$OUTDIR" ]] || { usage; fail "--outdir is required."; }

URL_BASE="${URL_BASE%/}"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd -P)"

workdir="${OUTDIR}/.${BUNDLE}.download"
mkdir -p "$workdir"

sums="${workdir}/${BUNDLE}.SHA256SUMS"
download_file "${URL_BASE}/${BUNDLE}.SHA256SUMS" "$sums"

mapfile -t part_files < <(awk -v prefix="${BUNDLE}.tar.gz.part-" '{ name=$2; sub(/^.*\//, "", name); if (index(name, prefix) == 1) print name }' "$sums")
[[ "${#part_files[@]}" -gt 0 ]] || fail "No archive parts found in checksum file: $sums"

for part in "${part_files[@]}"; do
  echo "Downloading: $part" >&2
  download_file "${URL_BASE}/${part}" "${workdir}/${part}"
done

(
  cd "$workdir"
  grep -F "${BUNDLE}.tar.gz.part-" "${BUNDLE}.SHA256SUMS" > "${BUNDLE}.parts.SHA256SUMS"
  sha256_check "${BUNDLE}.parts.SHA256SUMS"
)

echo "Extracting bundle into: $OUTDIR" >&2
part_paths=()
for part in "${part_files[@]}"; do
  part_paths+=("${workdir}/${part}")
done
cat "${part_paths[@]}" | tar -xz -C "$OUTDIR"

if [[ "$KEEP_PARTS" -eq 0 ]]; then
  rm -rf "$workdir"
else
  echo "Kept downloaded parts in: $workdir"
fi

echo "Database bundle extracted: $OUTDIR"
