#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cat >&2 <<'EOF'
This legacy downloader only fetched FASTA sequences and did not build the
metadata files minitax needs. Use the new database builder instead.

Example:
  scripts/build_ncbi_minitax_db.sh --outdir /mnt/d/data/databases/all_NCBI_genomes --threads 32

Forwarding your arguments to scripts/build_ncbi_minitax_db.sh ...
EOF

exec "${script_dir}/scripts/build_ncbi_minitax_db.sh" "$@"
