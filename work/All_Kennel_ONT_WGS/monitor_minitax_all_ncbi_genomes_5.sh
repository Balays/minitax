#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CONFIG_FILE="${project_dir}/minitax_all_ncbi_genomes_5.config.tsv"

exec "${project_dir}/monitor_minitax_5.sh" "$@"
