#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${project_dir}/../.." && pwd)"
config="${CONFIG_FILE:-${project_dir}/minitax_bacteria_refseq_5.config.tsv}"
mode="${1:-all}"

case "$mode" in
  all|--all) mode="all" ;;
  validate|--validate-only) mode="validate" ;;
  dry-run|--dry-run) mode="dry-run" ;;
  map|--map-only) mode="map" ;;
  classify|--classify-only) mode="classify" ;;
  *)
    echo "Usage: $0 [--all|--validate-only|--dry-run|--map-only|--classify-only]" >&2
    exit 2
    ;;
esac

RSCRIPT="${RSCRIPT:-$(command -v Rscript || true)}"
[[ -n "$RSCRIPT" ]] || { echo "Rscript not found in PATH; activate the minitax environment or set RSCRIPT." >&2; exit 1; }

cd "$project_dir"

config_indir="$(awk -F '\t' '$1 == "indir" { print $2; exit }' "$config")"
if [[ "$config_indir" == "fastq_5" || "$config_indir" == "./fastq_5" || "$config_indir" == "${project_dir}/fastq_5" ]]; then
  "${project_dir}/prepare_minitax_5_inputs.sh"
fi
"${project_dir}/place_minitax_5_output_on_d.sh" "$config"

case "$mode" in
  validate)
    "$repo_dir/minitax_validate.sh" "$config" --step all
    ;;
  dry-run)
    "$repo_dir/minitax_validate.sh" "$config" --step all
    "$repo_dir/minitax.sh" "$config" --dry-run
    ;;
  map)
    "$repo_dir/minitax_validate.sh" "$config" --step map
    "$repo_dir/minitax.sh" "$config"
    ;;
  classify)
    "$repo_dir/minitax_validate.sh" "$config" --step classify
    "$RSCRIPT" "$repo_dir/minitax.complete.R" "$config"
    ;;
  all)
    "$repo_dir/minitax_validate.sh" "$config" --step all
    "$repo_dir/minitax.sh" "$config"
    "$repo_dir/minitax_validate.sh" "$config" --step classify
    "$RSCRIPT" "$repo_dir/minitax.complete.R" "$config"
    ;;
esac
