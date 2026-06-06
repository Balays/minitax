#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="/mnt/c/GitHub/minitax"
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

export PATH="/home/kakuk/anaconda3/envs/minitax/bin:${PATH}"
export RSCRIPT="/home/kakuk/mambaforge/envs/R_4.4/bin/Rscript"

cd "$project_dir"

"${project_dir}/prepare_minitax_5_inputs.sh"
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
