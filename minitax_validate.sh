#!/usr/bin/env bash
set -euo pipefail

resolve_script_path() {
  local source="${BASH_SOURCE[0]}"
  while [[ -L "$source" ]]; do
    local dir
    dir="$(cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd)"
    source="$(readlink "$source")"
    [[ "$source" != /* ]] && source="$dir/$source"
  done
  cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd
}

script_dir="$(resolve_script_path)"
rscript_bin="${RSCRIPT:-}"

if [[ -z "$rscript_bin" ]]; then
  rscript_bin="$(command -v Rscript || true)"
fi

if [[ -z "$rscript_bin" ]]; then
  cat >&2 <<'EOF'
ERROR: Rscript was not found in PATH.

Activate your R environment first, for example:
  mamba activate minitax

Or pass an explicit Rscript through RSCRIPT:
  RSCRIPT=/path/to/Rscript ./minitax_validate.sh minitax_config.txt
EOF
  exit 1
fi

exec "$rscript_bin" "${script_dir}/scripts/validate_minitax_config.R" "$@"
