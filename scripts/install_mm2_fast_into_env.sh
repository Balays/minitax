#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Install a locally built mm2-fast binary into a conda/mamba environment.

Usage:
  scripts/install_mm2_fast_into_env.sh --env-prefix ENV_PREFIX [options]

Options:
  --env-prefix DIR  Target conda/mamba environment prefix. Required unless
                    CONDA_PREFIX is set.
  --source PATH     mm2-fast binary or source directory containing minimap2.
                    Default: /mnt/c/ubuntu/programs/mm2-fast/minimap2
  --copy            Copy the binary instead of symlinking it.
  -h, --help        Show this help.

The script creates these environment-local commands:
  ENV_PREFIX/bin/mm2-fast
  ENV_PREFIX/bin/minimap2-mm2fast
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

ENV_PREFIX="${CONDA_PREFIX:-}"
SOURCE="${MM2_FAST_SOURCE:-/mnt/c/ubuntu/programs/mm2-fast/minimap2}"
COPY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-prefix) ENV_PREFIX="${2:-}"; shift 2 ;;
    --source) SOURCE="${2:-}"; shift 2 ;;
    --copy) COPY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$ENV_PREFIX" ]] || { usage; fail "--env-prefix is required when CONDA_PREFIX is unset."; }
[[ -d "$ENV_PREFIX/bin" ]] || fail "Environment bin directory not found: $ENV_PREFIX/bin"

if [[ -d "$SOURCE" ]]; then
  SOURCE="${SOURCE%/}/minimap2"
fi

[[ -x "$SOURCE" ]] || fail "mm2-fast minimap2 binary is not executable: $SOURCE"

for name in mm2-fast minimap2-mm2fast; do
  target="${ENV_PREFIX}/bin/${name}"
  if [[ "$COPY" -eq 1 ]]; then
    cp "$SOURCE" "$target"
    chmod +x "$target"
  else
    ln -sfn "$SOURCE" "$target"
  fi
done

echo "Installed mm2-fast command links:"
echo "  ${ENV_PREFIX}/bin/mm2-fast -> ${SOURCE}"
echo "  ${ENV_PREFIX}/bin/minimap2-mm2fast -> ${SOURCE}"
echo
"${ENV_PREFIX}/bin/mm2-fast" --version
