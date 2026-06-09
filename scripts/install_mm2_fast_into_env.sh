#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Install mm2-fast into a conda/mamba environment.

Usage:
  scripts/install_mm2_fast_into_env.sh --env-prefix ENV_PREFIX [options]

Options:
  --env-prefix DIR  Target conda/mamba environment prefix. Required unless
                    CONDA_PREFIX is set.
  --source PATH     Existing mm2-fast binary or source directory containing
                    minimap2. If omitted, the script clones/builds mm2-fast.
  --build-dir DIR   Source/build directory when cloning mm2-fast.
                    Default: $HOME/src/mm2-fast
  --repo-url URL    mm2-fast git URL.
                    Default: https://github.com/bwa-mem2/mm2-fast.git
  --copy            Copy the binary instead of symlinking it.
  -h, --help        Show this help.

The script creates these environment-local commands:
  ENV_PREFIX/bin/mm2-fast
  ENV_PREFIX/bin/minimap2-mm2fast
  ENV_PREFIX/bin/mm2-fast.{avx2,avx,sse42,sse41}
  ENV_PREFIX/bin/minimap2-mm2fast.{avx2,avx,sse42,sse41}

mm2-fast is not distributed as a conda package in this repo. It is built from
source and then linked/copied into the active environment.
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

ENV_PREFIX="${CONDA_PREFIX:-}"
SOURCE="${MM2_FAST_SOURCE:-}"
BUILD_DIR="${MM2_FAST_BUILD_DIR:-${HOME}/src/mm2-fast}"
REPO_URL="${MM2_FAST_REPO_URL:-https://github.com/bwa-mem2/mm2-fast.git}"
COPY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-prefix) ENV_PREFIX="${2:-}"; shift 2 ;;
    --source) SOURCE="${2:-}"; shift 2 ;;
    --build-dir) BUILD_DIR="${2:-}"; shift 2 ;;
    --repo-url) REPO_URL="${2:-}"; shift 2 ;;
    --copy) COPY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown argument: $1" ;;
  esac
done

[[ -n "$ENV_PREFIX" ]] || { usage; fail "--env-prefix is required when CONDA_PREFIX is unset."; }
[[ -d "$ENV_PREFIX/bin" ]] || fail "Environment bin directory not found: $ENV_PREFIX/bin"

# Make conda-provided headers/libraries visible to compilers that do not
# automatically search the active environment. This is required for zlib.h when
# building mm2-fast with conda compilers on some Linux/WSL systems.
[[ -f "$ENV_PREFIX/include/zlib.h" ]] || fail "zlib.h not found in $ENV_PREFIX/include. Install zlib/libzlib in this environment."
export CPATH="$ENV_PREFIX/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="$ENV_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$ENV_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$ENV_PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"

if [[ -z "$SOURCE" ]]; then
  command -v git >/dev/null 2>&1 || fail "git not found. Install git in the active environment."
  command -v make >/dev/null 2>&1 || fail "make not found. Install make in the active environment."
  command -v g++ >/dev/null 2>&1 || fail "g++ not found. Install gxx_linux-64 or a system compiler."

  if [[ ! -d "$BUILD_DIR/.git" ]]; then
    mkdir -p "$(dirname "$BUILD_DIR")"
    git clone --recursive "$REPO_URL" "$BUILD_DIR"
  else
    git -C "$BUILD_DIR" submodule update --init --recursive
  fi

  make -C "$BUILD_DIR"
  SOURCE="${BUILD_DIR%/}/minimap2"
elif [[ -d "$SOURCE" ]]; then
  if [[ ! -x "${SOURCE%/}/minimap2" ]]; then
    command -v make >/dev/null 2>&1 || fail "make not found. Install make in the active environment."
    make -C "$SOURCE"
  fi
  SOURCE="${SOURCE%/}/minimap2"
fi

[[ -x "$SOURCE" ]] || fail "mm2-fast minimap2 binary is not executable: $SOURCE"

SOURCE_DIR="$(cd "$(dirname "$SOURCE")" && pwd)"
SOURCE_BASE="$(basename "$SOURCE")"
SOURCE="${SOURCE_DIR}/${SOURCE_BASE}"

shopt -s nullglob
SOURCE_VARIANTS=("${SOURCE}".*)
shopt -u nullglob
[[ "${#SOURCE_VARIANTS[@]}" -gt 0 ]] || fail "No mm2-fast SIMD executables found next to ${SOURCE}; expected files like ${SOURCE}.avx2"

install_one() {
  local src="$1"
  local dest="$2"
  if [[ "$COPY" -eq 1 ]]; then
    cp "$src" "$dest"
    chmod +x "$dest"
  else
    ln -sfn "$src" "$dest"
  fi
}

for name in mm2-fast minimap2-mm2fast; do
  target="${ENV_PREFIX}/bin/${name}"
  install_one "$SOURCE" "$target"

  for variant in "${SOURCE_VARIANTS[@]}"; do
    suffix="${variant#"$SOURCE"}"
    [[ -x "$variant" ]] || continue
    install_one "$variant" "${target}${suffix}"
  done
done

echo "Installed mm2-fast command links:"
echo "  ${ENV_PREFIX}/bin/mm2-fast -> ${SOURCE}"
echo "  ${ENV_PREFIX}/bin/minimap2-mm2fast -> ${SOURCE}"
echo "  SIMD variants installed with matching command prefixes"
echo
"${ENV_PREFIX}/bin/mm2-fast" --version
