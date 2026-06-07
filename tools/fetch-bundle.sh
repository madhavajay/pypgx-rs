#!/usr/bin/env bash
#
# Fetch the pypgx-bundle reference data and assemble a self-contained bundle for
# pypgx-rs: the 1000 Genomes phasing panels (Beagle ref=) plus the CNV models.
#
# PyPGx ships its CNV callers as pickled scikit-learn models, which cannot load
# in Rust. By default this script converts each to the Rust-native weight form
# (`data.json`, see tools/convert_cnv_models_all.py) so pypgx-rs can run CNV with
# no Python at runtime. Use --raw to keep the pickles (no Python needed).
#
# The resulting layout mirrors the upstream bundle so `predict_cnv` /
# `estimate_phase_beagle` resolve it directly:
#   $DEST/1kgp/{GRCh37,GRCh38}/<gene>.vcf.gz(.tbi)
#   $DEST/cnv/{GRCh37,GRCh38}/<gene>.zip
#   $DEST/VERSION
#
# Usage:   tools/fetch-bundle.sh [DEST]            (default DEST=./pypgx-bundle)
# Env:     PYPGX_BUNDLE_VERSION  bundle git branch/tag        (default 0.26.0)
#          BUNDLE_REPO           bundle git URL               (default sbslee/pypgx-bundle)
#          PYTHON                python with pypgx+sklearn    (default python3)
# Flags:   --raw            keep pickled CNV models (skip conversion)
#          --no-cnv         omit CNV models entirely (panels only)
#          --assembly A     only GRCh37 or GRCh38 (default: both)
set -euo pipefail

DEST="${1:-./pypgx-bundle}"
[[ "${DEST:0:1}" == "-" ]] && DEST="./pypgx-bundle"   # first arg was a flag
VERSION="${PYPGX_BUNDLE_VERSION:-0.26.0}"
REPO="${BUNDLE_REPO:-https://github.com/sbslee/pypgx-bundle}"
PYTHON="${PYTHON:-python3}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONVERT=1; WITH_CNV=1; ONLY_ASM=""
for arg in "$@"; do
  case "$arg" in
    --raw) CONVERT=0 ;;
    --no-cnv) WITH_CNV=0 ;;
    --assembly) NEXT_IS_ASM=1 ;;
    GRCh37|GRCh38) [[ "${NEXT_IS_ASM:-0}" == "1" ]] && ONLY_ASM="$arg" && NEXT_IS_ASM=0 ;;
  esac
done
ASSEMBLIES=(GRCh37 GRCh38)
[[ -n "$ONLY_ASM" ]] && ASSEMBLIES=("$ONLY_ASM")

echo ">> pypgx-bundle $VERSION -> $DEST  (cnv: convert=$CONVERT with_cnv=$WITH_CNV; assemblies: ${ASSEMBLIES[*]})"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
echo ">> cloning $REPO @ $VERSION (shallow)"
git clone --depth 1 --branch "$VERSION" "$REPO" "$TMP/bundle"

mkdir -p "$DEST"
for asm in "${ASSEMBLIES[@]}"; do
  echo ">> copying 1kgp/$asm panels"
  mkdir -p "$DEST/1kgp/$asm"
  cp -f "$TMP/bundle/1kgp/$asm"/*.vcf.gz "$DEST/1kgp/$asm/" 2>/dev/null || true
  cp -f "$TMP/bundle/1kgp/$asm"/*.tbi    "$DEST/1kgp/$asm/" 2>/dev/null || true
done

if [[ "$WITH_CNV" == "1" ]]; then
  if [[ "$CONVERT" == "1" ]]; then
    echo ">> converting CNV models to Rust weights (this unpickles sklearn in Python)"
    # Convert per assembly so --assembly limits the work.
    for asm in "${ASSEMBLIES[@]}"; do
      mkdir -p "$TMP/src/$asm"
      cp -f "$TMP/bundle/cnv/$asm"/*.zip "$TMP/src/$asm/" 2>/dev/null || true
    done
    "$PYTHON" "$SCRIPT_DIR/convert_cnv_models_all.py" "$TMP/src" "$DEST/cnv" --verify
  else
    echo ">> copying raw (pickled) CNV models"
    for asm in "${ASSEMBLIES[@]}"; do
      mkdir -p "$DEST/cnv/$asm"
      cp -f "$TMP/bundle/cnv/$asm"/*.zip "$DEST/cnv/$asm/" 2>/dev/null || true
    done
  fi
fi

cp -f "$TMP/bundle/VERSION" "$DEST/VERSION" 2>/dev/null || echo "$VERSION" > "$DEST/VERSION"

echo ">> done. bundle size:"
du -sh "$DEST"
find "$DEST" -maxdepth 2 -type d | sort
