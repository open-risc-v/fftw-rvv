#!/usr/bin/env bash
# 真实包用例：fftw 1024 点复 FFT，输出量化结果 + 误差指标。
set -euo pipefail
: "${STAGE:=run}" "${SIDE:=A}" "${BIN:=build/native}" "${DATA:=build/ab-data}" "${OUT:=build/ab-data/out}"
mkdir -p "$DATA" "$OUT"
case "$STAGE" in
  gen) : > "$DATA/fftw-input.txt" ;;
  run) "$BIN/fftw_harness" "$OUT/fftw.bin" > "$OUT/hash.txt" ;;
  *) echo "unknown STAGE=$STAGE" >&2; exit 2 ;;
esac
