#!/usr/bin/env bash
# 交叉编译基线：容器内把 baseline 编成 riscv64，探针用 qemu-user 跑，比对两侧哈希。
# 上游是 autotools（fftw）：--host=riscv64-linux-gnu。
set -euo pipefail
cd "$(dirname "$0")/.."

IMAGE=${IMAGE:-rvv-cross:2}
MARCH=${MARCH:-rv64gc}
OUT=build/riscv64-b
HARNESS=${HARNESS:-fftw_harness}

docker run --rm -v "$PWD":/work -w /work "${IMAGE}" bash -lc "
set -e
rm -rf ${OUT}
mkdir -p ${OUT}
HOST=riscv64-linux-gnu \
CC=riscv64-linux-gnu-gcc \
CFLAGS='-march=${MARCH}' \
./cmake/pkg-build.sh /work/baseline /work/${OUT}/pkg > /tmp/pkg.log 2>&1 || { tail -30 /tmp/pkg.log; exit 1; }

cmake -S . -B ${OUT}/harness -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=/work/cmake/toolchain-riscv64-linux.cmake \
  -DCMAKE_C_FLAGS='-march=${MARCH}' \
  -DPKG_PREFIX=/work/${OUT}/pkg/prefix > /tmp/h.log 2>&1 || { tail -30 /tmp/h.log; exit 1; }
cmake --build ${OUT}/harness -j > /tmp/hb.log 2>&1 || { tail -20 /tmp/hb.log; exit 1; }

mkdir -p ${OUT}/out
qemu-riscv64 -L /usr/riscv64-linux-gnu ${OUT}/harness/${HARNESS} ${OUT}/out/fftw.bin > ${OUT}/out/hash.txt \
  || { echo 'qemu 运行失败（或 max_rel 超限）'; cat ${OUT}/out/hash.txt; exit 1; }
echo '--- riscv64（${MARCH}）---'; cat ${OUT}/out/hash.txt
echo '--- native 基线 ---'; cat build/ab-data/out/01_fftw/a/hash.txt
"
