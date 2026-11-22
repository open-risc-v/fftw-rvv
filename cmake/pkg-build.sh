#!/usr/bin/env bash
# 构建上游 fftw（autotools，C）→ <out>/prefix。
# 基线用官方发布包 fftw-3.3.10.tar.gz（sha256 56c93254…，含预生成 codelet）；
# git 树的 codelet 要 OCaml genfft + --enable-maintainer-mode，成本高于收益，故不走 git。
# A/B 靠 CFLAGS（-march 与 -DRVV_KERNELS）驱动：autotools 只认环境 CFLAGS，
# 拼成 configure 参数会被吞掉（speexdsp/mpg123 踩过这个坑）。
# 关闭 threads / fortran / shared：基线只要确定性的 double 变换，多线程会让 plan 不可复现。
set -euo pipefail
SRC=$1
OUT=$2
rm -rf "${OUT}"
mkdir -p "${OUT}/prefix" "${OUT}/src"
cp -R "${SRC}" "${OUT}/src/pkg"
cd "${OUT}/src/pkg"
HOST_ARG=""
[ -n "${HOST:-}" ] && HOST_ARG="--host=${HOST}"
MARCH=${3:-rv64gc}
KERN=${4:-0}
export CFLAGS="-O2 -march=${MARCH}"
if [ "${KERN}" = "1" ]; then
  export CFLAGS="${CFLAGS} -DRVV_KERNELS=1"
fi
./configure --prefix="${OUT}/prefix" ${HOST_ARG} \
  --disable-shared --enable-static --disable-fortran --disable-threads \
  --disable-openmp >/dev/null
# 容器里没有 texinfo，doc 目标会去调 makeinfo 直接失败；把它置空（我们不要 info 手册）
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" MAKEINFO=true >/dev/null
make install MAKEINFO=true >/dev/null
echo "== 上游 fftw 构建完成：${OUT}/prefix"
