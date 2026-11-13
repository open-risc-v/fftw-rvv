#!/usr/bin/env bash
# 构建上游 fftw（autotools，C）→ <out>/prefix。
# 基线用官方发布包 fftw-3.3.10.tar.gz（sha256 56c93254…，含预生成 codelet）；
# git 树的 codelet 要 OCaml genfft + --enable-maintainer-mode，成本高于收益，故不走 git。
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
./configure --prefix="${OUT}/prefix" ${HOST_ARG} \
  --disable-shared --enable-static --disable-fortran --disable-threads \
  --disable-openmp >/dev/null
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" >/dev/null
make install >/dev/null
echo "== 上游 fftw 构建完成：${OUT}/prefix"
