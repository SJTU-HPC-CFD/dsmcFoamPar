#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${ROOT_DIR}/metis_mpi_minimal"
SRC="${ROOT_DIR}/metis_mpi_minimal.cpp"
: "${PARMETIS_DIR:=/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install}"

if command -v mpiicpc >/dev/null 2>&1; then
  MPICXX=mpiicpc
elif command -v mpiicpx >/dev/null 2>&1; then
  MPICXX=mpiicpx
else
  echo "error: neither mpiicpc nor mpiicpx found in current environment" >&2
  exit 1
fi

export LD_LIBRARY_PATH="${PARMETIS_DIR}/lib:${LD_LIBRARY_PATH:-}"

"${MPICXX}" \
  -O2 -g -std=c++11 \
  -I"${PARMETIS_DIR}/include" \
  "${SRC}" \
  -L"${PARMETIS_DIR}/lib" \
  -Wl,-rpath,"${PARMETIS_DIR}/lib" \
  -lmetis -lGKlib \
  -o "${OUT}"

echo "built ${OUT} with MPICXX=${MPICXX}"
