#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${ROOT_DIR}/metis_mpi_minimal"
CASE_DIR="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8"
OWNER="${CASE_DIR}/constant/polyMesh/owner"
NEIGHBOUR="${CASE_DIR}/constant/polyMesh/neighbour"
LOG="${ROOT_DIR}/log.metis_mpi_minimal_zb64_$(date +%Y%m%d_%H%M%S)"

: "${MPIEXEC:=mpirun}"
: "${NP:=64}"
: "${PARMETIS_DIR:=/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install}"

export LD_LIBRARY_PATH="${PARMETIS_DIR}/lib:${LD_LIBRARY_PATH:-}"

{
  echo "date=$(date --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "bin=${BIN}"
  echo "owner=${OWNER}"
  echo "neighbour=${NEIGHBOUR}"
  echo "np=${NP}"
  echo "mpiexec=${MPIEXEC}"
  echo "parmetis_dir=${PARMETIS_DIR}"
  echo "ldd:"
  ldd "${BIN}" | grep -E 'metis|mpi|GKlib' || true
  echo
  /usr/bin/time -p "${MPIEXEC}" -np "${NP}" "${BIN}" \
    --owner "${OWNER}" \
    --neighbour "${NEIGHBOUR}" \
    --parts "${NP}" \
    --fpe \
    --verbose
} 2>&1 | tee "${LOG}"

echo "log saved to ${LOG}"
