#!/usr/bin/env bash

if [ -z "${BASH_VERSION:-}" ]; then
    echo "error: run this script with bash, not sh" >&2
    exit 1
fi

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    echo "error: do not source this script; run it with 'bash makeParMETIS.sh'" >&2
    return 1
fi

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${WORK_DIR:-$PWD}"
BASE_DIR="$(cd "${WORK_DIR}" && pwd)"
cd "${BASE_DIR}"

: "${INSTALL_PREFIX:=${BASE_DIR}/parmetis-install}"
: "${JOBS:=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)}"
: "${CLEAN_SOURCE:=1}"
: "${BACKUP_INSTALL:=1}"

require_env()
{
    local name="$1"
    if [[ -z "${!name:-}" ]]; then
        echo "error: required environment variable '${name}' is not set" >&2
        exit 1
    fi
}


require_cmd()
{
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "error: required command '${cmd}' not found in PATH" >&2
        exit 1
    fi
}


find_tarball()
{
    local name
    for name in "$@"; do
        if [[ -f "${BASE_DIR}/${name}" ]]; then
            printf '%s\n' "${BASE_DIR}/${name}"
            return 0
        fi
    done

    return 1
}


tar_root_dir()
{
    tar -tzf "$1" | head -1 | cut -d/ -f1
}


backup_path()
{
    local path="$1"
    local stamp="$2"
    if [[ -e "${path}" ]]; then
        mv "${path}" "${path}.bak.${stamp}"
    fi
}


run_make_config()
{
    local src_dir="$1"
    shift

    (
        cd "${src_dir}"
        make distclean >/dev/null 2>&1 || true
        make config "$@"
    )
}


run_make_build_install()
{
    local src_dir="$1"

    (
        cd "${src_dir}"
        make -j"${JOBS}"
        make install
    )
}


require_env WM_CC
require_env WM_CXX
require_cmd tar
require_cmd make
require_cmd cmake

CC_BIN="$(command -v "${WM_CC}")"
CXX_BIN="$(command -v "${WM_CXX}")"
require_cmd "${WM_CC}"
require_cmd "${WM_CXX}"

CC_NAME="$(basename "${CC_BIN}")"
MPI_CC_NAME="mpi${CC_NAME}"

if command -v "${MPI_CC_NAME}" >/dev/null 2>&1; then
    MPICC_BIN="$(command -v "${MPI_CC_NAME}")"
elif command -v mpiicc >/dev/null 2>&1 \
    && [[ "$(mpiicc -show 2>/dev/null | awk '{print $1}')" == "${CC_NAME}" ]]; then
    MPICC_BIN="$(command -v mpiicc)"
elif command -v mpicc >/dev/null 2>&1 \
    && [[ "$(mpicc -show 2>/dev/null | awk '{print $1}')" == "${CC_NAME}" ]]; then
    MPICC_BIN="$(command -v mpicc)"
else
    echo "error: cannot find an MPI C wrapper that matches WM_CC=${WM_CC}" >&2
    echo "tried '${MPI_CC_NAME}', then mpiicc/mpicc compiler-front matching" >&2
    exit 1
fi

GKLIB_TAR="$(find_tarball gklib-master.tar.gz GKlib-master.tar.gz)"
METIS_TAR="$(find_tarball metis-master.tar.gz METIS-master.tar.gz)"
PARMETIS_TAR="$(find_tarball parmetis-main.tar.gz ParMETIS-main.tar.gz)"

if [[ -z "${GKLIB_TAR:-}" || -z "${METIS_TAR:-}" || -z "${PARMETIS_TAR:-}" ]]; then
    echo "error: missing one or more required tarballs in ${BASE_DIR}" >&2
    echo "required: gklib-master.tar.gz, metis-master.tar.gz, parmetis-main.tar.gz" >&2
    exit 1
fi

GKLIB_SRC="${BASE_DIR}/$(tar_root_dir "${GKLIB_TAR}")"
METIS_SRC="${BASE_DIR}/$(tar_root_dir "${METIS_TAR}")"
PARMETIS_SRC="${BASE_DIR}/$(tar_root_dir "${PARMETIS_TAR}")"

if [[ "${CLEAN_SOURCE}" == "1" ]]; then
    rm -rf "${GKLIB_SRC}" "${METIS_SRC}" "${PARMETIS_SRC}"
fi

tar -xzf "${GKLIB_TAR}"
tar -xzf "${METIS_TAR}"
tar -xzf "${PARMETIS_TAR}"

STAMP="$(date +%Y%m%d_%H%M%S)"
if [[ "${BACKUP_INSTALL}" == "1" ]]; then
    backup_path "${INSTALL_PREFIX}" "${STAMP}"
fi
mkdir -p "${INSTALL_PREFIX}"

LABEL_SIZE="${WM_LABEL_SIZE:-32}"
if [[ "${WM_OPTIONS:-}" == *Int64* ]]; then
    LABEL_SIZE=64
fi

REAL_SIZE=32
if [[ -n "${WM_DP:-}" || "${WM_OPTIONS:-}" == *DP* || "${WM_PRECISION_OPTION:-}" == "DP" ]]; then
    REAL_SIZE=64
fi

METIS_WIDTH_FLAGS=()
if [[ "${LABEL_SIZE}" == "64" ]]; then
    METIS_WIDTH_FLAGS+=(i64=1)
fi
if [[ "${REAL_SIZE}" == "64" ]]; then
    METIS_WIDTH_FLAGS+=(r64=1)
fi

echo "=== OpenFOAM toolchain summary ==="
echo "SCRIPT_DIR=${SCRIPT_DIR}"
echo "WORK_DIR=${WORK_DIR}"
echo "WM_CC=${WM_CC}"
echo "WM_CXX=${WM_CXX}"
echo "WM_OPTIONS=${WM_OPTIONS:-unknown}"
echo "CC_BIN=${CC_BIN}"
echo "CXX_BIN=${CXX_BIN}"
echo "MPICC_BIN=${MPICC_BIN}"
echo "LABEL_SIZE=${LABEL_SIZE}"
echo "REAL_SIZE=${REAL_SIZE}"
echo "INSTALL_PREFIX=${INSTALL_PREFIX}"
echo

echo "=== Extracted sources ==="
echo "GKLIB_SRC=${GKLIB_SRC}"
echo "METIS_SRC=${METIS_SRC}"
echo "PARMETIS_SRC=${PARMETIS_SRC}"
echo

echo "=== Building GKlib ==="
run_make_config \
    "${GKLIB_SRC}" \
    cc="${CC_BIN}" \
    shared=1 \
    prefix="${INSTALL_PREFIX}"
run_make_build_install "${GKLIB_SRC}"

if [[ -d "${INSTALL_PREFIX}/lib64" && ! -e "${INSTALL_PREFIX}/lib" ]]; then
    ln -s lib64 "${INSTALL_PREFIX}/lib"
fi

echo "=== Building METIS ==="
run_make_config \
    "${METIS_SRC}" \
    cc="${CC_BIN}" \
    shared=1 \
    prefix="${INSTALL_PREFIX}" \
    gklib_path="${INSTALL_PREFIX}" \
    "${METIS_WIDTH_FLAGS[@]}"
run_make_build_install "${METIS_SRC}"

if [[ -d "${INSTALL_PREFIX}/lib64" && ! -e "${INSTALL_PREFIX}/lib" ]]; then
    ln -s lib64 "${INSTALL_PREFIX}/lib"
fi

echo "=== Building ParMETIS ==="
run_make_config \
    "${PARMETIS_SRC}" \
    cc="${MPICC_BIN}" \
    shared=1 \
    prefix="${INSTALL_PREFIX}" \
    gklib_path="${INSTALL_PREFIX}" \
    metis_path="${INSTALL_PREFIX}"
run_make_build_install "${PARMETIS_SRC}"

if [[ -d "${INSTALL_PREFIX}/lib64" && ! -e "${INSTALL_PREFIX}/lib" ]]; then
    ln -s lib64 "${INSTALL_PREFIX}/lib"
fi

LIB_DIR="${INSTALL_PREFIX}/lib"
if [[ ! -d "${LIB_DIR}" && -d "${INSTALL_PREFIX}/lib64" ]]; then
    LIB_DIR="${INSTALL_PREFIX}/lib64"
fi

echo
echo "=== Installation complete ==="
echo "Headers:"
ls "${INSTALL_PREFIX}/include" | grep -E "metis|parmetis"
echo
echo "Libraries (${LIB_DIR}):"
ls "${LIB_DIR}" | grep -E "metis|parmetis|GKlib"
echo
echo "metis.h widths:"
grep -E "^#define (IDXTYPEWIDTH|REALTYPEWIDTH)" "${INSTALL_PREFIX}/include/metis.h"
echo
echo "ldd libmetis.so:"
ldd "${LIB_DIR}/libmetis.so" | grep -E "GKlib|imf|svml|irng|intlc" || true
echo
echo "Done."
