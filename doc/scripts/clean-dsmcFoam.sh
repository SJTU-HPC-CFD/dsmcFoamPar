#!/bin/bash
# Clean dsmcFoam+ build artifacts from hyStrath_dlb local directory.
# Usage:
#   source doc/scripts/clean-dsmcFoam.sh
#   source doc/scripts/clean-dsmcFoam.sh --force
#   bash doc/scripts/clean-dsmcFoam.sh --dry-run

set -e

HYSTRATH_DLB="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FORCE=0
DRY_RUN=0

for arg in "$@"
do
    case "$arg" in
        -f|--force)
            FORCE=1
            ;;
        -n|--dry-run)
            DRY_RUN=1
            ;;
        -h|--help)
            echo "Usage: source ${BASH_SOURCE[0]} [--force|--dry-run]"
            return 0 2>/dev/null || exit 0
            ;;
        *)
            echo "Unknown argument: $arg" 1>&2
            return 2 2>/dev/null || exit 2
            ;;
    esac
done

# Source the same base environment used by build-dsmcFoam.sh.
if [ -f "$HOME/intel/oneapi/setvars.sh" ]
then
    source "$HOME/intel/oneapi/setvars.sh" --force 2>/dev/null || true
fi

source "$HOME/code/OpenFoam/OF-1706/OpenFOAM-v1706/etc/bashrc"

# Override paths to keep cleanup scoped to this checkout.
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/bin
export PARMETIS_DIR=${PARMETIS_DIR:-/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install}
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$PARMETIS_DIR/lib:$LD_LIBRARY_PATH

run_cmd()
{
    echo "+ $*"
    if [ "$DRY_RUN" -eq 0 ]
    then
        "$@"
    fi
}

remove_glob()
{
    local pattern="$1"
    local matches

    shopt -s nullglob
    matches=( $pattern )
    shopt -u nullglob

    if [ "${#matches[@]}" -eq 0 ]
    then
        echo "  not found: $pattern"
        return 0
    fi

    run_cmd rm -f "${matches[@]}"
}

echo "=== Clean targets ==="
echo "  LIBBIN: $FOAM_USER_LIBBIN"
echo "  APPBIN: $FOAM_USER_APPBIN"
echo "  USER_DIR: $WM_PROJECT_USER_DIR"
echo ""

if [ "$DRY_RUN" -eq 0 ] && [ "$FORCE" -eq 0 ]
then
    echo "This removes local hyStrath_dlb wmake objects, lnInclude dirs, libraries, and binaries."
    echo "Type CLEAN to continue:"
    read confirm
    if [ "$confirm" != "CLEAN" ]
    then
        echo "Aborted."
        return 1 2>/dev/null || exit 1
    fi
fi

# Clean in reverse build order.
echo "=== Cleaning dsmcInitialise+ utility ==="
run_cmd wclean "$HYSTRATH_DLB/applications/utilities/preProcessing/dsmc/dsmcInitialise+"

echo "=== Cleaning dsmcFoam+ solver ==="
run_cmd wclean "$HYSTRATH_DLB/applications/solvers/discreteMethods/dsmc/dsmcFoam+"

echo "=== Cleaning libdsmcFoam+ ==="
run_cmd wclean libso "$HYSTRATH_DLB/src/lagrangian/dsmc"

echo "=== Cleaning libgeneralMolecule ==="
run_cmd wclean libso "$HYSTRATH_DLB/src/lagrangian/molecularDynamics/general"

echo "=== Cleaning liblagrangian+ ==="
run_cmd wclean libso "$HYSTRATH_DLB/src/lagrangian/basic"

echo "=== Removing local binaries ==="
remove_glob "$FOAM_USER_APPBIN/dsmcFoam+"
remove_glob "$FOAM_USER_APPBIN/dsmcInitialise+"

echo "=== Removing local libraries ==="
remove_glob "$FOAM_USER_LIBBIN/libdsmcFoam+.*"
remove_glob "$FOAM_USER_LIBBIN/libdsmcFoam+.so"
remove_glob "$FOAM_USER_LIBBIN/libgeneralMolecule.*"
remove_glob "$FOAM_USER_LIBBIN/libgeneralMolecule.so"
remove_glob "$FOAM_USER_LIBBIN/liblagrangian+.*"
remove_glob "$FOAM_USER_LIBBIN/liblagrangian+.so"

echo ""
echo "=== Done ==="
echo "Cleaned dsmcFoam+ hyStrath_dlb build artifacts."

cd "$HYSTRATH_DLB"
