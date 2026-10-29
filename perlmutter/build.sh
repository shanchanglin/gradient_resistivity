#!/bin/bash
# Build Athena++ with the reconnection_gradient problem generator on Perlmutter (CPU nodes).
# usage: ./perlmutter/build.sh [ATHENA_DIR]     (default: $SCRATCH/athena; no angle brackets)
set -ex
ATHENA_DIR=${1:-$SCRATCH/athena}
HERE=$(cd "$(dirname "$0")/.." && pwd)

source "$HERE/perlmutter/modules.sh"
# Embed library paths (incl. HDF5) in the binary so it also runs when the job's modules differ.
export CRAY_ADD_RPATH=yes

if [ ! -f "$ATHENA_DIR/configure.py" ]; then
  if [ -e "$ATHENA_DIR" ] && [ -n "$(ls -A "$ATHENA_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $ATHENA_DIR exists but is not an Athena++ checkout (no configure.py)."
    echo "Use an empty/new directory, e.g. ./perlmutter/build.sh \$SCRATCH/athena"; exit 1
  fi
  git clone https://github.com/PrincetonUniversity/athena.git "$ATHENA_DIR"
fi
cp "$HERE/pgen/reconnection_gradient.cpp" "$ATHENA_DIR/src/pgen/"
cd "$ATHENA_DIR"
# CC is the Cray compiler wrapper (adds MPI and HDF5 from the loaded modules automatically).
python3 configure.py --prob=reconnection_gradient --coord=cartesian --eos=adiabatic \
    --flux=hlld -b -mpi -omp -hdf5 -h5double --cxx=g++ --mpiccmd=CC
make clean
make -j 16
ldd bin/athena | grep -i hdf5
echo "built: $ATHENA_DIR/bin/athena"
