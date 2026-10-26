#!/bin/bash
# Build Athena++ with the reconnection_gradient problem generator on Perlmutter (CPU nodes).
# usage: ./build.sh [ATHENA_DIR]   (default: $SCRATCH/athena)
set -e
ATHENA_DIR=${1:-$SCRATCH/athena}
HERE=$(cd "$(dirname "$0")/.." && pwd)

module load PrgEnv-gnu cray-hdf5-parallel
[ -d "$ATHENA_DIR" ] || git clone https://github.com/PrincetonUniversity/athena.git "$ATHENA_DIR"
cp "$HERE/pgen/reconnection_gradient.cpp" "$ATHENA_DIR/src/pgen/"
cd "$ATHENA_DIR"
# 'CC' is the Cray MPI-aware C++ wrapper: --cxx=g++ picks GNU flags, then CC overrides the compiler.
python3 configure.py --prob=reconnection_gradient --coord=cartesian --eos=adiabatic \
    --flux=hlld -b -mpi -omp -hdf5 --cxx=g++ --include=$HDF5_DIR/include --lib_path=$HDF5_DIR/lib \
    --lib=hdf5 --mpiccmd=CC -h5double 2>/dev/null || \
python3 configure.py --prob=reconnection_gradient --coord=cartesian --eos=adiabatic \
    --flux=hlld -b -mpi -omp -hdf5 --cxx=g++ --mpiccmd=CC
make clean && make -j 16
echo "built: $ATHENA_DIR/bin/athena"
