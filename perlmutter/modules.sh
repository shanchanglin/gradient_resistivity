# Sourced by build.sh and run.sbatch so the build and run environments match.
# Explicit CPU toolchain (Perlmutter CPU nodes are AMD Milan). Loading an already-loaded
# module is harmless.
module load PrgEnv-gnu
module load craype-x86-milan        # CPU target for the CC wrapper
module load craype-network-ofi      # network target (brings in the PMI/MPI link paths)
module load cray-mpich
module load cray-pmi
module load libfabric 2>/dev/null || true   # may already come with craype-network-ofi
module load cray-hdf5-parallel
module list
# Runtime library paths for the Cray MPI stack (libfabric etc.) as set by the modules above.
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH:+$CRAY_LD_LIBRARY_PATH:}$LD_LIBRARY_PATH
