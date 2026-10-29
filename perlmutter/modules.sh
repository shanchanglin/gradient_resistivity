# Sourced by build.sh and run.sbatch so the build and run environments match.
module load cpu              # CPU-node toolchain defaults (Perlmutter logins default to GPU)
module load PrgEnv-gnu
module load cray-hdf5-parallel
