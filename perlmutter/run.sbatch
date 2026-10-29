#!/bin/bash
#SBATCH -A m0000            # <-- your NERSC project
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 04:00:00
#SBATCH -J recon
# usage: sbatch --export=RUN=T1 run.sbatch      (inputs from ../inputs/athinput.$RUN)
# 512^2 with 64^2 MeshBlocks = 64 blocks -> 64 MPI ranks x 2 OpenMP threads on one node.
HERE=${SLURM_SUBMIT_DIR:-$PWD}; [ -f "$HERE/perlmutter/modules.sh" ] || HERE=$(dirname "$HERE")
source "$HERE/perlmutter/modules.sh"
RUN=${RUN:?set RUN=T1|T2|...|E5}
ATHENA=${ATHENA:-$SCRATCH/athena/bin/athena}
OUT=$SCRATCH/recon/$RUN; mkdir -p $OUT; cd $OUT
export OMP_NUM_THREADS=2 OMP_PLACES=threads OMP_PROC_BIND=spread
# Fail once (not on all 64 ranks) if a shared library is missing.
if ldd $ATHENA | grep -q "not found"; then
  echo "ERROR: missing libraries for $ATHENA:"; ldd $ATHENA | grep "not found"
  echo "Rebuild with perlmutter/build.sh (same modules as modules.sh)."; exit 1
fi
RESTART=$(ls -t $OUT/*.rst 2>/dev/null | head -1)
if [ -n "$RESTART" ]; then ARGS="-r $RESTART"; else ARGS="-i $HERE/inputs/athinput.$RUN"; fi
srun -n 64 -c 4 --cpu-bind=cores $ATHENA $ARGS output5/file_type=rst output5/dt=1.0
