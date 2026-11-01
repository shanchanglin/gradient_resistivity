#!/bin/bash
#SBATCH -A m2114            # <-- your NERSC project
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J recon

module load python
for r in T1 T2 T3 T4 T5 T6; do
  python3 analysis/reconnection_rate.py $SCRATCH/recon/$r $r --athena $SCRATCH/athena
done
