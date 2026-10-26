# gradient_resistivity

Athena++ setup reproducing Lin, Liu & Li, *Fast Magnetic Reconnection induced by Resistivity
Gradients in 2D MHD* (arXiv:2109.07526).

- `pgen/reconnection_gradient.cpp` – problem generator (force-free sheet, uniform ρ and p,
  user resistivity profile: tanh(y) Eq. 6, exp(-r) Eq. 12, or uniform). Outputs Jz, eta, Ez.
- `scripts/gen_inputs.py` – writes `inputs/athinput.{T1..T6,U1,E1..E5}` (Tables I, II); all
  boundaries are `outflow` (zero gradient). Runs with eta >= 1e-2 use RKL2 super-time-stepping.
- `perlmutter/build.sh`, `perlmutter/run.sbatch` – build and run on Perlmutter CPU nodes.
- `analysis/reconnection_rate.py` – R = Ez(x-line)/(B0 VA) vs time, plus a summary figure.

## Perlmutter
    ./perlmutter/build.sh                       # clones/builds Athena++ in $SCRATCH/athena
    sbatch --export=RUN=T1 perlmutter/run.sbatch
    python3 analysis/reconnection_rate.py $SCRATCH/recon/T1 T1 --athena $SCRATCH/athena

## Validation (serial, 256^2, half the paper's resolution)
T1 reaches R ≈ 0.05 (paper: ≈0.04 simulated, 0.07 averaged-equations), x-line at y ≈ 0.96,
on the high-eta side. Unperturbed uniform eta (U1) has no x-line; use `--pert 1` for U1.

## Differences from the paper / caveats
- Athena++ here uses 2 ghost cells (paper: 4); integrator vl2 + PLM + HLLD, CFL 0.3.
  The paper's exact integrator/Riemann solver and initial perturbation are not stated.
- Large-eta runs (E3-E5) may need tuning of STS (`sts_max_dt_ratio`) – not tested here.
- The MPI/HDF5 build in `perlmutter/build.sh` was NOT tested (no Perlmutter access here);
  check module names and the configure line.
