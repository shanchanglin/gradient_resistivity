# gradient_resistivity

Athena++ setup to reproduce the 2D resistive-MHD reconnection runs T1-T6 of

> S.-C. Lin, Y.-H. Liu, X. Li, *Fast Magnetic Reconnection induced by Resistivity Gradients
> in 2D Magnetohydrodynamics*, arXiv:2109.07526

(resistivity profile `eta(y) = 0.5*eta1*(1 - tanh((y-1)/l_eta)) + eta2`, Eq. 6; Table I).
Boundaries are zero-gradient (`outflow`) on all four sides. Written for NERSC Perlmutter
(CPU nodes), but the serial/local instructions work anywhere with a C++11 compiler.

## Contents

| Path | Purpose |
|---|---|
| `pgen/reconnection_gradient.cpp` | Athena++ problem generator (force-free sheet, uniform rho and p, user resistivity profile, user outputs `Jz`, `eta`, `Ez`) |
| `scripts/gen_inputs.py` | Writes `inputs/athinput.<run>` for T1-T6 (and U1, E1-E5) |
| `inputs/` | Generated Athena++ input files (512^2, HDF5 output) |
| `perlmutter/modules.sh` | Module list shared by the build and run scripts |
| `perlmutter/build.sh` | Clones and builds Athena++ with the problem generator |
| `perlmutter/run.sbatch` | Slurm job script (one node, 64 MPI ranks x 2 OpenMP threads) |
| `analysis/reconnection_rate.py` | Reconnection rate R(t) at the x-line from the HDF5 output |
| `analysis/plot_rates.py` | R(t) for T1-T6 and R vs l_eta (cf. Figs. 3 and 4a) |
| `run_ana.sbatch` | Example job that runs the analysis for T1-T6 |
| `rate/`, `figs/` | Results from the 512^2 Perlmutter runs |

## Physical setup

- Domain `[0,2] x [0,2]`, 512 x 512 cells, 2D, `gamma = 5/3`.
- Force-free sheet centred at `x = 1`:
  `B = B0 tanh((x-1)/lambda) yhat + B0 sech((x-1)/lambda) zhat`, `B0 = 1`, `lambda = 0.04`.
  Because `|B| = B0` everywhere, density and pressure are uniform: `rho = 1`,
  `p = beta*B0^2/2` with `beta = 0.1`, so `V_A = 1`.
- No initial perturbation (`pert_amp = 0`); the resistivity gradient itself breaks the symmetry.
- Resistivity is Ohmic and set in the problem generator. `eta1 = eta_bottom - eta_top`,
  `eta2 = eta_top`.
- Reconnection rate: `R = Ez(x-line) / (B0 V_A)`, with the x-line taken on `x = 1` at the
  zero crossing of `Bx` nearest the peak `|Jz|`.

| Run | l_eta | eta_bottom = eta1+eta2 | eta_top = eta2 | end time |
|---|---|---|---|---|
| T1 | 0.05 | 1e-3 | 2e-4 | 10 |
| T2 | 0.1 | 1e-3 | 2e-4 | 10 |
| T3 | 0.2 | 1e-3 | 2e-4 | 10 |
| T4 | 0.4 | 1e-3 | 2e-4 | 10 |
| T5 | 0.05 | 2e-3 | 1.2e-3 | 5 |
| T6 | 0.05 | 5e-3 | 4.2e-3 | 5 |

## 1. Install Athena++ and build (Perlmutter)

Run on a login node.

```bash
git clone https://github.com/shanchanglin/gradient_resistivity.git
cd gradient_resistivity
git checkout claude/brave-thompson-hqw82e      # or the branch you merged this into
module purge                                     # start from a clean module environment
./perlmutter/build.sh $SCRATCH/athena
```

`build.sh`:
1. loads the modules in `perlmutter/modules.sh` (`PrgEnv-gnu`, `craype-x86-milan`,
   `craype-network-ofi`, `cray-mpich`, `cray-pmi`, `libfabric`, `cray-hdf5-parallel`);
2. clones Athena++ (https://github.com/PrincetonUniversity/athena) into the directory you
   give, which must be new/empty or already an Athena++ checkout;
3. copies `pgen/reconnection_gradient.cpp` into `src/pgen/`;
4. configures and builds with MPI, OpenMP and parallel HDF5, using the Cray `CC` wrapper;
5. sets `CRAY_ADD_RPATH=yes` so the HDF5 path is stored in the binary, and prints the linked
   HDF5 library at the end.

The equivalent manual steps:

```bash
source perlmutter/modules.sh
export CRAY_ADD_RPATH=yes
git clone https://github.com/PrincetonUniversity/athena.git $SCRATCH/athena
cp pgen/reconnection_gradient.cpp $SCRATCH/athena/src/pgen/
cd $SCRATCH/athena
python3 configure.py --prob=reconnection_gradient --coord=cartesian --eos=adiabatic \
    --flux=hlld -b -mpi -omp -hdf5 -h5double --cxx=g++ --mpiccmd=CC
make -j 16
```

Check the result: `source perlmutter/modules.sh && ldd $SCRATCH/athena/bin/athena | grep -E "hdf5|not found"`
should show an HDF5 path and no `not found`.

### Build elsewhere (laptop / workstation, serial)

```bash
git clone https://github.com/PrincetonUniversity/athena.git && cd athena
cp /path/to/gradient_resistivity/pgen/reconnection_gradient.cpp src/pgen/
python3 configure.py --prob=reconnection_gradient --coord=cartesian --eos=adiabatic \
    --flux=hlld -b            # add -mpi -omp for parallel, -hdf5 -h5double for HDF5 output
make -j 4
```

Without HDF5 generate VTK inputs (`gen_inputs.py --format vtk`, see below). VTK output is
one file per MeshBlock, so use a single MeshBlock (`--mb` equal to `--nx`) for the analysis
script. This is practical only at reduced resolution (a 256^2 T1 run to t = 10 takes
about 10 minutes on one core).

## 2. Generate the input files

```bash
python3 scripts/gen_inputs.py T1 T2 T3 T4 T5 T6      # 512^2, HDF5 output, into inputs/
```

Useful options: `--nx 256` (resolution), `--mb 64` (MeshBlock size; 512/64 = 64 blocks, one per
rank), `--format vtk`, `--dtout 0.25` (output interval), `--tlim` (override end time),
`--pert A` (seed perturbation amplitude, default 0). Each input already contains the
history, `prim`, `bcc` and `uov` (Jz, eta, Ez) outputs and a restart block (`rst`, every 1.0).

## 3. Run

Edit `perlmutter/run.sbatch` once and set your NERSC project in `#SBATCH -A` (the file ships
with a placeholder). Then, from the repository root:

```bash
for r in T1 T2 T3 T4 T5 T6; do
  sbatch --export=RUN=$r,ATHENA=$SCRATCH/athena/bin/athena perlmutter/run.sbatch
done
```

- Output goes to `$SCRATCH/recon/<run>/` (`<run>.hst`, `<run>.out2.*.athdf` = prim,
  `<run>.out3.*.athdf` = B, `<run>.out4.*.athdf` = Jz/eta/Ez, `*.rst` restarts).
- The script checks the binary with `ldd` before `srun` and stops with a message if a shared
  library is missing.
- If a job hits its time limit, resubmit the same command: it resumes from the newest `.rst`.
- On the 512^2 runs, each of T1-T4 used about one minute of CPU time on 64 ranks
  (about 17,000 steps to t = 10). A healthy log shows `cycle=... time=... dt=...` lines
  advancing (dt about 1e-3 at the start, about 5e-4 later) and ends with `Terminating on time limit`.
  The `MPICH: freeing inactive persistent request` lines printed at the end are harmless.

To run one case locally instead (serial build, small test):

```bash
mkdir -p run_T1 && cd run_T1
/path/to/athena/bin/athena -i /path/to/gradient_resistivity/inputs/athinput.T1 \
    time/tlim=1.0 mesh/nx1=128 mesh/nx2=128 meshblock/nx1=64 meshblock/nx2=64
```

(Athena++ command-line overrides use `block/name=value` and only work for names that already
exist in the input file.)

## 4. Analyze

Needs `numpy`, `h5py` and `matplotlib` (on Perlmutter `module load python` provides them).

```bash
# reconnection rate for each run; writes <run>_rate.txt and <run>_summary.png
for r in T1 T2 T3 T4 T5 T6; do
  python3 analysis/reconnection_rate.py $SCRATCH/recon/$r $r --athena $SCRATCH/athena
done
# or as a job:  sbatch run_ana.sbatch     (set your project in #SBATCH -A first)

# R(t) for all six runs and R vs l_eta; point it at a directory holding T*_rate.txt
mkdir -p rate figs && cp $SCRATCH/recon/T?/T?_rate.txt rate/
python3 analysis/plot_rates.py rate figs/R_summary.png
```

`--athena` is the Athena++ checkout (used only for `vis/python/athena_read.py`).
Use `--format vtk` for VTK output.

## Results (512^2, Perlmutter)

R is the mean over the last 1.5 time units; plots are in `figs/`.

| Run | R | x-line y |
|---|---|---|
| T1 | 0.043 | 0.958 |
| T2 | 0.036 | 0.918 |
| T3 | 0.029 | 0.832 |
| T4 | 0.023 | 0.610 |
| T5 | 0.022 | 0.956 |
| T6 | 0.023 (still decaying at t = 5) | 0.604 |

- T1 gives R close to the roughly 0.04 reported in the paper, and the x-line sits just on the
  high-eta side of y = 1.
- R falls with l_eta as R ~ l_eta^-0.30 (fit to T1-T4), as in the paper, and not the l_eta^-0.5
  of Eq. 10, which also lies about 2.3x above the data.
- A higher background eta (T1, T5, T6) lowers R, and T6 develops a thick current sheet.
- A 256^2 run gave values within about 5-15% of these.
- The early part of R(t) (t < 2) is not meaningful: the x-line finder can switch between
  zero crossings of Bx while the sheet forms.

## Differences from the paper and caveats

- The paper does not state the time integrator, reconstruction, Riemann solver, CFL number or
  any seed perturbation. These runs use `vl2`, PLM (`xorder = 2`), HLLD, CFL 0.3 and no
  perturbation. Athena++ here has 2 ghost cells (the paper used 4). I expect the trends to be
  robust, but the absolute values are not guaranteed to match.
- The uniform-resistivity comparison run U1 needs a seed perturbation (without one there is no
  x-line). It is generated with `gen_inputs.py U1 --pert 1`, but the perturbation shape and
  amplitude are my choice and were not run on Perlmutter.
- `gen_inputs.py` also writes the E1-E5 inputs (exponential, localized resistivity, Table II);
  they were only partly tested (E3-E5 use super-time-stepping, `sts_integrator = rkl2`).
- `pgen/reconnection_gradient.cpp` requires `<problem>/eta_ohm` to be non-zero in the input
  (Athena++ only enables Ohmic diffusion then); its value is not used.

## Troubleshooting (problems seen on Perlmutter)

| Symptom | Cause / fix |
|---|---|
| `cannot open shared object file: libhdf5_parallel_gnu...so` | Binary built with different modules than the job. Rebuild with `build.sh` (sets rpath) and use the same `modules.sh` |
| `libfabric.so.1 => not found` in `ldd` | Run `source perlmutter/modules.sh` first; `modules.sh` also exports `CRAY_LD_LIBRARY_PATH` |
| `cannot find -lpmi` / `No supported cpu target is set` at link | Load the full module set in `modules.sh` (`craype-x86-milan`, `craype-network-ofi`, `cray-pmi`, ...), ideally after `module purge` |
| `cp: cannot create regular file .../src/pgen/` | The target directory is not an Athena++ checkout; give `build.sh` a new directory |
| `Block name 'outputN' on command line not found` | Command-line overrides need the block to exist in the input file |
| `Unrecognized file format = 'athdf'` | The Athena++ keyword is `file_type = hdf5` (it writes `.athdf` files) |
| `Output variable 'B' not implemented` | Use `variable = bcc` (stored in an HDF5 dataset named `B`) |
