#!/usr/bin/env python3
"""Generate Athena++ input files for the runs in Lin, Liu & Li (arXiv:2109.07526).

Table I  (hyperbolic-tangent eta(y), Eq. 6):   T1-T6, U1 (uniform)
Table II (exponential eta(r) at the x-line, Eq. 12): E1-E5

eta1 and eta2 are the *profile parameters* of Eqs. 6/12, i.e.
    tanh :  eta_bottom = eta1 + eta2 ,  eta_top = eta2
    exp  :  eta_xline  = eta1 + eta2 ,  eta_background = eta2
Table I/II list eta_bottom (or eta_xline) and eta_top (or eta_background), so
eta1 = eta_bottom - eta_top.

usage: gen_inputs.py [--nx 512] [--format athdf|vtk] [--outdir inputs] [runs ...]
"""
import argparse
import os

# name: (profile, eta_at_peak (eta1+eta2), eta2, l_eta, tlim)
RUNS = {
    "T1": (0, 1e-3, 2e-4, 0.05, 10.0),
    "T2": (0, 1e-3, 2e-4, 0.1, 10.0),
    "T3": (0, 1e-3, 2e-4, 0.2, 10.0),
    "T4": (0, 1e-3, 2e-4, 0.4, 10.0),
    "T5": (0, 2e-3, 1.2e-3, 0.05, 5.0),
    "T6": (0, 5e-3, 4.2e-3, 0.05, 5.0),
    "U1": (2, 1e-3, 1e-3, 0.0, 5.0),      # uniform eta = 1e-3
    "E1": (1, 5e-4, 1e-4, 0.05, 5.0),
    "E2": (1, 1e-3, 1e-4, 0.05, 5.0),
    "E3": (1, 1e-2, 1e-4, 0.05, 5.0),
    "E4": (1, 5e-2, 1e-4, 0.05, 5.0),
    "E5": (1, 1.0, 1e-4, 0.05, 5.0),
}

TEMPLATE = """<comment>
problem   = 2D reconnection with a resistivity gradient (run {name})
reference = Lin, Liu & Li, arXiv:2109.07526
configure = --prob=reconnection_gradient --flux=hlld --eos=adiabatic -b [-mpi -omp -hdf5]

<job>
problem_id = {name}

<output1>
file_type = hst
dt        = 0.01

<output2>
file_type   = {fmt}
variable    = prim
dt          = {dtout}

<output3>
file_type   = {fmt}
variable    = {bvar}
dt          = {dtout}

<output4>
file_type   = {fmt}
variable    = uov
dt          = {dtout}

<output5>
file_type = rst
dt        = 1.0

<time>
cfl_number     = 0.3
nlim           = -1
tlim           = {tlim}
integrator     = vl2
xorder         = 2
ncycle_out     = 100
sts_integrator = {sts}

<mesh>
nx1    = {nx}
x1min  = 0.0
x1max  = 2.0
ix1_bc = outflow
ox1_bc = outflow

nx2    = {nx}
x2min  = 0.0
x2max  = 2.0
ix2_bc = outflow
ox2_bc = outflow

nx3    = 1
x3min  = -0.5
x3max  = 0.5
ix3_bc = periodic
ox3_bc = periodic

<meshblock>
nx1 = {mb}
nx2 = {mb}
nx3 = 1

<hydro>
gamma = 1.6666666666667

<problem>
eta_ohm    = {eta_peak:e}   # must be non-zero to enable Ohmic diffusion; value unused
iprofile   = {iprof}        # 0: tanh(y), 1: exp(-r), 2: uniform
eta1       = {eta1:e}
eta2       = {eta2:e}
l_eta      = {leta}
B0         = 1.0
lambda     = 0.04
beta       = 0.1
d0         = 1.0
xc         = 1.0
yc         = 1.0
pert_amp   = {pert}         # flux-perturbation amplitude (units of B0*lambda); 0 = none
pert_width = 0.2
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*", default=list(RUNS))
    ap.add_argument("--nx", type=int, default=512)
    ap.add_argument("--mb", type=int, default=64, help="meshblock size (cells per side)")
    ap.add_argument("--format", default="athdf", choices=["athdf", "vtk"])
    ap.add_argument("--dtout", type=float, default=0.25)
    ap.add_argument("--pert", type=float, default=0.0)
    ap.add_argument("--tlim", type=float, default=None, help="override end time")
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "inputs"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    for name in a.runs:
        iprof, eta_peak, eta2, leta, tlim = RUNS[name]
        eta1 = eta_peak - eta2
        # Explicit Ohmic dt ~ dx^2/(4 eta) becomes prohibitive for large eta: use RKL2 STS
        sts = "rkl2" if eta_peak >= 1e-2 else "none"
        txt = TEMPLATE.format(
            name=name, fmt="hdf5" if a.format == "athdf" else "vtk", dtout=a.dtout, tlim=a.tlim or tlim, sts=sts,
            bvar="bcc", nx=a.nx, mb=min(a.mb, a.nx),
            eta_peak=eta_peak, iprof=iprof, eta1=eta1, eta2=eta2, leta=leta, pert=a.pert)
        path = os.path.join(a.outdir, f"athinput.{name}")
        with open(path, "w") as f:
            f.write(txt)
        print("wrote", os.path.normpath(path))


if __name__ == "__main__":
    main()
