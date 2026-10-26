#!/usr/bin/env python3
"""Reconnection rate R = Ez(x-line)/(B0*VA) and snapshots for runs of reconnection_gradient.

The x-line is taken on the symmetry column x = xc as the zero crossing of Bx(y) closest
to the position of peak |Jz|. There v ~ 0 and By = 0, so Ez = eta*Jz.

usage: reconnection_rate.py RUNDIR PROBLEM_ID [--athena PATH_TO_athena] [--format athdf|vtk]
  Outputs expected (from gen_inputs.py): out2 = prim, out3 = B/bcc, out4 = uov.
  (vtk output is one file per MeshBlock; use it only for single-block serial test runs.)
Writes RUNDIR/PROBLEM_ID_rate.txt and RUNDIR/PROBLEM_ID_summary.png
"""
import argparse
import glob
import os
import re
import sys

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("pid")
    ap.add_argument("--athena", default=os.environ.get("ATHENA_DIR", "athena"),
                    help="Athena++ checkout (for vis/python/athena_read.py)")
    ap.add_argument("--format", default="athdf", choices=["athdf", "vtk"])
    ap.add_argument("--xc", type=float, default=1.0)
    ap.add_argument("--yc", type=float, default=1.0)
    a = ap.parse_args()
    sys.path.insert(0, os.path.join(a.athena, "vis", "python"))
    import athena_read as ar

    def files(n):
        pat = (f"{a.pid}.block0.out{n}.*.vtk" if a.format == "vtk"
               else f"{a.pid}.out{n}.*.athdf")
        return sorted(glob.glob(os.path.join(a.rundir, pat)))

    def read(fn):
        if a.format == "vtk":
            xf, yf, zf, d = ar.vtk(fn)
            x, y = 0.5 * (xf[1:] + xf[:-1]), 0.5 * (yf[1:] + yf[:-1])
            out = {}
            for k, v in d.items():
                v = np.asarray(v)
                if v.ndim == 4:                      # vector field
                    for c in range(v.shape[-1]):
                        out[f"{k}{c+1}"] = v[0, :, :, c]
                else:
                    out[k] = v[0]
            with open(fn, "rb") as f:
                t = float(re.search(rb"time=\s*([-+0-9.eE]+)", f.readline() + f.readline()).group(1))
            return t, x, y, out
        d = ar.athdf(fn)
        out = {k: d[k][0] for k in d if isinstance(d[k], np.ndarray) and d[k].ndim == 3}
        return float(d["Time"]), d["x1v"], d["x2v"], out

    fb = files(3)
    fu = files(4)
    if not fb or len(fb) != len(fu):
        sys.exit(f"no/mismatched output files in {a.rundir} for {a.pid}")
    rows = []
    for fB, fU in zip(fb, fu):
        t, x, y, B = read(fB)
        t2, _, _, U = read(fU)
        bx = next(B[k] for k in B if k.lower() in ("bcc1", "bcc"))
        i = np.argmin(abs(x - a.xc))
        cols = [i - 1, i] if x[i] > a.xc else [i, i + 1]
        bxc = bx[:, cols].mean(1)
        jz = U["Jz"][:, cols].mean(1)
        ez = U["Ez"][:, cols].mean(1)
        yj = y[np.argmax(abs(jz))]
        s = np.where(np.sign(bxc[:-1]) != np.sign(bxc[1:]))[0]
        if len(s) == 0:
            rows.append((t, np.nan, np.nan, np.nan)); continue
        s = s[np.argmin(abs(y[s] - yj))]
        w = bxc[s] / (bxc[s] - bxc[s + 1])
        yx = y[s] + w * (y[s + 1] - y[s])
        rows.append((t, yx, ez[s] + w * (ez[s + 1] - ez[s]), jz[s] + w * (jz[s + 1] - jz[s])))
    rows = np.array(rows)
    np.savetxt(os.path.join(a.rundir, f"{a.pid}_rate.txt"), rows,
               header="time y_xline R=Ez/(B0 VA) Jz_xline")
    print("time  y_xline  R  Jz_xline")
    for r in rows[:: max(1, len(rows) // 20)]:
        print("%7.3f %8.4f %9.5f %10.3f" % tuple(r))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t, x, y, B = read(fb[-1])
    _, _, _, U = read(fu[-1])
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].plot(rows[:, 0], rows[:, 2]); ax[0].set(xlabel="t", ylabel="R = Ez/(B0 VA)")
    im = ax[1].pcolormesh(x, y, U["Jz"], cmap="RdBu_r", vmin=-100, vmax=100, shading="auto")
    ax[1].set(title="Jz (final)", xlabel="x", ylabel="y", xlim=(0.5, 1.5)); fig.colorbar(im, ax=ax[1])
    ax[2].plot(U["eta"][:, len(x) // 2], y); ax[2].set(xlabel="eta(x=xc,y)", ylabel="y")
    fig.tight_layout(); fig.savefig(os.path.join(a.rundir, f"{a.pid}_summary.png"), dpi=130)


if __name__ == "__main__":
    main()
