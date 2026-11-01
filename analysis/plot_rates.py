#!/usr/bin/env python3
"""R(t) for T1-T6 and R vs l_eta (cf. Fig. 3 and 4a of arXiv:2109.07526).
usage: plot_rates.py [RATE_DIR] [OUT_PNG]   (RATE_DIR holds T*_rate.txt from reconnection_rate.py)
R vs l_eta uses T1-T4, averaged over the last 1.5 time units; the Eq. 10 prediction
R = sqrt(eta_bottom/(2 l_eta)) (mu = VA = 1) and the fitted power law are overplotted."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rd = sys.argv[1] if len(sys.argv) > 1 else "rate"
out = sys.argv[2] if len(sys.argv) > 2 else "figs/R_summary.png"
# run: (l_eta, eta_bottom, eta_top)
P = {"T1": (0.05, 1e-3, 2e-4), "T2": (0.1, 1e-3, 2e-4), "T3": (0.2, 1e-3, 2e-4),
     "T4": (0.4, 1e-3, 2e-4), "T5": (0.05, 2e-3, 1.2e-3), "T6": (0.05, 5e-3, 4.2e-3)}
col = {"T1": "#1f77b4", "T2": "#ff7f0e", "T3": "#2ca02c", "T4": "#d62728",
       "T5": "#9467bd", "T6": "#8c564b"}
D = {r: np.loadtxt(f"{rd}/{r}_rate.txt") for r in P}

def late(d, w=1.5):
    m = (d[:, 0] >= d[:, 0].max() - w) & np.isfinite(d[:, 2])
    return d[m, 2].mean()

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for r, d in D.items():
    ok = np.isfinite(d[:, 2]) & (d[:, 0] > 0.3)   # drop t=0 and the initial transient
    ax[0].plot(d[ok, 0], d[ok, 2], color=col[r], label=f"{r}  (l$_\\eta$={P[r][0]}, "
               f"$\\eta_{{b}}$={P[r][1]:g}, $\\eta_{{t}}$={P[r][2]:g})")
ax[0].set(xlabel="t [system size / Alfvén speed]", ylabel=r"R = E$_z$(x-line)/(B$_0$V$_A$)",
          title="(a) Reconnection rate", ylim=(0, 0.065))
ax[0].legend(fontsize=7.5, frameon=False)

runs = ["T1", "T2", "T3", "T4"]
l = np.array([P[r][0] for r in runs]); R = np.array([late(D[r]) for r in runs])
s, c = np.polyfit(np.log(l), np.log(R), 1)
lf = np.geomspace(0.04, 0.5, 50)
ax[1].loglog(l, R, "o", color="#1f77b4", ms=7, label="simulation (T1-T4)")
ax[1].loglog(lf, np.sqrt(1e-3/(2*lf)), color="#ff7f0e", label=r"Eq. 10: $\sqrt{\eta_{bottom}/2l_\eta}$")
ax[1].loglog(lf, np.exp(c)*lf**s, color="#2ca02c", label=f"fit  R $\\propto$ l$_\\eta^{{{s:.2f}}}$")
ax[1].set(xlabel=r"l$_\eta$", ylabel="R (late-time mean)", title="(b) Reconnection rate vs gradient scale")
ax[1].legend(frameon=False, fontsize=8)
for a in ax:
    a.grid(alpha=.25, which="both")
fig.tight_layout(); fig.savefig(out, dpi=140)
print("slope %.3f" % s, " R:", dict(zip(runs, np.round(R, 4))))
