#!/usr/bin/env python
"""
End-to-end multi-colored clock-noise generation from an ADEV specification
==========================================================================

Worked example for three representative clocks:

* a bench OCXO (drift-dominated, no white floor anywhere in the span),
* the AccuBeat USO (Fig. 1 dataset of De Marchi et al. 2024),
* the Orolia RAFS (Fig. 2 dataset of De Marchi et al. 2024).

For each clock the chain is:

1. ADEV specification -> piecewise power-law PSD Sy(f) = h_i f^alpha_i
   (allantools.adev2psd_piecewise_approx)
2. exact round trip PSD -> ADEV (NIST SP 1065 Eq. 65,
   allantools.psd_piecewise_to_adev) -- must reproduce the input
3. time-domain noise realisation (Timmer-Koenig synthesis,
   allantools.noise.timmer_koenig_from_psd)
4. stochastic validation: overlapping ADEV of the generated noise
   (allantools.oadev)

Practical parameter choices encoded below:

* timestep = tau_min / 2, so the synthesis Nyquist frequency
  f_H = 1/(2 dt) covers the highest reconstructed PSD node -- power above
  f_H is truncated (not aliased), so a coarser dt would simply remove the
  short-tau stability from the generated noise;
* duration = 100 * tau_max, so f_1 = 1/T lies well below the lowest PSD
  node -- the synthesis contains no power below 1/T, and the ADEV estimate
  is biased low for tau > T/10, so a shorter record makes the long-tau end
  of the noise ADEV droop below the specification (the redder the
  low-frequency PSD, the longer the record needs to be).

Reference: F. De Marchi, M. K. Plumaris, E. A. Burt and L. Iess,
"An Algorithm to Estimate the Power Spectral Density From Allan Deviation",
IEEE Trans. UFFC 71(4):506-515, 2024.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import allantools as at

CLOCKS = {
    "OCXO": {
        "taus": np.array([1, 10, 1e2, 1e3, 5e3, 1e4]),
        "adevs": np.array([1e-11, 2e-11, 5.5e-11, 1.3e-10, 2.85e-10, 4e-10]),
        "color": "tab:green",
    },
    "AccuBeat USO": {
        "taus": np.array([4, 8, 16, 32, 64, 128, 256, 512, 1.02e3, 2.05e3,
                          4.1e3, 8.19e3, 1.64e4, 3.28e4, 6.55e4]),
        "adevs": np.array([1.07e-13, 1.04e-13, 9.84e-14, 1.04e-13, 1.11e-13,
                           1.21e-13, 1.33e-13, 1.49e-13, 1.60e-13, 2.17e-13,
                           3.40e-13, 6.60e-13, 1.12e-12, 1.97e-12, 3.27e-12]),
        "color": "tab:blue",
    },
    "Orolia RAFS": {
        "taus": np.array([10, 30, 120, 480, 1.92e3, 1.54e4, 6.14e4,
                          1.23e5, 2.46e5, 4.92e5, 9.83e5]),
        "adevs": np.array([1e-12, 3.99e-13, 1.66e-13, 8.85e-14, 4.15e-14,
                           1.62e-14, 8.99e-15, 9.73e-15, 1.17e-14, 1.26e-14,
                           1.36e-14]),
        "color": "tab:pink",
        # tau_max is ~1e6 s: cap the record at 20 tau_max to keep the demo
        # quick; expect a slight droop of the noise ADEV at the last taus.
        "duration_factor": 20,
    },
}

N_SEEDS = 2


def run_clock(name, spec):
    taus, adevs = spec["taus"], spec["adevs"]

    # 1) ADEV -> piecewise power-law PSD
    f_nodes, Sy_nodes, h, alpha = at.adev2psd_piecewise_approx(
        adevs, taus, vartype="adev")

    # 2) exact round trip on a dense tau grid
    tau_dense = np.geomspace(taus[0], taus[-1], 25)
    adev_rt = at.psd_piecewise_to_adev(h, alpha, f_nodes, tau_dense)
    adev_rt_in = at.psd_piecewise_to_adev(h, alpha, f_nodes, taus)
    err = np.abs(adev_rt_in / adevs - 1.0)
    print("%-14s %2d PSD segments | round-trip error at input nodes: "
          "max %.1f%%, median %.1f%%"
          % (name, h.size, 100 * err.max(), 100 * np.median(err)))

    # 3) + 4) noise synthesis and OADEV of the realisations
    timestep = taus[0] / 2.0                       # Nyquist covers highest node
    duration = spec.get("duration_factor", 100) * taus[-1]
    acc = []
    for seed in range(1, N_SEEDS + 1):
        x = at.noise.timmer_koenig_from_psd(
            f_nodes, h, alpha, duration, timestep, output="phase", seed=seed)
        t_noise, adev_noise, _, _ = at.oadev(
            x, rate=1.0 / timestep, data_type="phase", taus=taus)
        acc.append(adev_noise)
    return tau_dense, adev_rt, t_noise, np.mean(acc, axis=0)


def main():
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    for name, spec in CLOCKS.items():
        tau_dense, adev_rt, t_noise, adev_noise = run_clock(name, spec)
        color = spec["color"]
        ax.loglog(spec["taus"], spec["adevs"], "o", ms=6, color=color, zorder=4)
        ax.loglog(tau_dense, adev_rt, "-", lw=1.8, color=color, zorder=3)
        ax.loglog(t_noise, adev_noise, "--", lw=1.2, color=color, alpha=0.85,
                  zorder=2)
        idx = int(0.5 * (len(spec["taus"]) - 1))
        ax.annotate(name, (spec["taus"][idx], spec["adevs"][idx] * 2.4),
                    color=color, fontsize=10, fontweight="bold", ha="center")

    style_handles = [
        Line2D([], [], marker="o", ls="none", color="0.2",
               label="input ADEV specification"),
        Line2D([], [], ls="-", lw=1.8, color="0.2",
               label=r"PSD $\rightarrow$ ADEV round trip (exact integral)"),
        Line2D([], [], ls="--", lw=1.2, color="0.2",
               label="OADEV of synthesised noise (mean of %d seeds)" % N_SEEDS),
    ]
    ax.legend(handles=style_handles, loc="lower left", fontsize=9)
    ax.set_xlabel(r"integration time $\tau$ [s]")
    ax.set_ylabel(r"$\sigma_y(\tau)$ [-]")
    ax.set_title(r"ADEV $\rightarrow$ PSD $\rightarrow$ multi-colored noise:"
                 " loop closure across clock types")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig("multi-color-noise-demo.png", dpi=150)
    print("figure saved to multi-color-noise-demo.png")
    plt.show()


if __name__ == "__main__":
    main()
