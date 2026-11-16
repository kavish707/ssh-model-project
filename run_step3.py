#!/usr/bin/env python3
"""
Project A -- Step 3 driver: localisation and coherence, measured separately.

    python run_step3.py            # full grid (~8 min)
    python run_step3.py --quick    # full ED only, no L = 16/18 Lanczos

Outputs
-------
data/step3_grid.csv            one row per (L, delta, V, method): z_0, xi_W,
                               n_eff, window status, sum-rule error
data/step3_window_scan.csv     z_0 against the window factor 0.6 … 0.9
data/step3_lanczos_conv.csv    100 / 200 / 400 Lanczos vectors vs full ED
data/step3_lines_*.csv         lines and weights at the edge and centre sites
data/step3_envelope_*.csv      W_j profiles with the sublattice split
data/step3_metadata.json
figures/fig10 … fig14
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

import manybody as mb
import spectral as spx
import ssh_free as sf
from plot_style import CAT, LINESTYLES, MARKERS, REFERENCE, use_style, style_axes

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"

L_FULL = [8, 10, 12]                 # full spectrum in the addition sector
L_LANCZOS = [16, 18]
V_BASE = [0.0, 0.5, 1.0, 2.0, 3.0]
V_BIG = [0.0, 1.0, 3.0]              # L = 18 representative points
DELTAS = [-0.30, +0.30]
M_PROD = 200                         # production Lanczos vectors
M_TESTS = (100, 200, 400)            # the guide's convergence test
N_FIT = 8


def tag(L, d, V):
    return f"L{L}_d{d:+.2f}_V{V:.2f}".replace(".", "p")


# ---------------------------------------------------------------- one point
def analyse(spec: spx.AdditionSpectrum, n_fit: int = N_FIT) -> dict:
    """Everything Step 3 asks for, at one (L, delta, V)."""
    L = spec.L
    j_edge, j_centre = 0, L // 2
    wr = spx.window_report(spec, j_edge)
    sc = spx.window_scan(spec, j_edge)
    f_safe = spx.safe_factor(spec, j_edge)
    fit_e = spx.fit_weight_envelope(spec, 0, n_fit=n_fit)
    fit_safe = spx.fit_weight_envelope(spec, 0, factor=f_safe, n_fit=n_fit)
    fit_o = spx.fit_weight_envelope(spec, 1, n_fit=n_fit)
    fit_or = spx.fit_weight_envelope(spec, 1, n_fit=n_fit, from_right=True)
    wr_c = spx.window_report(spec, j_centre) if j_centre in spec.sites else {}

    row = {
        "L": L, "N": spec.N, "delta": spec.delta, "V": spec.V,
        "method": spec.method, "boundary": "open",
        "delta_plus_ring": spec.delta_plus_ring,
        "addition_gap_open": spec.addition_gap_open,
        "window_hi": wr["window_hi"], "window_status": sc["status"],
        "window_z_rel_spread": sc["z_rel_spread"],
        "n_lines_edge": wr["n_lines"],
        "W_0": wr["W_j"], "z_0": wr["z"], "n_eff_0": wr["n_eff"],
        "isolation_gap": wr["isolation_gap"],
        # the widest window in 0.6-0.9 that has not swallowed another line
        "factor_safe": f_safe,
        "z_0_safe": spec.z(j_edge, f_safe),
        "n_eff_safe": spec.n_eff(j_edge, f_safe),
        "n_lines_safe": spec.n_lines(j_edge, f_safe),
        "xi_W_safe": fit_safe.get("xi_W", np.nan),
        "xi_W_safe_r2": fit_safe.get("r2", np.nan),
        "window_choice_shifts_z": abs(spec.z(j_edge, f_safe) - wr["z"]),
        "one_minus_n0": 1.0 - spec.n_expect[0],
        "z_centre": wr_c.get("z", np.nan),
        "n_lines_centre": wr_c.get("n_lines", np.nan),
        "xi_W_even": fit_e.get("xi_W", np.nan),
        "xi_W_even_r2": fit_e.get("r2", np.nan),
        "xi_W_even_nfit": fit_e.get("n_fit", 0),
        "xi_W_even_site_max": fit_e.get("fit_site_max", -1),
        "xi_W_even_omitted_sub": fit_e.get("weight_omitted_frac_sublattice", np.nan),
        "xi_W_even_omitted_tot": fit_e.get("weight_omitted_frac_total", np.nan),
        "odd_decays_from_left": fit_o.get("ok", False),
        "xi_W_odd_from_right": fit_or.get("xi_W", np.nan),
        "W_total_all_sites": float(spec.weight_profile().sum()),
        "sum_rule_error": spec.sum_rule_error,
        "gs_parity_bare": spec.gs_parity, "gs_residual": spec.gs_residual,
        "gs_gap": spec.gs_gap, "gs_cluster": spec.gs_cluster,
        "fit_ok": fit_e.get("ok", False),
        "fit_reason": fit_e.get("reason", ""),
    }
    if spec.V == 0.0:
        ref = spx.free_reference(L, spec.delta)
        row["z0_free_reference"] = ref["z0"] if spec.delta < 0 else 0.0
        row["xi_free_reference"] = ref["xi"] if spec.delta < 0 else np.nan
    else:
        row["z0_free_reference"] = np.nan
        row["xi_free_reference"] = np.nan
    return row


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(
                f"{r[c]:.12g}" if isinstance(r[c], (int, float, np.floating))
                and not isinstance(r[c], bool) else str(r[c]) for c in cols) + "\n")


# ---------------------------------------------------------------- figures
def fig_spectrum(specs: dict) -> None:
    """Lines first, broadening only for presentation."""
    use_style()
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 6.2))
    L, d = 12, -0.30
    for col, V in enumerate([0.0, 1.0, 3.0]):
        s = specs[(L, d, V)]
        om, w = s.omega_of(0), s.w_of(0)
        hi = 0.8 * s.delta_plus_ring

        ax = axes[0, col]
        keep = w > spx.WEIGHT_FLOOR
        ax.axvspan(0, hi, color=CAT[2], alpha=0.10, lw=0)
        ax.vlines(om[keep], 0, w[keep], color=CAT[0], lw=1.4)
        inw = keep & (om < hi)
        ax.vlines(om[inw], 0, w[inw], color=CAT[1], lw=2.4)
        ax.set_xlim(-0.1, 3.2)
        ax.set_ylim(0, 0.42)
        ax.text(hi, 0.40, r"  $0.8\Delta_+$", fontsize=8, color="#55554f", va="top")
        style_axes(ax, title=rf"$V={V:g}$   $z_0={s.z(0):.3f}$,  "
                             rf"$n_{{\rm eff}}={s.n_eff(0):.2f}$",
                   ylabel=r"$w_{0m}$" if col == 0 else None)

        ax = axes[1, col]
        grid = np.linspace(-0.2, 3.2, 1400)
        for i, eta in enumerate(spx.ETAS):
            ax.plot(grid, spx.lorentzian(grid, om, w, eta), lw=1.5,
                    ls=LINESTYLES[i], color=CAT[i], label=rf"$\eta={eta}$")
        ax.axvspan(0, hi, color=CAT[2], alpha=0.10, lw=0)
        ax.set_xlim(-0.2, 3.2)
        style_axes(ax, xlabel=r"$\omega$",
                   ylabel=r"$A^+_0(\omega)$" if col == 0 else None)
        if col == 0:
            ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(rf"Addition spectrum at the edge site, $L={L}$, $\delta={d:+.2f}$"
                 "\ntop: the lines themselves (orange = inside the edge window)   ·   "
                 "bottom: the same lines broadened — $\\eta$ is a drawing choice, "
                 "not a lifetime",
                 x=0.012, y=0.995, ha="left", va="top", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig10_addition_spectrum")


def fig_edge_vs_centre(specs: dict, scan_rows: list[dict]) -> None:
    use_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.3))

    ax = axes[0]
    L, d, V = 12, -0.30, 1.0
    s = specs[(L, d, V)]
    hi = 0.8 * s.delta_plus_ring
    ax.axvspan(0, hi, color=CAT[2], alpha=0.10, lw=0)
    for i, (j, lbl) in enumerate(((0, "edge, $j=0$"), (L // 2, f"centre, $j={L//2}$"))):
        om, w = s.omega_of(j), s.w_of(j)
        keep = w > spx.WEIGHT_FLOOR
        ax.vlines(om[keep], 0, w[keep] * (1 if i == 0 else -1), color=CAT[i], lw=1.6,
                  label=f"{lbl}   $z={s.z(j):.3f}$")
    ax.axhline(0, color=REFERENCE, lw=0.8)
    ax.set_xlim(-0.1, 3.2)
    lim = max(abs(np.array(ax.get_ylim()))) * 1.08
    ax.set_ylim(-lim, lim)
    ax.text(hi, lim * 0.95, r"  $0.8\Delta_+$", fontsize=8,
            color="#55554f", va="top")
    style_axes(ax, title=rf"Edge vs. centre addition spectrum  "
                         rf"($L={L}$, $\delta={d:+.2f}$, $V={V:g}$)",
               xlabel=r"$\omega$", ylabel=r"$w_{jm}$   (centre plotted downwards)")
    ax.legend(loc="upper right", fontsize=8)

    ax = axes[1]
    for i, V in enumerate(V_BASE):
        rows = [r for r in scan_rows
                if r["L"] == 12 and r["delta"] == -0.30 and r["V"] == V]
        rows.sort(key=lambda r: r["factor"])
        ax.plot([r["factor"] for r in rows], [r["z"] for r in rows], lw=1.6,
                ls=LINESTYLES[i % len(LINESTYLES)], marker=MARKERS[i % len(MARKERS)],
                ms=5, color=CAT[i % len(CAT)], label=f"$V={V:g}$")
    ax.set_ylim(0, 0.85)
    ax.axvline(0.8, color=REFERENCE, lw=1.0, ls=":")
    ax.text(0.8, 0.03, "  guide's starting factor", fontsize=8, color="#55554f")
    style_axes(ax, title="Window stability: $z_0$ against the window factor\n"
                         r"flat ⇒ the boundary sits in a spectral gap wherever it is put",
               xlabel=r"factor $f$ in $0\leq\omega<f\Delta_+$", ylabel="$z_0$")
    ax.legend(loc="lower left", ncol=2, fontsize=8)

    # (c) the case where the default window is NOT isolated -- shown, not hidden
    ax = axes[2]
    key = (18, -0.30, 3.0)
    if key in specs:
        s = specs[key]
        om, w = s.omega_of(0), s.w_of(0)
        keep = w > spx.WEIGHT_FLOOR
        f_safe = spx.safe_factor(s, 0)
        ax.axvspan(0, f_safe * s.delta_plus_ring, color=CAT[2], alpha=0.14, lw=0)
        ax.vlines(om[keep], 0, w[keep], color=CAT[0], lw=1.5)
        inw = keep & (om < f_safe * s.delta_plus_ring)
        ax.vlines(om[inw], 0, w[inw], color=CAT[1], lw=2.6)
        for f, c, lbl in ((f_safe, CAT[2], rf"$f={f_safe:g}$ (isolated)"),
                          (0.8, CAT[1], r"$f=0.8$ (not isolated)")):
            ax.axvline(f * s.delta_plus_ring, color=c, lw=1.3, ls="--", label=lbl)
        ax.set_xlim(-0.1, 2.6)
        ax.set_ylim(0, 0.26)
        above = keep & (om >= f_safe * s.delta_plus_ring)
        x_int = float(om[above].min()) if above.any() else 1.18
        y_int = float(w[om == x_int][0]) if above.any() else 0.042
        ax.annotate("bottom of the open-chain\naddition continuum:\n"
                    r"inside $0.8\Delta_+$, outside $0.7\Delta_+$",
                    xy=(x_int, y_int), xytext=(1.45, 0.155), fontsize=7.8,
                    color="#55554f",
                    arrowprops=dict(arrowstyle="->", color="#8a8a82", lw=0.9))
        style_axes(ax, title=rf"Where the default window fails  "
                             rf"($L=18$, $\delta=-0.30$, $V=3$)" "\n"
                             rf"$\Delta_+$ outgrows the open chain's own continuum",
                   xlabel=r"$\omega$", ylabel=r"$w_{0m}$")
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig11_edge_vs_centre_and_window")


def fig_envelope(specs: dict) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    L, d = 18, -0.30

    ax = axes[0]
    s = specs[(L, d, 0.0)]
    W = s.weight_profile()
    sites = np.array(s.sites)
    ax.semilogy(sites[sites % 2 == 0], np.maximum(W[sites % 2 == 0], 1e-16), lw=1.5,
                marker="o", ms=5, color=CAT[0], label="even sites (left edge)")
    ax.semilogy(sites[sites % 2 == 1], np.maximum(W[sites % 2 == 1], 1e-16), lw=1.5,
                ls="--", marker="s", ms=5, mfc="none", color=CAT[1],
                label="odd sites (partner at the far edge)")
    f = spx.fit_weight_envelope(s, 0, n_fit=N_FIT)
    jj = np.linspace(0, f["fit_site_max"] + 1, 50)
    ax.semilogy(jj, np.exp(f["intercept"] + f["slope"] * jj), color="#1c1c1a",
                lw=1.0, ls=(0, (1, 2)), zorder=6,
                label=rf"fit $e^{{-2j/\xi_W}}$,  $\xi_W={f['xi_W']:.3f}$")
    ax.axvspan(f["fit_site_min"] - 0.4, f["fit_site_max"] + 0.4, color=CAT[2],
               alpha=0.10, lw=0)
    ax.set_ylim(1e-6, 2)
    style_axes(ax, title=rf"Window weight profile  ($L={L}$, $\delta={d:+.2f}$, $V=0$)"
                         "\n" rf"shaded: fit range;  exact $\xi={sf.xi_analytic(d):.4f}$",
               xlabel="site $j$", ylabel=r"$W_j=\sum_{m\in\mathcal{W}}w_{jm}$")
    ax.legend(loc="lower left", fontsize=7.8)

    ax = axes[1]
    for i, V in enumerate(V_BIG):
        s = specs[(L, d, V)]
        f = spx.safe_factor(s, 0)
        W = s.weight_profile(f)
        sites = np.array(s.sites)
        ev = sites % 2 == 0
        fit = spx.fit_weight_envelope(s, 0, factor=f, n_fit=N_FIT)
        ax.semilogy(sites[ev], np.maximum(W[ev], 1e-16), lw=1.5,
                    ls=LINESTYLES[i], marker=MARKERS[i], ms=5, color=CAT[i],
                    label=rf"$V={V:g}$  ($f={f:g}$),  $\xi_W={fit['xi_W']:.3f}$")
    ax.set_ylim(1e-6, 2)
    style_axes(ax, title="Even sublattice at three interaction strengths\n"
                         "each at its widest isolated window $f$",
               xlabel="site $j$", ylabel=r"$W_j$")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    save(fig, "fig12_weight_envelope")


def fig_xi_and_z(grid_rows: list[dict]) -> None:
    """The figure the guide asks for: xi_W(V) beside z_0(V)."""
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9))
    d = -0.30
    Ls = sorted({r["L"] for r in grid_rows if r["delta"] == d})

    def draw(ax, key_safe, key_default, ok_filter=None):
        """Plot the isolated-window value, and wherever the guide's default factor
        0.8 is NOT isolated, show where it would have landed instead."""
        for i, L in enumerate(Ls):
            rows = sorted([r for r in grid_rows
                           if r["L"] == L and r["delta"] == d
                           and (ok_filter is None or ok_filter(r))],
                          key=lambda r: r["V"])
            if not rows:
                continue
            col = CAT[i % len(CAT)]
            ax.plot([r["V"] for r in rows], [r[key_safe] for r in rows], lw=1.6,
                    ls=LINESTYLES[i % len(LINESTYLES)],
                    marker=MARKERS[i % len(MARKERS)], ms=5, color=col,
                    label=f"L = {L}")
            for r in rows:
                if r["factor_safe"] < 0.8 - 1e-9 and np.isfinite(r[key_default]):
                    ax.plot([r["V"], r["V"]], [r[key_safe], r[key_default]],
                            lw=1.1, color=col, alpha=0.8, zorder=2)
                    ax.plot([r["V"]], [r[key_default]], lw=0, marker="o", ms=8,
                            mfc="none", mew=1.4, color=col, zorder=3)

    ax = axes[0]
    ax.axhline(sf.xi_analytic(d), color=REFERENCE, lw=1.3,
               label=rf"Step 1  $\xi={sf.xi_analytic(d):.3f}$")
    draw(ax, "xi_W_safe", "xi_W_even", ok_filter=lambda r: r["fit_ok"])
    ax.set_ylim(2.4, 4.4)
    style_axes(ax, title=r"LOCALISATION:  $\xi_W(V)$" "\n"
                         "how far the edge excitation extends",
               xlabel="$V$", ylabel=r"$\xi_W$  (sites)")
    ax.legend(loc="upper left", ncol=2, fontsize=8)

    ax = axes[1]
    draw(ax, "z_0_safe", "z_0")
    r = (1 + d) / (1 - d)
    ax.plot([0], [1 - r**2], lw=0, marker="*", ms=13, color="#1c1c1a",
            label=rf"$V=0$ exact: $1-r^2={1-r**2:.3f}$")
    ax.set_ylim(0, 0.85)
    style_axes(ax, title=r"COHERENCE:  $z_0(V)$" "\n"
                         "how much of one added electron goes into it",
               xlabel="$V$", ylabel=r"$z_0(\mathcal{W})$")
    ax.legend(loc="lower left", ncol=2, fontsize=8)

    fig.suptitle(r"They answer different questions: at $\delta=-0.30$ the edge stays "
                 r"the same size while the electron creates it with falling probability"
                 "\nfilled: widest isolated window   ·   hollow, joined by a bar: "
                 r"where the guide's default $0.8\Delta_+$ lands when it is not isolated",
                 x=0.012, y=0.995, ha="left", va="top", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig13_xi_W_and_z0")


def fig_lanczos(conv_rows: list[dict]) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))

    ax = axes[0]
    for i, L in enumerate(sorted({r["L"] for r in conv_rows})):
        rows = sorted([r for r in conv_rows if r["L"] == L and r["V"] == 1.0],
                      key=lambda r: r["m"])
        ax.semilogy([r["m"] for r in rows],
                    [max(abs(r["z_0"] - r["z_0_reference"]), 2e-18)
                     if np.isfinite(r["z_0_reference"]) else np.nan for r in rows],
                    lw=1.6, ls=LINESTYLES[i], marker=MARKERS[i], ms=6, color=CAT[i],
                    label=f"L = {L}")
    ax.set_xticks(list(M_TESTS))
    ax.set_ylim(1e-18, 1e-2)
    style_axes(ax, title="Lanczos $z_0$ against the reference\n"
                         "$L=12$: full ED;  $L=16,18$: the 400-vector run",
               xlabel="number of Lanczos vectors", ylabel=r"$|z_0-z_0^{\rm ref}|$")
    ax.legend(loc="upper right")

    ax = axes[1]
    for i, L in enumerate(sorted({r["L"] for r in conv_rows})):
        rows = sorted([r for r in conv_rows if r["L"] == L and r["V"] == 1.0],
                      key=lambda r: r["m"])
        ax.plot([r["m"] for r in rows], [r["weight_in_lowest_10"] for r in rows],
                lw=1.6, ls=LINESTYLES[i], marker=MARKERS[i], ms=6, color=CAT[i],
                label=f"L = {L}")
    ax.axhline(1.0, color=REFERENCE, lw=1.2, ls="--")
    ax.text(0.97, 0.80, "total weight is 1 by construction;\n"
                        "the lowest few lines never carry it all",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            color="#55554f")
    ax.set_xticks(list(M_TESTS))
    ax.set_ylim(0, 1.15)
    style_axes(ax, title="Fraction of the sum rule held by the 10 lowest lines",
               xlabel="number of Lanczos vectors",
               ylabel=r"$\sum_{m<10}w_{0m}\,/\,(1-\langle n_0\rangle)$")
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, "fig14_lanczos_convergence")


def save(fig, name: str) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="full ED only")
    args = ap.parse_args(argv)
    DATA.mkdir(parents=True, exist_ok=True)

    print("Project A -- Step 3: localisation and coherence\n" + "=" * 72)
    specs: dict = {}
    grid: list[dict] = []
    scan_rows: list[dict] = []

    print("1. full spectrum in the addition sector (L <= 12)")
    for L in L_FULL:
        for d in DELTAS:
            for V in V_BASE:
                s = spx.addition_spectrum_full(L, d, V)
                specs[(L, d, V)] = s
                grid.append(analyse(s))
                for f in spx.FACTORS:
                    wr = spx.window_report(s, 0, f)
                    scan_rows.append({"L": L, "delta": d, "V": V, "method": s.method,
                                      **wr})
        print(f"   L = {L:2d}: dim(N+1) = {mb.Sector.make(L, L // 2 + 1).dim}, "
              f"max sum-rule error "
              f"{max(r['sum_rule_error'] for r in grid if r['L'] == L):.1e}")

    if not args.quick:
        print(f"\n2. Lanczos spectral measure at L = 16, 18 ({M_PROD} vectors)")
        for L in L_LANCZOS:
            vs = V_BASE if L == 16 else V_BIG
            sites = tuple(range(L // 2 + 1))
            for d in DELTAS:
                for V in vs:
                    s = spx.addition_spectrum_lanczos(L, d, V, M_PROD, sites=sites)
                    specs[(L, d, V)] = s
                    grid.append(analyse(s))
                    for f in spx.FACTORS:
                        scan_rows.append({"L": L, "delta": d, "V": V,
                                          "method": s.method,
                                          **spx.window_report(s, 0, f)})
            print(f"   L = {L:2d}: dim(N+1) = {mb.Sector.make(L, L // 2 + 1).dim}, "
                  f"{len(vs)} V values x 2 deltas, sites 0..{sites[-1]}")

    grid.sort(key=lambda r: (r["L"], r["delta"], r["V"]))
    write_csv(DATA / "step3_grid.csv", grid)
    write_csv(DATA / "step3_window_scan.csv", scan_rows)
    print(f"\n   wrote data/step3_grid.csv ({len(grid)} rows) and "
          f"step3_window_scan.csv ({len(scan_rows)} rows)")

    # per-point line lists and envelopes
    for (L, d, V), s in specs.items():
        if V not in V_BIG or d != -0.30:
            continue
        j_c = L // 2
        om0, w0 = s.omega_of(0), s.w_of(0)
        keep = w0 > spx.WEIGHT_FLOOR
        rows = [{"omega": float(o), "w_edge": float(w),
                 "w_centre": float(s.w_of(j_c)[i]) if s.omegas.ndim == 1 else np.nan,
                 "in_window": bool(o < 0.8 * s.delta_plus_ring and o >= 0)}
                for i, (o, w) in enumerate(zip(om0, w0)) if keep[i]]
        write_csv(DATA / f"step3_lines_{tag(L, d, V)}.csv", rows)
        W = s.weight_profile()
        write_csv(DATA / f"step3_envelope_{tag(L, d, V)}.csv",
                  [{"j": int(j), "W_j": float(W[i]),
                    "sublattice": int(j % 2)} for i, j in enumerate(s.sites)])

    print("\n3. Lanczos convergence test (100 / 200 / 400 vectors)")
    conv: list[dict] = []
    if not args.quick:
        for L in (12, 16, 18):
            for V in (0.0, 1.0, 3.0):
                ref = (specs[(12, -0.30, V)].z(0) if L == 12 else np.nan)
                entries = []
                for m in M_TESTS:
                    s = spx.addition_spectrum_lanczos(L, -0.30, V, m, sites=(0,))
                    om, w = s.omega_of(0), s.w_of(0)
                    order = np.argsort(om)
                    entries.append({
                        "L": L, "delta": -0.30, "V": V, "m": m,
                        "z_0": s.z(0), "n_eff_0": s.n_eff(0),
                        "n_lines": s.n_lines(0),
                        "sum_rule_error": s.sum_rule_error,
                        "weight_in_lowest_10":
                            float(w[order][:10].sum() / (1 - s.n_expect[0])),
                        "z_0_reference": ref,
                    })
                if not np.isfinite(ref):            # L = 16, 18: use the 400 run
                    for e in entries:
                        e["z_0_reference"] = entries[-1]["z_0"]
                conv.extend(entries)
        write_csv(DATA / "step3_lanczos_conv.csv", conv)
        worst = max(abs(e["z_0"] - e["z_0_reference"]) for e in conv
                    if e["L"] == 12)
        print(f"   L = 12 vs full ED: worst |z_0 - z_0_exact| = {worst:.1e}")
        print(f"   lowest 10 lines hold "
              f"{min(e['weight_in_lowest_10'] for e in conv):.3f}-"
              f"{max(e['weight_in_lowest_10'] for e in conv):.3f} of the sum rule "
              f"-- not all of it, as the guide warns")

    print("\n4. figures")
    fig_spectrum(specs)
    fig_edge_vs_centre(specs, scan_rows)
    if not args.quick:
        fig_envelope(specs)
    fig_xi_and_z(grid)
    if conv:
        fig_lanczos(conv)

    # ---- headline summary
    print("\n5. result at delta = -0.30 (weak end bonds)")
    print(f"   {'L':>3} {'V':>5} {'z_0':>8} {'n_eff':>7} {'xi_W':>8} {'f':>5} "
          f"{'window@0.8':>11}")
    for r in grid:
        if r["delta"] != -0.30 or r["L"] not in (12, 18):
            continue
        print(f"   {r['L']:>3} {r['V']:>5.2f} {r['z_0_safe']:>8.4f} "
              f"{r['n_eff_safe']:>7.3f} {r['xi_W_safe']:>8.4f} "
              f"{r['factor_safe']:>5.1f} {r['window_status']:>11}")
    print("   (z_0, n_eff and xi_W quoted at the widest isolated window f; "
          "column 'window@0.8' is the status at the guide's default factor)")
    ctl = [r for r in grid if r["delta"] > 0]
    print(f"   control, delta = +0.30 (strong end bonds): window status "
          f"{ {r['window_status'] for r in ctl} }, z_0 = "
          f"{max(r['z_0'] for r in ctl):.1e} at every point")

    meta = {
        "project": "A -- can a topological edge survive strong interactions?",
        "step": 3,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "scipy": scipy.__version__, "platform": platform.platform(),
        "command": " ".join([Path(sys.argv[0]).name] + sys.argv[1:]),
        "conventions": {
            "inherited_from": "Steps 1 and 2 (ssh_free.py, manybody.py)",
            "w_jm": "|<m,N+1|c^dag_j|Omega_N>|^2",
            "omega_m": "E_m^{N+1} - E_0^N",
            "sum_rule": "sum_m w_jm = 1 - <n_j> (= 1/2 at half filling)",
            "z": "cluster weight sum_{m in W} w_jm / (1 - <n_j>)",
            "edge_window": "0 <= omega < f * Delta_plus, f = 0.8 by default, "
                           "scanned over 0.6-0.9; Delta_plus measured on a ring "
                           "at the same parameters",
            "xi_W": "slope of ln W_j vs j is -2/xi_W (W_j is a weight)",
            "broadening": "Lorentzian eta in (0.02, 0.01, 0.005) is PRESENTATION "
                          "ONLY; every reported number comes from the lines",
            "state_selection": "parity+ (Step 2 rule)",
            "random_seed": None,
        },
        "grid": {"L_full": L_FULL, "L_lanczos": [] if args.quick else L_LANCZOS,
                 "delta": DELTAS, "V_base": V_BASE, "V_big": V_BIG,
                 "lanczos_production": M_PROD, "lanczos_tests": list(M_TESTS),
                 "n_fit": N_FIT},
        "worst": {
            "sum_rule_error": max(r["sum_rule_error"] for r in grid),
            "gs_residual": max(r["gs_residual"] for r in grid),
            "window_z_rel_spread": max(r["window_z_rel_spread"] for r in grid),
        },
        "caveats": [
            "A finite system always has discrete lines; a Lorentzian width is a "
            "drawing choice and never a lifetime.",
            "z_j is a normalised electron-addition weight in a stated window, not "
            "the relativistic LSZ field-renormalisation constant.",
            "The Lanczos total weight equals the sum rule by construction; the "
            "lowest few lines do not, and no full sum rule is claimed from them.",
            "Only size dependence can say whether the dominant line survives or "
            "fragments in a larger system; L <= 18 here.",
        ],
    }
    (DATA / "step3_metadata.json").write_text(json.dumps(meta, indent=2))
    print("\n  wrote data/step3_metadata.json")
    print("\nStep 3 complete.  Run `python test_step3.py` for the gate checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
