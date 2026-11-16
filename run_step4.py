#!/usr/bin/env python3
"""
Project A -- Step 4 driver: watch an edge particle evolve.

    python run_step4.py            # full programme (~12 min)
    python run_step4.py --quick    # L = 8, 12, 16 only; no L = 18/20 ladder

Order of work follows the guide:
  1. Krylov evolution at L = 8, 12, 16, dt = 0.05, 0 <= t <= 8; full
     diagonalisation at L = 8 (and 12); tolerance halved everywhere.
  2. Controls: injection near the centre, and the opposite dimerisation.
  3. Norm, energy, sum_j Delta n_j = 1; density profiles, not only Q_edge.
  4. Boundary reflections and the two-edge beating; the finite-size-free interval
     is certified by comparing each L with the next larger one, and the time
     range is extended to t = 16 only at L = 18, 20.
  5. Filtered preparation onto the edge window, with z_0(W).

Outputs
-------
data/step4_runs.csv               one row per run: checks, windows, late means
data/step4_finite_size.csv        certified finite-size-free intervals
data/step4_tolerance_scan.csv     Krylov error against tolerance
data/step4_beating.csv            splitting-induced beating inspection
data/step4_far_edge.csv           t = 0 far-edge response, selected vs broken
data/step4_trajectories.npz       every trajectory: t, G(t), Delta n_j(t), norm, E
data/step4_traces_*.csv           representative time series
data/step4_metadata.json
figures/fig15 … fig19
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time as _time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

import dynamics as dy
import manybody as mb
import spectral as spx
import ssh_free as sf
from plot_style import (CAT, LINESTYLES, MARKERS, REFERENCE, INK, INK_MUTED,
                        use_style, style_axes)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"

L_MAIN = [8, 12, 16]
L_LADDER = [18, 20]
DELTAS = [-0.30, +0.30]
V_BASE = [0.0, 0.5, 1.0, 2.0, 3.0]
V_BIG = [0.0, 1.0, 3.0]
T_EXT = 16.0
EPS_VISIBLE = 0.01           # "visible" departure between sizes
EPS_TIGHT = 1e-3
TOLS_SCAN = (1e-3, 1e-5, 1e-7, 1e-9, 1e-10, 5e-11)
LANCZOS_SPEC_M = 400


def centre_site(L: int) -> int:
    """Even site nearest the centre, so every L injects into the same local
    environment (A sublattice, strong bond to the left, weak to the right)."""
    return 2 * (L // 4)


def key(L, d, V, site, kind="electron"):
    return (L, round(d, 2), round(V, 2), "edge" if site == 0 else "centre", kind)


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- one run
def run_one(L, d, V, site, t_max=dy.T_MAX, halve=True, exact=None,
            spectral=None):
    """Krylov run plus every check the guide asks for at this point."""
    exact = (L <= 12) if exact is None else exact
    prep = dy.prepare(L, d, V, site)
    k = dy.krylov_evolve(prep, t_max=t_max)
    row = {
        "L": L, "delta": d, "V": V, "site": site,
        "injection": "edge" if site == 0 else "centre",
        "kind": "electron", "t_max": t_max, "tol": k.tol,
        "dim_Np1": prep.H.shape[0], "prep_weight": prep.prep_weight,
        "gs_parity_bare": prep.gs_parity, "gs_residual": prep.gs_residual,
        "neutral_gap_N": prep.gs_gap,
        "krylov_dim_max": k.krylov_dim_max, "krylov_vectors": k.krylov_vectors_total,
        "krylov_splits": k.n_splits, "krylov_err_estimate": k.err_estimate,
        "norm_drift": float(np.max(np.abs(k.norm - 1.0))),
        "energy_drift": float(np.max(np.abs(k.energy - k.energy[0]))),
        "sum_dn_minus_1": float(np.max(np.abs(k.total_excess() - 1.0))),
    }
    # exact t = 0 identities
    pred = -prep.corr_ref / (1.0 - prep.n_ref[site])
    pred[site] = 1.0 - prep.n_ref[site]
    row["t0_identity_err"] = float(np.max(np.abs(k.dn[0] - pred)))

    if halve:
        kh = dy.krylov_evolve(prep, t_max=t_max, tol=k.tol / 2)
        c = dy.compare(k, kh)
        row.update({"halving_amp": c["amp"], "halving_P": c["P_return"],
                    "halving_dn": c["dn"],
                    "halving_extra_vectors": kh.krylov_vectors_total
                    - k.krylov_vectors_total})
    else:
        row.update({"halving_amp": np.nan, "halving_P": np.nan,
                    "halving_dn": np.nan, "halving_extra_vectors": np.nan})

    if exact:
        e = dy.exact_evolve(prep, t_max=t_max)
        c = dy.compare(k, e)
        row.update({"exact_amp": c["amp"], "exact_P": c["P_return"],
                    "exact_dn": c["dn"]})
    else:
        row.update({"exact_amp": np.nan, "exact_P": np.nan, "exact_dn": np.nan})

    if spectral is not None:
        G, _ = dy.spectral_amplitude(L, d, V, site, k.times,
                                     lanczos_m=None if spectral == "full"
                                     else LANCZOS_SPEC_M)
        row.update({"spectral_method": spectral,
                    "spectral_amp": float(np.max(np.abs(k.amp - G))),
                    "spectral_P": float(np.max(np.abs(k.P_return - np.abs(G) ** 2)))})
    else:
        row.update({"spectral_method": "", "spectral_amp": np.nan,
                    "spectral_P": np.nan})

    if V == 0.0:
        Gf, dnf = dy.free_fermion_trajectory(L, d, site, k.times)
        row.update({"free_fermion_amp": float(np.max(np.abs(k.amp - Gf))),
                    "free_fermion_dn": float(np.max(np.abs(k.dn - dnf)))})
    else:
        row.update({"free_fermion_amp": np.nan, "free_fermion_dn": np.nan})

    row.update({
        "Q2_0": float(k.Q_local(2)[0]), "Q4_0": float(k.Q_local(4)[0]),
        "Qright4_0": float(k.Q_right(4)[0]), "dn_far_0": float(k.dn[0, -1]),
    })
    return row, k, prep


def run_filtered(prep, t_max=dy.T_MAX, evolve=True):
    f = dy.prepare_filtered(prep)
    if f is None:
        return {"L": prep.L, "delta": prep.delta, "V": prep.V, "site": prep.site,
                "injection": "edge", "kind": "filtered", "window_empty": True,
                "z0W": 0.0}, None, None
    kf = dy.krylov_evolve(f, t_max=t_max if evolve else 1.0)
    row = {
        "L": prep.L, "delta": prep.delta, "V": prep.V, "site": prep.site,
        "injection": "edge", "kind": "filtered", "window_empty": False,
        "z0W": f.prep_weight, "z0W_step3": float(f.window["z_step3"]),
        "window_factor": f.window["factor"], "window_lines": f.window["n_lines"],
        "P_return_min": float(kf.P_return.min()),
        "dn_static_dev": float(np.max(np.abs(kf.dn - kf.dn[0]))),
        "q2": float(kf.Q_left(2)[0]), "q4": float(kf.Q_left(4)[0]),
        "q_left_half": float(kf.dn[0, : prep.L // 2].sum()),
        "q_right_half": float(kf.dn[0, prep.L // 2:].sum()),
        "norm_drift": float(np.max(np.abs(kf.norm - 1.0))),
    }
    return row, kf, f


# ---------------------------------------------------------------- analysis
def certify(trajs: dict) -> list[dict]:
    """Finite-size-free interval for every run that has a larger partner."""
    rows = []
    ladder = sorted({k[0] for k in trajs})
    for (L, d, V, inj, kind), tr in trajs.items():
        if kind != "electron":
            continue
        bigger = [Lb for Lb in ladder if Lb > L and (Lb, d, V, inj, kind) in trajs]
        if not bigger:
            continue
        Lb = bigger[0]
        big = trajs[(Lb, d, V, inj, kind)]
        for eps in (EPS_VISIBLE, EPS_TIGHT):
            r = dy.finite_size_time(tr, big, eps=eps)
            rows.append({"L": L, "L_ref": Lb, "delta": d, "V": V, "injection": inj,
                         "eps": eps, "t_fs": r["t_fs"], "reached": r["reached"],
                         "t_common_max": r["t_common_max"],
                         "first_by_P": r["by"]["P_return"],
                         "first_by_Q2": r["by"]["Q2"], "first_by_Q4": r["by"]["Q4"],
                         "static_offset_Q2": r["static_offset"]["Q2"],
                         "static_offset_Q4": r["static_offset"]["Q4"],
                         "static_offset_max": r["static_offset_max"],
                         "initial_state_ok": r["initial_state_ok"]})
    return rows


def late_window(t_fs: float, t_max: float, reached: bool):
    """Late-time averaging interval inside the certified window: the second half
    of [0, min(t_fs, t_max)], starting no earlier than t = 2."""
    t_b = min(t_fs, t_max) if reached else t_max
    t_a = max(2.0, 0.5 * t_b)
    return (t_a, t_b) if t_b - t_a >= 1.0 else (np.nan, np.nan)


# ---------------------------------------------------------------- figures
DIVERGING = LinearSegmentedColormap.from_list(
    "projA_div", ["#0072B2", "#e7e7e3", "#D55E00"])


def save(fig, name):
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    log(f"  wrote figures/{name}.png / .pdf")


def fig_validation(trajs, scan_rows, runs):
    use_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.3))

    ax = axes[0]
    for i, V in enumerate(V_BIG):
        tr = trajs[key(8, -0.30, V, 0)]
        e = dy.exact_evolve(tr.prep)
        G, _ = dy.spectral_amplitude(8, -0.30, V, 0, tr.times)
        ax.semilogy(tr.times, np.maximum(np.abs(tr.amp - e.amp), 1e-17), lw=1.4,
                    color=CAT[i], label=rf"$V={V:g}$: Krylov − exact")
        ax.semilogy(tr.times, np.maximum(np.abs(tr.amp - G), 1e-17), lw=1.2,
                    ls="--", color=CAT[i], label=rf"$V={V:g}$: Krylov − spectral")
    ax.axhline(dy.TOL * dy.T_MAX, color=REFERENCE, lw=1.1, ls=":")
    ax.text(0.2, dy.TOL * dy.T_MAX * 1.6, r"accumulated bound tol$\cdot t_{\max}$",
            fontsize=8, color="#55554f")
    ax.set_ylim(1e-17, 1e-7)
    style_axes(ax, title=r"$L=8$: complex amplitude $G(t)$, three routes",
               xlabel="$t$", ylabel=r"$|G_{\rm Krylov}-G|$")
    ax.legend(loc="center right", fontsize=7.2)

    ax = axes[1]
    for i, L in enumerate(sorted({r["L"] for r in scan_rows})):
        rr = sorted([r for r in scan_rows if r["L"] == L], key=lambda r: r["tol"])
        ax.loglog([r["tol"] for r in rr], [max(r["err_amp"], 1e-17) for r in rr],
                  lw=1.6, marker=MARKERS[i], ms=6, color=CAT[i],
                  label=f"L = {L} vs {rr[0]['reference']}")
    tt = np.array(TOLS_SCAN)
    ax.loglog(tt, tt * dy.T_MAX, color=REFERENCE, lw=1.2, ls="--",
              label=r"bound tol$\cdot t_{\max}$")
    ax.axvline(dy.TOL, color=REFERENCE, lw=0.9, ls=":")
    ax.text(dy.TOL * 1.3, 1e-2, "production", fontsize=8, color="#55554f",
            rotation=90, va="top")
    ax.set_ylim(1e-16, 1e-1)
    style_axes(ax, title="Tolerance scan: error falls with the tolerance\n"
                         "and stays ~10³ below the a-posteriori bound",
               xlabel="Krylov tolerance (error per unit time)",
               ylabel=r"$\max_t |G-G_{\rm ref}|$")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[2]
    tr = trajs[key(16, -0.30, 3.0, 0)]
    for i, (lbl, y) in enumerate((
            (r"$|\langle\Psi|\Psi\rangle-1|$", np.abs(tr.norm - 1)),
            (r"$|E(t)-E(0)|$", np.abs(tr.energy - tr.energy[0])),
            (r"$|\sum_j\Delta n_j-1|$", np.abs(tr.total_excess() - 1)))):
        ax.semilogy(tr.times, np.maximum(y, 1e-17), lw=1.4, ls=LINESTYLES[i],
                    color=CAT[i], label=lbl)
    tr0 = trajs[key(16, -0.30, 0.0, 0)]
    _, dnf = dy.free_fermion_trajectory(16, -0.30, 0, tr0.times)
    ax.semilogy(tr0.times, np.maximum(np.abs(tr0.dn - dnf).max(1), 1e-17), lw=1.2,
                color=CAT[3], ls=LINESTYLES[3],
                label=r"$V=0$: $\max_j|\Delta n_j-|\chi_j|^2|$ (one-body route)")
    ax.set_ylim(1e-17, 1e-8)
    style_axes(ax, title=r"$L=16$: conservation ($V=3$) and the" "\n"
                         r"independent free-fermion route ($V=0$)",
               xlabel="$t$", ylabel="deviation")
    ax.legend(loc="upper left", fontsize=7.6)
    fig.tight_layout()
    save(fig, "fig15_dynamics_validation")


def _cells(tr):
    """Excess charge per two-site cell: removes the sublattice stripes."""
    return tr.dn[:, 0::2] + tr.dn[:, 1::2]


TIME_RAMP = ["#9cc5e3", "#5fa2d1", "#0072B2", "#00548a", "#003660"]   # one hue, ordered


def fig_spacetime(trajs, fs_lookup, v_free, speeds):
    """Top: where the charge GOES -- change since t = 0, per cell.
    Bottom: the density profile itself, Delta n_j(t), at a few times."""
    from matplotlib.colors import TwoSlopeNorm
    use_style()
    fig = plt.figure(figsize=(14.4, 8.6))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.045],
                          height_ratios=[1.2, 1.0], hspace=0.42, wspace=0.2)

    tops = []
    for V in V_BIG:
        tr = trajs.get(key(20, -0.30, V, 0)) or trajs[key(16, -0.30, V, 0)]
        ch = _cells(tr) - _cells(tr)[0]
        tops.append((V, tr, ch))
    # scale set by the bulk cells: the injection cell's own depletion would
    # otherwise wash out the propagating front (it saturates, which is fine)
    vmax = float(np.percentile(np.abs(np.concatenate([c[:, 1:].ravel()
                                                       for *_, c in tops])), 99))
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    im = None
    for col, (V, tr, ch) in enumerate(tops):
        ax = fig.add_subplot(gs[0, col])
        L = tr.prep.L
        x_edges = 2.0 * np.arange(L // 2 + 1) - 0.5
        t_edges = np.append(tr.times, tr.times[-1] + dy.DT_OUT) - 0.5 * dy.DT_OUT
        im = ax.pcolormesh(x_edges, t_edges, ch, cmap=DIVERGING, norm=norm,
                           shading="flat", rasterized=True)
        tt = np.linspace(0, tr.times[-1], 300)
        sp = speeds.get(V)
        v_line = sp["v"] if (sp is not None and not sp["bound"]) else v_free
        out = v_line * tt
        ax.plot(np.where(out <= L - 1, out, np.nan), tt, color=INK, lw=0.9,
                ls=(0, (3, 2)), alpha=0.75)
        back = 2 * (L - 1) - out
        ax.plot(np.where((out > L - 1) & (back >= 0), back, np.nan), tt, color=INK,
                lw=0.9, ls=(0, (3, 2)), alpha=0.75)
        t_fs = fs_lookup.get((L, -0.30, V, "edge"))
        if t_fs is not None and np.isfinite(t_fs) and t_fs < tr.times[-1]:
            ax.axhline(t_fs, color=INK, lw=1.4)
            ax.text(L - 0.8, t_fs, f"certified to t = {t_fs:.1f}", color=INK,
                    fontsize=7.6, ha="right", va="bottom")
        ax.set_xlim(-0.5, L - 0.5)
        ax.set_ylim(0, tr.times[-1])
        ax.grid(False)
        for sd in ("top", "right"):
            ax.spines[sd].set_visible(False)
        ax.set_title(rf"edge injection, $V={V:g}$, $L={L}$;  dashed $v={v_line:.1f}$",
                     loc="left", fontsize=9.6)
        ax.set_xlabel("site $j$ (two-site cells)")
        if col == 0:
            ax.set_ylabel("time $t$")
    cax = fig.add_subplot(gs[0, 3])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$\Delta N_{\rm cell}(t)-\Delta N_{\rm cell}(0)$"
                 "\n(injection cell saturates)")

    bottoms = (
        (key(16, -0.30, 1.0, 0), r"edge injection, $\delta=-0.30$"),
        (key(16, -0.30, 1.0, centre_site(16)), r"centre injection, $\delta=-0.30$"),
        (key(16, +0.30, 1.0, 0), r"edge injection, $\delta=+0.30$ (no edge state)"),
    )
    snaps = (0.0, 1.0, 2.0, 4.0, 8.0)
    for col, (k_, title) in enumerate(bottoms):
        ax = fig.add_subplot(gs[1, col])
        tr = trajs[k_]
        L = tr.prep.L
        for i, t in enumerate(snaps):
            idx = int(round(t / dy.DT_OUT))
            ax.plot(np.arange(L), tr.dn[idx], lw=1.3, marker=MARKERS[i % 5], ms=3.6,
                    color=TIME_RAMP[i], label=f"t = {t:g}")
        if col == 0 and key(16, -0.30, 1.0, 0, "filtered") in trajs:
            f = trajs[key(16, -0.30, 1.0, 0, "filtered")]
            ax.plot(np.arange(L), f.dn[0], color=CAT[1], lw=1.4, ls="--",
                    label="filtered (static)")
        ax.axhline(0, color=REFERENCE, lw=0.8)
        ax.set_ylim(-0.2, 0.55)
        style_axes(ax, title=title + r", $V=1$, $L=16$", xlabel="site $j$",
                   ylabel=r"$\Delta n_j(t)$" if col == 0 else None)
        if col == 0:
            ax.legend(loc="upper center", fontsize=7.4, ncol=3)
    fig.suptitle(
        "Where the added electron goes.  Top: change of the excess charge per two-site "
        "cell since $t=0$ (largest $L$ available).  Dashed: a front and its reflection at "
        rf"the free maximal group velocity $2\min(J_w,J_s)={v_free:.1f}$ ($V=0$) or at the "
        "speed implied by the certified return ($V>0$)." "\n"
        "Bottom: the density profile itself — note the excess already sitting at the far "
        r"end at $t=0$, and the negative $\Delta n_j$ where the background rearranges.",
        x=0.012, y=0.995, ha="left", va="top", fontsize=10)
    fig.subplots_adjust(left=0.05, right=0.94, top=0.88, bottom=0.07)
    save(fig, "fig16_spacetime_density")


def fig_traces(trajs, filt_rows, fs_lookup, z0):
    use_style()
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 7.6), sharex=True)
    for col, V in enumerate(V_BIG):
        e16 = trajs[key(16, -0.30, V, 0)]
        c16 = trajs[key(16, -0.30, V, centre_site(16))]
        t16 = trajs[key(16, +0.30, V, 0)]
        f16 = trajs.get(key(16, -0.30, V, 0, "filtered"))
        e20 = trajs.get(key(20, -0.30, V, 0))
        fr = next(r for r in filt_rows if r["L"] == 16 and r["delta"] == -0.30
                  and r["V"] == V)

        ax = axes[0, col]
        if e20 is not None:
            ax.plot(e20.times, e20.P_return, color=CAT[0], lw=0.9, alpha=0.45,
                    label="edge, $L=20$ (to $t=16$)")
        ax.plot(e16.times, e16.P_return, color=CAT[0], lw=1.8, label="edge, $L=16$")
        ax.plot(c16.times, c16.P_return, color=CAT[1], lw=1.5, ls="--",
                label="centre injection")
        ax.plot(t16.times, t16.P_return, color=CAT[2], lw=1.5, ls="-.",
                label=r"edge, $\delta=+0.30$")
        if f16 is not None:
            ax.plot(f16.times, f16.P_return, color=CAT[4], lw=1.2, ls=":",
                    label="filtered (eigenstate)")
        ax.axhline(fr["z0W"] ** 2, color=REFERENCE, lw=1.1, ls=(0, (4, 2)))
        ax.text(T_EXT * 0.98, fr["z0W"] ** 2, rf"$z_0(\mathcal{{W}})^2={fr['z0W']**2:.3f}$",
                fontsize=8, color="#55554f", ha="right", va="bottom")
        ax.set_ylim(-0.03, 1.05)
        style_axes(ax, title=rf"$V={V:g}$:  return probability",
                   ylabel=r"$P_{\rm return}(t)$" if col == 0 else None)

        ax = axes[1, col]
        if e20 is not None:
            ax.plot(e20.times, e20.Q_left(4), color=CAT[0], lw=0.9, alpha=0.45)
        ax.plot(e16.times, e16.Q_left(4), color=CAT[0], lw=1.8,
                label=r"edge, $Q_{\rm edge}$ ($\ell_e=4$)")
        ax.plot(e16.times, e16.Q_left(2), color=CAT[0], lw=1.0, ls=(0, (1, 1)),
                label=r"edge, $Q_{\rm edge}$ ($\ell_e=2$)")
        ax.plot(e16.times, e16.Q_right(4), color=CAT[5], lw=1.4, ls=(0, (5, 1)),
                label="far edge, last 4 sites")
        ax.plot(c16.times, c16.Q_local(4), color=CAT[1], lw=1.5, ls="--",
                label="centre injection, 4 central sites")
        ax.plot(t16.times, t16.Q_left(4), color=CAT[2], lw=1.5, ls="-.",
                label=r"edge, $\delta=+0.30$")
        ax.axhline(fr["z0W"] * fr["q4"], color=REFERENCE, lw=1.1, ls=(0, (4, 2)))
        ax.text(T_EXT * 0.98, fr["z0W"] * fr["q4"],
                rf"$z_0\,q_4={fr['z0W']*fr['q4']:.3f}$", fontsize=8,
                color="#55554f", ha="right", va="top")
        ax.set_ylim(-0.05, 1.0)
        style_axes(ax, xlabel="$t$",
                   ylabel="excess charge" if col == 0 else None)
        t16 = fs_lookup.get((16, -0.30, V, "edge"), np.nan)
        t20 = fs_lookup.get((20, -0.30, V, "edge"), np.nan)
        t_end = T_EXT if e20 is not None else dy.T_MAX
        for a in (axes[0, col], axes[1, col]):
            if np.isfinite(t16) and t16 < dy.T_MAX:
                a.axvspan(t16, dy.T_MAX, color=INK_MUTED, alpha=0.16, lw=0)
            if e20 is not None and np.isfinite(t20) and t20 < T_EXT:
                a.axvspan(max(t20, dy.T_MAX), T_EXT, color=INK_MUTED, alpha=0.08, lw=0)
            a.set_xlim(0, t_end)
    axes[0, 0].legend(loc="upper right", fontsize=7.6)
    axes[1, 0].legend(loc="upper right", fontsize=7.3)
    cert = []
    for V in V_BIG:
        t16 = fs_lookup.get((16, -0.30, V, "edge"), np.nan)
        cert.append(f"V={V:g}: " + ("not certified" if not np.isfinite(t16) and
                                     (16, -0.30, V, "edge") not in fs_lookup
                                     else (f"t<{t16:.1f}" if t16 < dy.T_MAX
                                           else "all t<=8")))
    fig.suptitle(r"Time traces, $\delta=-0.30$ unless stated.  "
                 r"Grey dashed: what the edge eigenstate alone predicts "
                 r"($z_0^2$ for the return, $z_0 q_4$ for the edge charge)." "\n"
                 r"$L=16$ certified free of finite-size returns (vs $L=18$): "
                 + ", ".join(cert) + ".  Shaded: outside the certified interval; "
                 r"faint curves: $L=20$.",
                 x=0.012, y=0.995, ha="left", va="top", fontsize=9.8)
    fig.tight_layout(rect=[0, 0, 1, 0.935])
    save(fig, "fig17_time_traces")


def fig_finite_size(trajs, fs_rows, beat, speeds):
    use_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.4))

    ax = axes[0]
    pairs = [(8, 12), (12, 16), (16, 18), (18, 20)]
    for i, (a, b) in enumerate(pairs):
        ka, kb = key(a, -0.30, 1.0, 0), key(b, -0.30, 1.0, 0)
        if ka not in trajs or kb not in trajs:
            continue
        A, B = trajs[ka], trajs[kb]
        n = min(A.times.size, B.times.size)
        qa, qb = A.Q_left(4)[:n], B.Q_left(4)[:n]
        d = np.maximum(np.abs(A.P_return[:n] - B.P_return[:n]),
                       np.abs((qa - qa[0]) - (qb - qb[0])))
        ax.semilogy(A.times[:n], np.maximum(d, 1e-12), lw=1.5, ls=LINESTYLES[i],
                    color=CAT[i], label=f"L = {a} vs {b}")
    ax.axhline(EPS_VISIBLE, color=INK, lw=1.1)
    ax.text(0.2, EPS_VISIBLE * 1.4, r"visible: $\varepsilon=0.01$", fontsize=8,
            color=INK)
    ax.axhline(EPS_TIGHT, color=REFERENCE, lw=0.9, ls=":")
    ax.set_ylim(1e-10, 1)
    style_axes(ax, title=r"Departure between sizes ($\delta=-0.30$, $V=1$)" "\n"
                         r"$|\Delta P_{\rm return}|$ or the change of $Q_4$ since $t=0$",
               xlabel="$t$", ylabel="difference")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    from matplotlib.lines import Line2D
    for i, V in enumerate(V_BIG):
        rr = sorted([r for r in fs_rows if r["delta"] == -0.30 and r["V"] == V
                     and r["injection"] == "edge" and r["eps"] == EPS_VISIBLE],
                    key=lambda r: r["L"])
        dx = (i - 1) * 0.28                      # keep coincident points apart
        for r in rr:
            x = r["L"] + dx
            y = r["t_fs"] if r["reached"] else r["t_common_max"]
            if r.get("extrapolated"):
                ax.plot([x], [min(y, 17.6)], lw=0, marker="o" if r["reached"] else 6,
                        ms=9, mfc="none", mew=1.5, color=CAT[i])
            elif not r["initial_state_ok"]:
                ax.plot([x], [y], lw=0, marker="x", ms=8, mew=1.6, color=CAT[i])
            elif r["reached"]:
                ax.plot([x], [y], lw=0, marker="o", ms=8, color=CAT[i])
            else:
                ax.plot([x], [y], lw=0, marker=6, ms=11, color=CAT[i])
        if V in speeds and not speeds[V]["bound"]:
            xx = np.linspace(8, 20.5, 60)
            ax.plot(xx, (2 * (xx - 1) - 3) / speeds[V]["v"], color=CAT[i], lw=1.0,
                    ls=":")
    handles = [Line2D([], [], color=CAT[i], marker="o", lw=0, ms=7, label=f"$V={V:g}$")
               for i, V in enumerate(V_BIG)]
    handles += [Line2D([], [], color=INK_MUTED, marker=m, lw=0, ms=8, mfc=fc, mew=1.4,
                       label=lab) for m, fc, lab in (
                    ("o", None, "return seen"),
                    (6, None, "no return by end of run"),
                    ("x", None, "initial state not converged"),
                    ("o", "none", "L = 20, scaled from 18"))]
    ax.legend(handles=handles, loc="upper left", fontsize=7.4, ncol=2)
    ax.set_xticks([8, 12, 16, 18, 20])
    ax.set_ylim(0, 19.5)
    style_axes(ax, title=r"Return time vs $L$ (edge, $\delta=-0.30$)" "\n"
                         "dotted: round trip at the implied speed",
               xlabel="$L$", ylabel=r"$t_{\rm fs}$")

    ax = axes[2]
    for i, r in enumerate(beat):
        if r["V"] not in (0.0, 3.0):
            continue
        ax.semilogy(r["freqs"], np.maximum(r["spectrum_hann"], 1e-7), lw=1.2,
                    color=CAT[i], ls=LINESTYLES[i], label=f"$V={r['V']:g}$")
        ax.axvline(r["splitting"], color=CAT[i], lw=1.1, ls=":")
    ax.annotate("dotted: N-sector splitting $E_1-E_0$,\nwhere a two-edge beat "
                "would sit",
                xy=(beat[0]["splitting"], 0.3), xytext=(1.95, 0.25),
                fontsize=8, color="#55554f",
                arrowprops=dict(arrowstyle="->", color="#8a8a82", lw=0.9))
    ax.set_xlim(0, 3.0)
    ax.set_ylim(1e-7, 1)
    style_axes(ax, title=r"Beating check, $L=8$, exact to $t=300$" "\n"
                         r"spectrum of $Q_{\rm edge}(t)$ (Hann window)",
               xlabel=r"frequency $\nu$", ylabel=r"$|\widetilde{Q}_4(\nu)|$")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig18_finite_size_and_beating")


def fig_summary(runs, filt_rows, far_rows, bound):
    use_style()
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.5))

    def late(L, d, inj, col):
        rr = sorted([r for r in runs if r["L"] == L and r["delta"] == d
                     and r["injection"] == inj and r["kind"] == "electron"
                     and np.isfinite(r.get(col, np.nan))], key=lambda r: r["V"])
        return np.array([r["V"] for r in rr]), np.array([r[col] for r in rr])

    fz = sorted([r for r in filt_rows if r["L"] == 16 and r["delta"] == -0.30
                 and not r["window_empty"]], key=lambda r: r["V"])
    Vz = np.array([r["V"] for r in fz])
    z0 = np.array([r["z0W"] for r in fz])
    q4 = np.array([r["q4"] for r in fz])

    ax = axes[0]
    ax.plot(Vz, z0**2, color=INK, lw=1.4, ls="--", marker="*", ms=9,
            label=r"$z_0(\mathcal{W})^2$  (Step 3)")
    for i, (L, inj, d, lbl) in enumerate((
            (16, "edge", -0.30, "edge, L = 16"),
            (20, "edge", -0.30, "edge, L = 20 (window scaled)"),
            (16, "centre", -0.30, "centre, L = 16"),
            (16, "edge", +0.30, r"edge $\delta=+0.30$, L = 16"))):
        V, y = late(L, d, inj, "P_late")
        if V.size:
            ax.plot(V, y, lw=1.5, ls=LINESTYLES[i], marker=MARKERS[i], ms=6,
                    color=CAT[i], label=lbl)
    ax.set_ylim(-0.02, 0.62)
    style_axes(ax, title="COHERENCE: late return probability\n"
                         "tracks $z_0^2$ — the electron's overlap with one line",
               xlabel="$V$", ylabel=r"$\langle P_{\rm return}\rangle_{\rm late}$")
    ax.legend(loc="upper right", fontsize=7.8)

    ax = axes[1]
    ax.plot(Vz, z0 * q4, color=INK, lw=1.4, ls="--", marker="*", ms=9,
            label=r"$z_0\,q_4$: carried by the edge eigenstate")
    Vb = sorted(V for (d, V) in bound if d == -0.30)
    ax.plot(Vb, [bound[(-0.30, V)]["zq4"] for V in Vb], color=INK_MUTED, lw=1.3,
            ls=":", marker="P", ms=7,
            label=r"$\sum z\,q_4$ over all lines below $\Delta_+$")
    for i, (L, inj, d, lbl) in enumerate((
            (16, "edge", -0.30, "edge, L = 16"),
            (20, "edge", -0.30, "edge, L = 20 (window scaled)"),
            (16, "centre", -0.30, "centre, L = 16 (4 central sites)"),
            (16, "edge", +0.30, r"edge $\delta=+0.30$, L = 16"))):
        V, y = late(L, d, inj, "Q4_late")
        if V.size:
            ax.plot(V, y, lw=1.5, ls=LINESTYLES[i], marker=MARKERS[i], ms=6,
                    color=CAT[i], label=lbl)
    ax.set_ylim(-0.02, 0.56)
    style_axes(ax, title="LOCALISATION: late excess charge in 4 sites\n"
                         "the topological edge keeps it; the trivial edge does not",
               xlabel="$V$", ylabel=r"$\langle Q_4\rangle_{\rm late}$")
    ax.legend(loc="upper center", fontsize=7.0, ncol=2, columnspacing=0.8)

    ax = axes[2]
    for i, V in enumerate(V_BIG):
        rr = sorted([r for r in far_rows if r["V"] == V], key=lambda r: r["L"])
        ax.plot([r["L"] for r in rr], [r["selected"] for r in rr], lw=1.5,
                ls=LINESTYLES[i], marker=MARKERS[i], ms=6, color=CAT[i],
                label=f"selected $|\\Omega\\rangle$, $V={V:g}$")
        ax.plot([r["L"] for r in rr], [r["left_loaded"] for r in rr], lw=0,
                marker=MARKERS[i], ms=7, mfc="none", mew=1.2, color=CAT[i])
    r = (1 - 0.30) / (1 + 0.30)
    lim = (1 - (0.7 / 1.3) ** 2) ** 2 / 2
    ax.axhline(lim, color=REFERENCE, lw=1.1, ls="--")
    ax.text(19.8, lim, rf"$V=0$, $L\to\infty$: $(1-r^2)^2/2={lim:.3f}$  ",
            fontsize=8, color="#55554f", ha="right", va="bottom")
    ax.set_ylim(-0.02, 0.36)
    ax.set_xticks([8, 12, 16, 18, 20])
    style_axes(ax, title=r"Far-edge excess at $t=0$, $\Delta n_{L-1}(0)$" "\n"
                         "hollow: a left-loaded state from the same doublet",
               xlabel="$L$", ylabel=r"$\Delta n_{L-1}(0)$")
    ax.legend(loc="center right", fontsize=7.6)
    fig.tight_layout()
    save(fig, "fig19_summary_vs_V")


# ---------------------------------------------------------------- io
def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    with path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, (bool, np.bool_)):
                    vals.append(str(bool(v)))
                elif isinstance(v, (int, float, np.integer, np.floating)):
                    vals.append(f"{v:.12g}")
                else:
                    vals.append(str(v))
            f.write(",".join(vals) + "\n")


def trace_csv(path, tr):
    rows = [{"t": t, "P_return": tr.P_return[i], "ReG": tr.amp[i].real,
             "ImG": tr.amp[i].imag, "Q2": tr.Q_local(2)[i], "Q4": tr.Q_local(4)[i],
             "Qright2": tr.Q_right(2)[i], "Qright4": tr.Q_right(4)[i],
             "sum_dn": tr.total_excess()[i], "norm": tr.norm[i],
             "energy": tr.energy[i]} for i, t in enumerate(tr.times)]
    write_csv(path, rows)


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args(argv)
    DATA.mkdir(parents=True, exist_ok=True)
    t_start = _time.time()

    log("Project A -- Step 4: watch an edge particle evolve\n" + "=" * 72)
    trajs: dict = {}
    runs: list[dict] = []
    filt_rows: list[dict] = []

    # ---- 1-3: main grid, L = 8 first
    log("1. Krylov evolution at L = 8, 12, 16 (dt = 0.05, t <= 8), tolerance "
        f"{dy.TOL:g} and {dy.TOL/2:g}")
    for L in L_MAIN:
        spec_mode = "full" if L <= 12 else "lanczos"
        for d in DELTAS:
            for V in V_BASE:
                for site in (0, centre_site(L)):
                    row, k, prep = run_one(L, d, V, site, spectral=spec_mode)
                    runs.append(row)
                    trajs[key(L, d, V, site)] = k
                    if site == 0:
                        frow, kf, _ = run_filtered(prep)
                        filt_rows.append(frow)
                        if kf is not None:
                            trajs[key(L, d, V, 0, "filtered")] = kf
        rr = [r for r in runs if r["L"] == L]
        log(f"   L = {L:2d}: {len(rr)} runs | exact "
            f"{np.nanmax([r['exact_amp'] for r in rr]):.1e} | spectral "
            f"{np.nanmax([r['spectral_amp'] for r in rr]):.1e} | halving "
            f"{np.nanmax([r['halving_amp'] for r in rr]):.1e} | norm "
            f"{max(r['norm_drift'] for r in rr):.1e} | energy "
            f"{max(r['energy_drift'] for r in rr):.1e} | sum dn "
            f"{max(r['sum_dn_minus_1'] for r in rr):.1e}  ({_time.time()-t_start:.0f}s)")

    # ---- tolerance scan
    log("\n   tolerance scan")
    scan_rows = []
    p8 = trajs[key(8, -0.30, 1.0, 0)].prep
    e8 = dy.exact_evolve(p8)
    for tol in TOLS_SCAN:
        k = dy.krylov_evolve(p8, tol=tol)
        c = dy.compare(k, e8)
        scan_rows.append({"L": 8, "V": 1.0, "tol": tol, "reference": "exact",
                          "err_amp": c["amp"], "err_dn": c["dn"],
                          "krylov_vectors": k.krylov_vectors_total})
    p16 = trajs[key(16, -0.30, 1.0, 0)].prep
    ref16 = dy.krylov_evolve(p16, tol=1e-13)
    for tol in TOLS_SCAN:
        k = dy.krylov_evolve(p16, tol=tol)
        c = dy.compare(k, ref16)
        scan_rows.append({"L": 16, "V": 1.0, "tol": tol, "reference": "tol=1e-13",
                          "err_amp": c["amp"], "err_dn": c["dn"],
                          "krylov_vectors": k.krylov_vectors_total})
    for r in scan_rows:
        log(f"     L={r['L']:2d} tol={r['tol']:.0e}: max|dG| = {r['err_amp']:.1e}, "
            f"Krylov vectors {r['krylov_vectors']}")

    # ---- 4: size ladder
    if not args.quick:
        log("\n2. size ladder L = 18, 20 (certifies the finite-size-free interval)")
        for d in DELTAS:
            for V in V_BASE:
                for site in (0, centre_site(18)):
                    long = (d == -0.30 and site == 0 and V in V_BIG)
                    row, k, prep = run_one(18, d, V, site,
                                           t_max=T_EXT if long else dy.T_MAX,
                                           halve=(d == -0.30 and site == 0 and V == 1.0),
                                           exact=False)
                    runs.append(row)
                    trajs[key(18, d, V, site)] = k
                    if d == -0.30 and site == 0 and V in V_BIG:
                        frow, _, _ = run_filtered(prep, evolve=False)
                        filt_rows.append(frow)
        log(f"   L = 18 done ({_time.time()-t_start:.0f}s)")
        for V in V_BIG:
            row, k, prep = run_one(20, -0.30, V, 0, t_max=T_EXT,
                                   halve=(V == 1.0), exact=False)
            runs.append(row)
            trajs[key(20, -0.30, V, 0)] = k
            frow, _, _ = run_filtered(prep, evolve=False)
            filt_rows.append(frow)
            log(f"   L = 20, V = {V:g}: dim {row['dim_Np1']}, norm drift "
                f"{row['norm_drift']:.1e}, energy drift {row['energy_drift']:.1e} "
                f"({_time.time()-t_start:.0f}s)")

    # ---- finite-size windows
    log("\n3. finite-size-free intervals (departure from the next larger L)")
    v_free = dy.free_front_velocity(-0.30)
    fs_rows = certify(trajs)
    fs_lookup, fs_meta = {}, {}
    for r in fs_rows:
        if r["eps"] == EPS_VISIBLE:
            k_ = (r["L"], r["delta"], r["V"], r["injection"])
            fs_lookup[k_] = r["t_fs"] if r["reached"] else np.inf
            fs_meta[k_] = r
    # L = 20 has no larger partner.  Its interval is L = 18's certified time
    # scaled by the ratio of round-trip distances, i.e. the same speed assumed.
    speeds = {}
    for V in V_BIG:
        r18 = fs_meta.get((18, -0.30, V, "edge"))
        if r18 is None:
            continue
        scale = dy.round_trip(20) / dy.round_trip(18)
        bound = not r18["reached"]
        t18 = r18["t_common_max"] if bound else r18["t_fs"]
        speeds[V] = {"v": dy.return_velocity(18, t18), "bound": bound,
                     "t18": t18, "t20": t18 * scale}
        fs_lookup[(20, -0.30, V, "edge")] = np.inf if bound else t18 * scale
        row20 = {"L": 20, "L_ref": "scaled from L=18", "delta": -0.30, "V": V,
                 "injection": "edge", "eps": EPS_VISIBLE, "t_fs": t18 * scale,
                 "reached": not bound, "lower_bound": bound, "extrapolated": True,
                 "t_common_max": T_EXT, "static_offset_max": np.nan,
                 "initial_state_ok": r18["initial_state_ok"]}
        fs_rows.append(row20)
        fs_meta[(20, -0.30, V, "edge")] = row20

    # ---- late-time averages inside the certified window
    filt_lookup = {(r["L"], r["delta"], r["V"]): r for r in filt_rows}
    for r in runs:
        k4 = (r["L"], r["delta"], r["V"], r["injection"])
        tr = trajs[key(r["L"], r["delta"], r["V"], r["site"])]
        meta_ = fs_meta.get(k4)
        r["initial_state_ok"] = bool(meta_["initial_state_ok"]) if meta_ else False
        r["static_offset_max"] = meta_.get("static_offset_max", np.nan) if meta_ else np.nan
        r["t_fs_extrapolated"] = bool(meta_.get("extrapolated", False)) if meta_ else False
        if meta_ is None or not r["initial_state_ok"]:
            # no larger chain to certify against, or the initial state itself is
            # not converged at the 1% level: no size-independent late window
            t_a, t_b = np.nan, np.nan
            r["t_fs_visible"] = np.nan if meta_ is None else meta_["t_fs"]
        else:
            t_fs = fs_lookup[k4]
            t_a, t_b = late_window(t_fs, r["t_max"], np.isfinite(t_fs))
            r["t_fs_visible"] = t_fs if np.isfinite(t_fs) else r["t_max"]
        r["late_from"], r["late_to"] = t_a, t_b
        if np.isfinite(t_a):
            m = (tr.times >= t_a) & (tr.times <= t_b + 1e-9)
            r["P_late"] = float(tr.P_return[m].mean())
            r["Q2_late"] = float(tr.Q_local(2)[m].mean())
            r["Q4_late"] = float(tr.Q_local(4)[m].mean())
            r["Qright4_late"] = float(tr.Q_right(4)[m].mean())
        else:
            r["P_late"] = r["Q2_late"] = r["Q4_late"] = r["Qright4_late"] = np.nan
        f = filt_lookup.get((r["L"], r["delta"], r["V"]))
        if f is not None and r["injection"] == "edge":
            r["z0W"] = f["z0W"]
            r["pred_P_late"] = f["z0W"] ** 2
            if not f["window_empty"]:
                r["pred_Q2_late"] = f["z0W"] * f["q2"]
                r["pred_Q4_late"] = f["z0W"] * f["q4"]

    # ---- beating inspection
    log("\n4. two-edge beating: exact evolution to t = 300 at L = 8")
    beat = []
    for V in V_BIG:
        prep = trajs[key(8, -0.30, V, 0)].prep
        e = dy.exact_evolve(prep, t_max=300.0)
        vals, vecs = np.linalg.eigh(prep.H.toarray())
        wts = np.abs(vecs.T @ prep.psi0) ** 2
        om = vals - prep.E0_N
        ok = wts > 1e-10
        diffs = np.abs(om[ok][:, None] - om[ok][None, :])
        wprod = (wts[ok][:, None] * wts[ok][None, :])
        sel = diffs > 1e-9
        split = prep.gs_gap
        near = int(np.argmin(np.abs(diffs[sel] - split)))
        q = e.Q_left(4) - e.Q_left(4).mean()
        freqs = np.linspace(0.005, 3.0, 3000)
        ph = np.exp(1j * np.outer(freqs, e.times))
        spec = 2 * np.abs(ph @ q) / e.times.size
        a_split = float(2 * abs(np.mean(q * np.exp(1j * split * e.times))))
        hann = np.hanning(e.times.size)
        spec_h = 2 * np.abs(ph @ (q * hann)) / hann.sum()
        a_split_h = float(2 * abs(np.sum(q * hann * np.exp(1j * split * e.times)))
                          / hann.sum())
        beat.append({"V": V, "splitting": split, "period": 2 * np.pi / split,
                     "nearest_line_difference": float(diffs[sel][near]),
                     "weight_of_that_pair": float(wprod[sel][near]),
                     "fourier_at_split": a_split,
                     "fourier_at_split_hann": a_split_h,
                     "fourier_max_hann": float(spec_h.max()),
                     "fourier_max": float(spec.max()),
                     "fourier_max_at": float(freqs[int(np.argmax(spec))]),
                     "Q_left4_mean": float(e.Q_left(4).mean()),
                     "Q_right4_mean": float(e.Q_right(4).mean()),
                     "freqs": freqs, "spectrum": spec, "spectrum_hann": spec_h})
        log(f"   V = {V:g}: splitting {split:.4f} (period {2*np.pi/split:.1f}); "
            f"Fourier (Hann) at splitting {a_split_h:.1e} vs strongest "
            f"{spec_h.max():.1e}; plain {a_split:.1e} vs {spec.max():.1e} at "
            f"{freqs[int(np.argmax(spec))]:.3f}; long-time halves "
            f"{e.Q_left(4).mean():.4f} / {e.Q_right(4).mean():.4f}")
    periods = []
    for L in ([8, 12, 16] if args.quick else [8, 12, 16, 18, 20]):
        for V in V_BIG:
            k_ = key(L, -0.30, V, 0)
            if k_ in trajs:
                g = trajs[k_].prep.gs_gap
                periods.append({"L": L, "V": V, "splitting": g,
                                "beat_period": 2 * np.pi / g})

    # ---- far-edge response
    log("\n5. far-edge response at t = 0: selected vs symmetry-broken states")
    far_rows = []
    for L in ([8, 12, 16] if args.quick else [8, 12, 16, 18, 20]):
        for V in V_BIG:
            r = dy.far_edge_response(L, -0.30, V)
            far_rows.append({"L": L, "V": V,
                             "selected": r["selected"]["dn_far"],
                             "left_loaded": r["left_loaded"]["dn_far"],
                             "right_loaded": r["right_loaded"]["dn_far"],
                             "n0_selected": r["selected"]["n0"],
                             "n0_left_loaded": r["left_loaded"]["n0"],
                             "n0_right_loaded": r["right_loaded"]["n0"],
                             "doublet_splitting": r["doublet_splitting"]})
        log(f"   L = {L:2d}: " + "; ".join(
            f"V={x['V']:g}: {x['selected']:.4f} (broken {x['left_loaded']:.1e})"
            for x in far_rows if x["L"] == L))

    # ---- bound lines below the bulk threshold (what cannot leave the edge)
    log("\n5b. lines below the ring threshold: bound edge charge")
    line_rows, bound = [], {}
    for d in DELTAS:
        for V in V_BASE:
            prep = trajs[key(16, d, V, 0)].prep
            dp = spx.ring_addition_threshold(16, d, V)
            lines = dy.subthreshold_lines(prep, dp)
            for x in lines:
                line_rows.append({"L": 16, "delta": d, "V": V, "delta_plus": dp, **x})
            bound[(d, V)] = {"zq4": sum(x["z"] * x["q4"] for x in lines),
                             "zq2": sum(x["z"] * x["q2"] for x in lines),
                             "n": len(lines)}
            if d < 0 or V in (0.0, 3.0):
                log(f"   delta={d:+.2f} V={V:g}: {len(lines)} line(s) below Delta+="
                    f"{dp:.3f}; " + "; ".join(f"omega={x['omega']:.3f} z={x['z']:.3f} "
                                              f"q4={x['q4']:.3f}" for x in lines)
                    + f";  sum z*q4 = {bound[(d, V)]['zq4']:.4f}")
    for r in runs:
        b = bound.get((r["delta"], r["V"]))   # L = 16 lines, converged in L
        if b is not None and r["injection"] == "edge":
            r["bound_Q4"] = b["zq4"]
            r["bound_Q2"] = b["zq2"]
            r["n_bound_lines"] = b["n"]

    # ---- write data
    write_csv(DATA / "step4_runs.csv", runs + filt_rows)
    if line_rows:
        write_csv(DATA / "step4_bound_lines.csv", line_rows)
    write_csv(DATA / "step4_finite_size.csv", fs_rows)
    write_csv(DATA / "step4_tolerance_scan.csv", scan_rows)
    write_csv(DATA / "step4_beating.csv",
              [{k: v for k, v in b.items()
                if k not in ("freqs", "spectrum", "spectrum_hann")}
               for b in beat] + periods)
    write_csv(DATA / "step4_far_edge.csv", far_rows)
    np.savez_compressed(
        DATA / "step4_trajectories.npz",
        **{("_".join(str(x) for x in k)).replace(".", "p").replace("-", "m")
           .replace("+", "p"): np.column_stack(
               [tr.times, tr.amp.real, tr.amp.imag, tr.norm, tr.energy, tr.dn])
           for k, tr in trajs.items()})
    for V in V_BIG:
        for k_ in (key(16, -0.30, V, 0), key(16, -0.30, V, centre_site(16)),
                   key(16, +0.30, V, 0), key(16, -0.30, V, 0, "filtered"),
                   key(20, -0.30, V, 0)):
            if k_ in trajs:
                name = "_".join(str(x) for x in k_).replace(".", "p")
                trace_csv(DATA / f"step4_traces_{name}.csv", trajs[k_])
    log("\n   wrote data/step4_*.csv and step4_trajectories.npz")

    # ---- figures
    log("\n6. figures")
    fig_validation(trajs, scan_rows, runs)
    fig_spacetime(trajs, fs_lookup, v_free, speeds)
    fig_traces(trajs, filt_rows, fs_lookup, None)
    fig_finite_size(trajs, fs_rows, beat, speeds)
    fig_summary(runs, filt_rows, far_rows, bound)

    # ---- headline
    log("\n7. finite-size-free intervals, edge injection, delta = -0.30 "
        "(eps = 0.01 + 2 x static offset)")
    for L in [8, 12, 16, 18, 20]:
        cells = []
        for V in V_BASE:
            x = fs_meta.get((L, -0.30, V, "edge"))
            if x is None:
                cells.append("     —     ")
            elif not x["initial_state_ok"]:
                cells.append(f"s={x['static_offset_max']:.3f}!")
            elif x.get("extrapolated"):
                cells.append(f"t<~{x['t_fs']:.1f}" + ("+" if x.get("lower_bound") else "*"))
            elif x["reached"]:
                cells.append(f"t<{x['t_fs']:.2f}")
            else:
                cells.append(f"t<={x['t_common_max']:.0f}+")
        log(f"   L = {L:2d}: " + "  ".join(f"V={V:g}: {c:<11}" for V, c in zip(V_BASE, cells)))
    log("   s=..! : initial state not converged (static offset > 0.01 vs the next L);"
        "  + : no return seen to the end of the run;  * : scaled from L = 18")
    for V, f in speeds.items():
        log(f"   V = {V:g}: speed implied by the L = 18 return "
            f"{'<=' if f['bound'] else '='} {f['v']:.2f} sites/time "
            f"(free maximal group velocity {v_free:.2f})")

    log("\n8. late-time comparison with Step 3, L = 16, delta = -0.30, edge")
    log(f"   {'V':>4} {'z0(W)':>7} {'<P>late':>8} {'z0^2':>7} {'<Q4>late':>9} "
        f"{'z0*q4':>7} {'bound Q4':>8} {'centre Q4':>9} {'trivial Q4':>10}")
    for V in V_BASE:
        e = next(r for r in runs if r["L"] == 16 and r["delta"] == -0.30
                 and r["V"] == V and r["injection"] == "edge")
        c = next(r for r in runs if r["L"] == 16 and r["delta"] == -0.30
                 and r["V"] == V and r["injection"] == "centre")
        t = next(r for r in runs if r["L"] == 16 and r["delta"] == 0.30
                 and r["V"] == V and r["injection"] == "edge")
        log(f"   {V:>4.1f} {e['z0W']:>7.4f} {e['P_late']:>8.4f} {e['pred_P_late']:>7.4f} "
            f"{e['Q4_late']:>9.4f} {e['pred_Q4_late']:>7.4f} {e['bound_Q4']:>8.4f} "
            f"{c['Q4_late']:>9.4f} "
            f"{t['Q4_late']:>10.4f}")
    if not args.quick:
        log("   L = 20, extended to t = 16 (window scaled from L = 18):")
        for V in V_BIG:
            e = next(r for r in runs if r["L"] == 20 and r["V"] == V)
            log(f"   {V:>4.1f} {e['z0W']:>7.4f} {e['P_late']:>8.4f} "
                f"{e['pred_P_late']:>7.4f} {e['Q4_late']:>9.4f} "
                f"{e['pred_Q4_late']:>7.4f} {e['bound_Q4']:>8.4f}   "
                f"window [{e['late_from']:.2f}, {e['late_to']:.2f}]")

    meta = {
        "project": "A -- can a topological edge survive strong interactions?",
        "step": 4,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "scipy": scipy.__version__, "platform": platform.platform(),
        "command": " ".join([Path(sys.argv[0]).name] + sys.argv[1:]),
        "runtime_s": round(_time.time() - t_start, 1),
        "conventions": {
            "inherited_from": "Steps 1-3 (ssh_free, manybody, spectral)",
            "preparation": "|Psi(0)> = c^dag_j|Omega_N>/sqrt(1-<n_j>), Omega by the "
                           "Step 2 parity+ rule",
            "evolution": "exp(-i (H_{N+1} - E_0^N) t); the shift is kept so G(t) is "
                         "compared with Step 3 as a complex amplitude",
            "Q_edge": "sum_{j<l_e} Delta n_j, l_e = 2, 4",
            "centre_injection_site": "2*(L//4): even site nearest the centre",
            "krylov": f"short-iterative Lanczos, full reorthogonalisation, error "
                      f"estimate beta0*beta_m*|e_m^T exp(-iT dt) e_1| <= tol*dt, "
                      f"tol = {dy.TOL} per unit time, halved to {dy.TOL/2}",
            "finite_size_window": f"first t where P_return or Q_l (l=2,4) of L "
                                  f"departs from the next larger L by more than "
                                  f"{EPS_VISIBLE} ('visible'); {EPS_TIGHT} also "
                                  f"recorded",
            "late_window": "second half of the certified window, from t >= 2",
            "filtered": "projection onto Step 3's widest isolated window",
            "random_seed": None,
        },
        "grid": {"L_main": L_MAIN, "L_ladder": [] if args.quick else L_LADDER,
                 "delta": DELTAS, "V": V_BASE, "V_ladder": V_BIG,
                 "dt": dy.DT_OUT, "t_max": dy.T_MAX, "t_extended": T_EXT},
        "implied_return_speed": {str(V): {"v": f["v"], "upper_bound": f["bound"]}
                                 for V, f in speeds.items()},
        "free_max_group_velocity": v_free,
        "caveats": [
            "A revival in a short chain is not a long lifetime in the infinite system.",
            "Late-time averages use only the certified finite-size-free interval.",
            "The t = 0 far-edge excess is a ground-state correlation of the "
            "parity-selected |Omega>, not transport; it vanishes for a "
            "symmetry-broken choice from the same doublet.",
            "Spectral filtering onto one line manufactures a stationary state: "
            "P_return = 1 identically. It is reported as a filtered excitation.",
        ],
    }
    (DATA / "step4_metadata.json").write_text(json.dumps(meta, indent=2,
                                                         default=float))
    log(f"\n  wrote data/step4_metadata.json  (total {_time.time()-t_start:.0f}s)")
    log("\nStep 4 complete.  Run `python test_step4.py` for the gate checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
