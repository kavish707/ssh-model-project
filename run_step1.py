#!/usr/bin/env python3
"""
Project A -- Step 1 driver.  Runs the whole free-edge programme and writes
data tables, a run-metadata record and the figures.

    python run_step1.py                # full grid
    python run_step1.py --quick        # L = 40, 80 only

Outputs
-------
data/step1_summary.csv        one row per (L, delta): edge energies, splitting,
                              fitted and analytic xi, fit quality, edge weight
data/xi_vs_nfit.csv           fitted xi for n_fit = 6..10 (fit-window test)
data/splitting_vs_L.csv       end-state splitting over a sweep in L
data/ring_gaps.csv            band gap and addition cost vs the 4|d|, 2|d| laws
data/spectrum_L*_d*.csv       full single-particle spectra
data/edge_profile_L*_d*.csv   left-localised edge amplitudes and probabilities
data/run_metadata.json        versions, grid, conventions, benchmark check
figures/fig1..fig5 (.png/.pdf)
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

import ssh_free as sf
from plot_style import CAT, LINESTYLES, MARKERS, REFERENCE, use_style, style_axes

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"

L_GRID = [40, 80, 120]
DELTA_GRID = [-0.15, -0.30, -0.50, 0.15, 0.30, 0.50]
DELTA_NEG = [-0.15, -0.30, -0.50]
L_SWEEP = list(range(16, 161, 4))
N_FIT_DEFAULT = 8
XI_BENCHMARK = 3.230813126  # guide value at delta = -0.30


# ---------------------------------------------------------------- helpers
def tag(L: int, delta: float) -> str:
    return f"L{L}_d{delta:+.2f}".replace(".", "p")


def analyse(L: int, delta: float, n_fit: int = N_FIT_DEFAULT) -> dict:
    """All Step-1 quantities for one (L, delta)."""
    res = sf.localised_edge_states(L, delta)
    psi = res["psi_left"]
    weak_end = delta < 0

    row = {
        "L": L,
        "delta": delta,
        "N": L // 2,
        "boundary": "open",
        "J0": float(sf.bond_strengths(L, delta)[0]),
        "Jw": sf.weak_strong(delta)[0],
        "Js": sf.weak_strong(delta)[1],
        "weak_end_bond": bool(weak_end),
        "E_minus": float(res["energies"][0]),
        "E_plus": float(res["energies"][1]),
        "splitting": res["splitting"],
        "bulk_gap_ring": 4.0 * abs(delta),
        "x_left": res["x_left"],
        "x_right": res["x_right"],
        "ipr_left": float(np.sum(psi**4)),
        "even_weight": float(np.sum(psi[0::2] ** 2)),
        "weight_first4": float(np.sum(psi[:4] ** 2)),
    }

    if weak_end:
        xi_an = sf.xi_analytic(delta)
        fit = sf.fit_localisation_length(psi, n_fit=n_fit, sublattice=0)
        overlap = abs(float(psi @ sf.analytic_edge_profile(L, delta)))
        row.update(
            {
                "xi_analytic": xi_an,
                "xi_approx_1_over_delta": 1.0 / abs(delta),
                "xi_fit_amplitude": fit["xi_amp"],
                "xi_fit_probability": fit["xi_prob"],
                "xi_rel_err_amplitude": abs(fit["xi_amp"] - xi_an) / xi_an,
                "xi_rel_err_probability": abs(fit["xi_prob"] - xi_an) / xi_an,
                "fit_r2_amplitude": fit["r2_amp"],
                "fit_max_resid_amplitude": fit["max_resid_amp"],
                "fit_site_min": int(fit["fit_sites"][0]),
                "fit_site_max": int(fit["fit_sites"][-1]),
                "n_fit": n_fit,
                "overlap_with_analytic": overlap,
            }
        )
    else:
        row.update(
            {k: np.nan for k in
             ("xi_analytic", "xi_approx_1_over_delta", "xi_fit_amplitude",
              "xi_fit_probability", "xi_rel_err_amplitude",
              "xi_rel_err_probability", "fit_r2_amplitude",
              "fit_max_resid_amplitude", "overlap_with_analytic")}
        )
        row.update({"fit_site_min": -1, "fit_site_max": -1, "n_fit": 0})
    row["_psi_left"] = psi
    row["_psi_right"] = res["psi_right"]
    row["_prob_sum"] = res["prob_sum"]
    row["_spectrum"] = res["all_energies"]
    return row


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join("" if v is None else
                             (f"{v:.12g}" if isinstance(v, (int, float, np.floating)) else str(v))
                             for v in r) + "\n")


# ---------------------------------------------------------------- figures
def fig_spectra(results: dict, L: int) -> None:
    use_style()
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.8), sharex=True, sharey=True)
    for col, mag in enumerate([0.15, 0.30, 0.50]):
        for row_i, sgn in enumerate([-1.0, 1.0]):
            d = sgn * mag
            ax = axes[row_i, col]
            r = results[(L, d)]
            E = r["_spectrum"]
            n = np.arange(E.size)
            ax.plot(n, E, lw=0, marker="o", ms=2.2, color=CAT[0])
            i0, i1 = sf.edge_pair_indices(E)
            ax.plot([i0, i1], [E[i0], E[i1]], lw=0, marker="o", ms=6.5,
                    mfc="none", mew=1.6, color=CAT[1])
            ax.axhline(0.0, color=REFERENCE, lw=0.9, ls=":")
            for s in (+1, -1):
                ax.axhline(s * 2 * abs(d), color=REFERENCE, lw=0.9, ls="--")
            kind = "weak end bonds" if d < 0 else "strong end bonds"
            ax.set_title(rf"$\delta={d:+.2f}$  ({kind})")
            style_axes(ax)
            txt = (rf"$|E|_{{\min}}={abs(E[i0]):.2e}$" if d < 0
                   else rf"$|E|_{{\min}}={abs(E[i0]):.3f}$")
            ax.text(0.5, 0.06, txt, transform=ax.transAxes, ha="center",
                    fontsize=8.5, color="#55554f")
    for ax in axes[1]:
        ax.set_xlabel("eigenvalue index $n$")
    for ax in axes[:, 0]:
        ax.set_ylabel("$E_n$  (units of $J$)")
    fig.suptitle(f"Single-particle spectrum of the open SSH chain, L = {L}\n"
                 "orange rings: the two states nearest zero   ·   "
                 r"dashed: $\pm 2|\delta|$, the half bulk gap",
                 x=0.012, y=0.995, ha="left", va="top", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, "fig1_spectra")


def fig_edge_profiles(results: dict, L: int) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    for i, d in enumerate(DELTA_NEG):
        r = results[(L, d)]
        p = r["_psi_left"] ** 2
        j = np.arange(L)
        m = p > 1e-18
        ax.semilogy(j[m], p[m], lw=0, marker=MARKERS[i], ms=3.4, color=CAT[i],
                    label=rf"$\delta={d:+.2f}$   $\xi={r['xi_analytic']:.3f}$")
        xi = r["xi_analytic"]
        jj = np.linspace(0, min(L / 2, 10 * xi), 200)
        ax.semilogy(jj, (1 - (sf.weak_strong(d)[0] / sf.weak_strong(d)[1]) ** 2)
                    * np.exp(-2 * jj / xi), ls=LINESTYLES[i], lw=1.2,
                    color=REFERENCE, zorder=1)
    ax.set_ylim(1e-18, 5)
    ax.set_xlim(0, L / 2)
    style_axes(ax, title=rf"Left edge state, $|\psi_j|^2$   (L = {L})",
               xlabel="site $j$", ylabel=r"$|\psi_j|^2$")
    ax.legend(loc="upper right", title="grey curves: $ (1-r^2)e^{-2j/\\xi}$",
              title_fontsize=8)

    ax = axes[1]
    d, floor = -0.30, 1e-19
    r = results[(L, d)]
    j = np.arange(L)

    def _m(y):
        y = np.asarray(y, float)
        mask = y > floor
        return j[mask], y[mask]

    ax.semilogy(*_m(r["_prob_sum"]), lw=2.6, color=CAT[0], alpha=0.45,
                label=r"raw doublet  $|\psi_1|^2+|\psi_2|^2$  (rotation invariant)")
    ax.semilogy(*_m(r["_psi_left"] ** 2), lw=0, marker="s", ms=3.0,
                color=CAT[1], label=r"left-localised  $|\psi_L|^2$")
    ax.semilogy(*_m(r["_psi_right"] ** 2), lw=0, marker="^", ms=3.0,
                color=CAT[2], label=r"right-localised  $|\psi_R|^2$")
    ax.set_ylim(floor, 20)
    style_axes(ax, title=rf"Hybridised pair vs. position-diagonalised states  "
                         rf"($\delta={d:+.2f}$, L = {L})",
               xlabel="site $j$", ylabel=r"probability")
    ax.legend(loc="upper center", ncol=1)
    fig.tight_layout()
    save(fig, "fig2_edge_profiles")


def fit_window_study(L: int, delta: float, n_max: int = 18):
    """xi_fit vs n_fit, with and without the guide's 'exclude the opposite edge' rule.

    Returns (n_guide, xi_guide, n_naive, xi_naive).  The naive curve lets the fit
    window run past the chain midpoint, where the exponentially small tail of the
    *other* edge state and double-precision noise take over and the fitted length
    drifts away from the exact value.
    """
    psi = sf.localised_edge_states(L, delta)["psi_left"]
    n_g, xi_g, n_n, xi_n = [], [], [], []
    for nf in range(6, n_max + 1):
        try:
            xi_g.append(sf.fit_localisation_length(psi, n_fit=nf, max_site_frac=0.5)["xi_amp"])
            n_g.append(nf)
        except ValueError:
            pass
        try:
            xi_n.append(sf.fit_localisation_length(psi, n_fit=nf, max_site_frac=1.0,
                                                   floor=0.0)["xi_amp"])
            n_n.append(nf)
        except (ValueError, FloatingPointError):
            pass
    return np.array(n_g), np.array(xi_g), np.array(n_n), np.array(xi_n)


def fig_xi(summary: list[dict], nfit_rows: list[list]) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    dd = np.linspace(-0.85, -0.02, 400)
    ax.plot(np.abs(dd), [sf.xi_analytic(x) for x in dd], color=REFERENCE, lw=1.4,
            label=r"$\xi = 2/\ln(J_s/J_w)$")
    ax.plot(np.abs(dd), 1.0 / np.abs(dd), color=REFERENCE, lw=1.1, ls=":",
            label=r"$1/|\delta|$")
    for i, L in enumerate(sorted({r["L"] for r in summary})):
        xs = [abs(r["delta"]) for r in summary if r["L"] == L and r["delta"] < 0]
        ys = [r["xi_fit_amplitude"] for r in summary if r["L"] == L and r["delta"] < 0]
        ax.plot(xs, ys, lw=0, marker=MARKERS[i], ms=7, mfc="none", mew=1.6,
                color=CAT[i], label=f"fit, L = {L}")
    ax.set_xlim(0.05, 0.85)
    ax.set_ylim(0, 9)
    style_axes(ax, title=r"Fitted vs. analytic localisation length",
               xlabel=r"$|\delta|$", ylabel=r"$\xi$  (sites)")
    ax.legend(loc="upper right")

    ax = axes[1]
    Lw, dw = 40, -0.15                       # worst case: longest xi, shortest chain
    xi_ex = sf.xi_analytic(dw)
    n_g, xi_g, n_n, xi_n = fit_window_study(Lw, dw)
    ax.axhline(xi_ex, color=REFERENCE, lw=1.3,
               label=rf"exact $\xi={xi_ex:.4f}$")
    ax.plot(n_g, xi_g, lw=1.6, ls="-", marker="o", color=CAT[0],
            label=r"fit window kept to $j<L/2$ (as instructed)")
    ax.plot(n_n, xi_n, lw=1.6, ls="--", marker="s", color=CAT[1],
            label="window allowed past the midpoint")
    ax.axvline(len(n_g) + 5.5, color=REFERENCE, lw=0.8, ls=":")
    ax.text(0.42, 0.93, "window reaches\nthe far edge", transform=ax.transAxes,
            fontsize=8, color="#55554f", va="top")
    style_axes(ax, title=rf"Fit-window test  ($L={Lw}$, $\delta={dw:+.2f}$)",
               xlabel=r"number of even sites in the fit, $n_{\rm fit}$",
               ylabel=r"$\xi_{\rm fit}$  (sites)")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "fig3_xi")


def fig_splitting(sweep: dict) -> None:
    use_style()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for i, d in enumerate(DELTA_NEG):
        Ls = np.array(sorted(sweep[d]))
        sp = np.array([sweep[d][L] for L in Ls])
        m = sp > 1e-15
        ax.semilogy(Ls[m], sp[m], lw=1.4, ls=LINESTYLES[i], marker=MARKERS[i],
                    ms=4, color=CAT[i], label=rf"$\delta={d:+.2f}$")
        if (~m).any():
            ax.semilogy(Ls[~m], np.full((~m).sum(), 1e-15), lw=0, marker="x",
                        ms=4, color=CAT[i])
        xi = sf.xi_analytic(d)
        ref = sp[m][0] * np.exp(-(Ls[m] - Ls[m][0]) / xi)
        ax.semilogy(Ls[m], ref, color="#1c1c1a", lw=1.0, ls=(0, (1, 2)),
                    zorder=6, alpha=0.75)
    ax.axhline(1e-15, color=REFERENCE, lw=0.8)
    ax.text(150, 1.6e-15, "double precision floor", ha="right", fontsize=8,
            color="#8a8a82")
    ax.set_ylim(1e-16, 1)
    style_axes(ax, title="End-state splitting falls exponentially with L\n"
                         r"dotted: $e^{-L/\xi}$   ·   $\times$: below double precision",
               xlabel="chain length $L$", ylabel=r"$E_+ - E_-$")
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "fig4_splitting_vs_L")


def fig_ring(ring_rows: list[dict]) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    k = np.linspace(-np.pi, np.pi, 801)
    for i, d in enumerate(DELTA_NEG):
        Em, Ep = sf.bloch_bands(d, k)
        ax.plot(k, Ep, color=CAT[i], lw=1.6, ls=LINESTYLES[i],
                label=rf"$\delta={d:+.2f}$,  gap $=4|\delta|={4*abs(d):.1f}$")
        ax.plot(k, Em, color=CAT[i], lw=1.6, ls=LINESTYLES[i])
    # allowed momenta of a finite ring; the exact ED levels agree with these to
    # ~1e-15 (column bloch_max_abs_error in data/ring_gaps.csv)
    Lr, dr = 40, -0.30
    kk = 2 * np.pi * np.arange(Lr // 2) / (Lr // 2)
    kk = np.where(kk > np.pi, kk - 2 * np.pi, kk)
    Em_k, Ep_k = sf.bloch_bands(dr, kk)
    ax.plot(np.concatenate([kk, kk]), np.concatenate([Em_k, Ep_k]),
            lw=0, marker="o", ms=4, mfc="none", mew=1.2, color="#1c1c1a",
            label=rf"allowed $k$, ring L = {Lr}, $\delta={dr:+.2f}$")
    ax.set_xticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
    ax.set_ylim(-3.45, 2.35)
    style_axes(ax, title=r"Ring dispersion $E_\pm(k)=\pm\sqrt{J_w^2+J_s^2+2J_wJ_s\cos k}$",
               xlabel="$k$", ylabel="$E$")
    ax.legend(loc="lower center", fontsize=7.6, ncol=1)

    ax = axes[1]
    mags = np.array(sorted({abs(r["delta"]) for r in ring_rows}))
    gap = [np.mean([r["band_gap_numeric"] for r in ring_rows
                    if np.isclose(abs(r["delta"]), m)]) for m in mags]
    add = [np.mean([r["addition_cost_numeric"] for r in ring_rows
                    if np.isclose(abs(r["delta"]), m)]) for m in mags]
    xx = np.linspace(0, 0.6, 100)
    ax.plot(xx, 4 * xx, color=REFERENCE, lw=1.3, label=r"$4|\delta|$")
    ax.plot(xx, 2 * xx, color=REFERENCE, lw=1.3, ls="--", label=r"$2|\delta|$")
    ax.plot(mags, gap, lw=0, marker="o", ms=8, mfc="none", mew=1.8, color=CAT[0],
            label="band gap (ring ED)")
    ax.plot(mags, add, lw=0, marker="s", ms=8, mfc="none", mew=1.8, color=CAT[1],
            label="bulk addition cost (ring ED)")
    style_axes(ax, title="Bulk gaps against the analytic laws",
               xlabel=r"$|\delta|$", ylabel="energy (units of $J$)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "fig5_ring")


def save(fig, name: str) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="L = 40, 80 only")
    ap.add_argument("--n-fit", type=int, default=N_FIT_DEFAULT)
    args = ap.parse_args(argv)

    L_list = [40, 80] if args.quick else L_GRID
    DATA.mkdir(parents=True, exist_ok=True)

    print("Project A -- Step 1: free SSH edge problem")
    print(f"  grid: L = {L_list},  delta = {DELTA_GRID}")

    results, summary = {}, []
    for L in L_list:
        for d in DELTA_GRID:
            r = analyse(L, d, n_fit=args.n_fit)
            results[(L, d)] = r
            summary.append(r)
            write_csv(DATA / f"spectrum_{tag(L, d)}.csv", ["n", "E"],
                      [[n, e] for n, e in enumerate(r["_spectrum"])])
            if d < 0:
                psi = r["_psi_left"]
                write_csv(DATA / f"edge_profile_{tag(L, d)}.csv",
                          ["j", "psi_left", "prob_left", "prob_doublet_sum",
                           "psi_analytic"],
                          [[j, psi[j], psi[j] ** 2, r["_prob_sum"][j],
                            sf.analytic_edge_profile(L, d)[j]] for j in range(L)])

    cols = [k for k in summary[0] if not k.startswith("_")]
    write_csv(DATA / "step1_summary.csv", cols,
              [[r[k] for k in cols] for r in summary])
    print(f"  wrote data/step1_summary.csv  ({len(summary)} rows)")

    # ---- fit-window test
    nfit_rows = []
    for L in L_list:
        for d in DELTA_NEG:
            xi_an = sf.xi_analytic(d)
            for nf in range(6, 11):
                fit = sf.fit_localisation_length(results[(L, d)]["_psi_left"], n_fit=nf)
                nfit_rows.append([L, d, nf, (fit["xi_amp"] - xi_an) / xi_an,
                                  fit["xi_amp"], fit["xi_prob"], xi_an, fit["r2_amp"]])
    write_csv(DATA / "xi_vs_nfit.csv",
              ["L", "delta", "n_fit", "rel_err", "xi_fit_amp", "xi_fit_prob",
               "xi_analytic", "r2"], nfit_rows)
    print("  wrote data/xi_vs_nfit.csv")

    # ---- splitting sweep
    sweep = {d: {} for d in DELTA_NEG}
    rows = []
    for d in DELTA_NEG:
        for L in L_SWEEP:
            s = sf.localised_edge_states(L, d)["splitting"]
            sweep[d][L] = s
            rows.append([L, d, s, np.exp(-L / sf.xi_analytic(d))])
    write_csv(DATA / "splitting_vs_L.csv",
              ["L", "delta", "splitting", "exp(-L/xi)"], rows)
    print("  wrote data/splitting_vs_L.csv")

    # ---- ring diagnostics
    ring_rows = [sf.ring_gaps(L, d) for L in L_list for d in DELTA_GRID]
    rc = list(ring_rows[0].keys())
    write_csv(DATA / "ring_gaps.csv", rc, [[r[k] for k in rc] for r in ring_rows])
    print("  wrote data/ring_gaps.csv")

    # ---- figures
    Lfig = L_list[-1]
    fig_spectra(results, L_list[1] if len(L_list) > 1 else L_list[0])
    fig_edge_profiles(results, Lfig)
    fig_xi(summary, nfit_rows)
    fig_splitting(sweep)
    fig_ring(ring_rows)

    # ---- metadata + headline check
    bench = next(r for r in summary if r["L"] == Lfig and r["delta"] == -0.30)
    meta = {
        "project": "A -- can a topological edge survive strong interactions?",
        "step": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "command": " ".join([Path(sys.argv[0]).name] + sys.argv[1:]),
        "conventions": {
            "hbar": 1, "lattice_spacing": 1, "J": 1,
            "site_labels_start_at": 0,
            "J_j": "1 + (-1)^j * delta",
            "h[j,j+1]": "-J_j",
            "filling": "N = L/2",
            "V": 0.0,
            "random_seed": None,
        },
        "grid": {"L": L_list, "delta": DELTA_GRID, "n_fit": args.n_fit},
        "benchmark_delta_-0.30": {
            "xi_guide": XI_BENCHMARK,
            "xi_analytic": bench["xi_analytic"],
            "xi_fit_amplitude": bench["xi_fit_amplitude"],
            "xi_fit_probability": bench["xi_fit_probability"],
            "abs_err_vs_guide": abs(bench["xi_fit_amplitude"] - XI_BENCHMARK),
            "L": Lfig,
        },
    }
    (DATA / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    print("  wrote data/run_metadata.json")

    print("\n  benchmark, delta = -0.30:")
    print(f"    xi (guide)      {XI_BENCHMARK:.9f}")
    print(f"    xi (analytic)   {bench['xi_analytic']:.9f}")
    print(f"    xi (fit, amp)   {bench['xi_fit_amplitude']:.9f}   "
          f"rel. err {bench['xi_rel_err_amplitude']:.2e}")
    print(f"    xi (fit, prob)  {bench['xi_fit_probability']:.9f}   "
          f"rel. err {bench['xi_rel_err_probability']:.2e}")
    print("\nStep 1 complete.  Run `python test_step1.py` for the gate checks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
