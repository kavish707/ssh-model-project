#!/usr/bin/env python3
"""
Project A -- Step 2 driver: a many-body calculation you can trust.

    python run_step2.py              # validation suite + the first interacting grid
    python run_step2.py --quick      # L = 8, 12 only, no L = 18
    python run_step2.py --no-l18     # skip the L = 18 representative points

Order of work follows the guide: validate at L = 6, 8, 10 by full diagonalisation,
then finish L = 8 at every grid point, then enlarge to L = 12 and 16, then add
L = 18 at representative points, refining any V interval where something moves.

Outputs
-------
data/step2_validation.csv     small-size checks: Hermiticity, residuals, <H> vs E,
                              V=0 agreement with the Step 1 free results
data/step2_grid.csv           the interacting grid, one row per (L, delta, V)
data/step2_degeneracy.csv     neutral gap vs L, i.e. where state selection starts
                              to matter
data/step2_metadata.json      versions, grid, conventions, tolerances
figures/fig6..fig9 (.png/.pdf)
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
import scipy.sparse as sp

import manybody as mb
import ssh_free as sf
from plot_style import CAT, LINESTYLES, MARKERS, REFERENCE, use_style, style_axes

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGS = HERE / "figures"

L_SMALL = [6, 8, 10]                       # full diagonalisation
L_GRID = [8, 12, 16]                       # the first interacting grid
L_BIG = [18]                               # representative points only
DELTAS = [-0.30, +0.30]
V_BASE = [0.0, 0.5, 1.0, 2.0, 3.0]
V_BIG = [0.0, 1.0, 3.0]
REFINE_STEP = 0.25

RES_TOL_SMALL = 1e-10                      # the guide's tolerance for small tests
RES_TOL_PROD = 1e-9                        # ... and for production states
CONVENTION = "parity+"                     # see manybody.select_ground_state


# ---------------------------------------------------------------- one grid point
def grid_point(L: int, delta: float, V: float, k: int = 6) -> dict:
    """Everything Step 2 asks to be recorded, for one (L, delta, V)."""
    N = L // 2
    sec = mb.Sector.make(L, N)
    H = mb.build_hamiltonian(sec, delta, V)
    g = mb.select_ground_state(H, sec, k=k, convention=CONVENTION)

    # independently evaluated expectation value, not the eigenvalue
    e_expect = float(g.psi @ (H @ g.psi))

    # N+1 and N-1 sectors, built independently
    sec_p, sec_m = mb.Sector.make(L, N + 1), mb.Sector.make(L, N - 1)
    Hp = mb.build_hamiltonian(sec_p, delta, V)
    Hm = mb.build_hamiltonian(sec_m, delta, V)
    ep, vp = mb.lowest_states(Hp, k=2)
    em, vm = mb.lowest_states(Hm, k=2)

    # the same chain on a ring: the bulk reference Step 3 needs
    Hr = mb.build_hamiltonian(sec, delta, V, periodic=True)
    Hrp = mb.build_hamiltonian(sec_p, delta, V, periodic=True)
    er, _ = mb.lowest_states(Hr, k=2)
    erp, _ = mb.lowest_states(Hrp, k=1)

    C = mb.correlation_matrix(g.psi, L, N)
    n = mb.density_profile(g.psi, sec)
    occ = sec.occupations() - 0.5

    return {
        "L": L, "N": N, "delta": delta, "V": V, "boundary": "open",
        "dim_N": sec.dim, "dim_Np1": sec_p.dim, "dim_Nm1": sec_m.dim,
        "solver": "dense" if sec.dim <= 2000 else "lanczos",
        "E0": g.energy, "E1": float(g.energies[1]), "E2": float(g.energies[2]),
        "E3": float(g.energies[3]),
        "residual_E0": g.residual,
        "residual_E1": mb.residual(H, g.energies[1], _vec(H, g, 1)),
        "E_expect_minus_eigenvalue": abs(e_expect - g.energy),
        "neutral_gap": float(g.energies[1] - g.energies[0]),
        "bulk_gap": float(g.energies[2] - g.energies[0]),
        "cluster_size": g.cluster_size, "cluster_spread": g.cluster_spread,
        "parity": g.parity, "parity_bare": g.parity_bare,
        "state_convention": g.convention, "rotated_in_cluster": g.rotated,
        "E0_Np1": float(ep[0]), "E0_Nm1": float(em[0]),
        "residual_Np1": mb.residual(Hp, ep[0], vp[:, 0]),
        "residual_Nm1": mb.residual(Hm, em[0], vm[:, 0]),
        "addition_gap": float(ep[0] - g.energy),
        "removal_gap": float(em[0] - g.energy),
        "E0_ring": float(er[0]), "neutral_gap_ring": float(er[1] - er[0]),
        "addition_gap_ring": float(erp[0] - er[0]),
        "ring_L_div4_warning": bool(L % 4 == 0),
        "n0": float(n[0]), "n1": float(n[1]),
        "max_dev_half_filling": float(np.max(np.abs(n - 0.5))),
        "C01": float(C[0, 1]), "C03": float(C[0, 3]),
        "edge_interaction": float((occ[:, 0] * occ[:, 1] * g.psi**2).sum()),
        "trace_C": float(np.trace(C)),
    }


def _vec(H, g, i):
    """The i-th eigenvector, recomputed consistently with select_ground_state."""
    _, vecs = mb.lowest_states(H, k=max(4, i + 1))
    return vecs[:, i]


# ---------------------------------------------------------------- validation
def validation_suite() -> list[dict]:
    """L = 6, 8, 10 by full diagonalisation, with every check the guide lists."""
    rows = []
    for L in L_SMALL:
        N = L // 2
        for delta in DELTAS:
            for V in V_BASE:
                for periodic in (False, True):
                    sec = mb.Sector.make(L, N)
                    H = mb.build_hamiltonian(sec, delta, V, periodic=periodic)
                    Hd = np.asarray(H.todense())
                    vals, vecs = np.linalg.eigh(Hd)
                    psi = vecs[:, 0]
                    row = {
                        "L": L, "N": N, "delta": delta, "V": V,
                        "boundary": "periodic" if periodic else "open",
                        "dim": sec.dim, "method": "dense",
                        "hermiticity": float(np.max(np.abs(Hd - Hd.T))),
                        "E0": float(vals[0]),
                        "residual": mb.residual(H, vals[0], psi),
                        "E_expect_minus_eigenvalue":
                            abs(float(psi @ (Hd @ psi)) - vals[0]),
                        "number_conserved": True,   # fixed-N basis, verified in jw_check
                    }
                    if V == 0.0:
                        free = sf.ground_state_energy(L, delta, periodic=periodic)
                        C = mb.correlation_matrix(psi, L, N)
                        Cf = mb.free_correlation_matrix(L, delta, N, periodic=periodic)
                        row["free_energy_check"] = abs(vals[0] - free)
                        row["free_corr_check"] = float(np.max(np.abs(C - Cf)))
                    else:
                        row["free_energy_check"] = np.nan
                        row["free_corr_check"] = np.nan
                    rows.append(row)
    return rows


# ---------------------------------------------------------------- refinement
WATCH = ("addition_gap", "bulk_gap", "neutral_gap", "C01")


def needs_refinement(a: dict, b: dict) -> tuple[bool, str]:
    """Did anything move between two neighbouring V values?

    V = 2 is the uniform-chain transition; it is NOT assumed to be special for a
    dimerised chain, so no interval is refined just because it contains V = 2.
    """
    reasons = []
    for key in WATCH:
        x, y = a[key], b[key]
        scale = max(abs(x), abs(y), 1e-12)
        if abs(y - x) / scale > 0.25:
            reasons.append(f"{key} {x:.4g}->{y:.4g}")
    if a["parity_bare"] * b["parity_bare"] < 0:
        reasons.append("ground-state parity flips")
    if a["cluster_size"] != b["cluster_size"]:
        reasons.append("degeneracy changes")
    return bool(reasons), "; ".join(reasons)


# ---------------------------------------------------------------- figures
def fig_validation(val_rows, grid_rows) -> None:
    use_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))

    ax = axes[0]
    for i, (lbl, key) in enumerate((("N sector", "residual_E0"),
                                    ("N+1 sector", "residual_Np1"),
                                    ("N-1 sector", "residual_Nm1"))):
        dims = [r["dim_N"] for r in grid_rows]
        ax.loglog(dims, [max(r[key], 1e-18) for r in grid_rows], lw=0,
                  marker=MARKERS[i], ms=5, mfc="none", mew=1.4, color=CAT[i],
                  label=lbl)
    ax.axhline(RES_TOL_PROD, color=REFERENCE, lw=1.2, ls="--")
    ax.text(0.98, 0.52, r"guide tolerance $10^{-9}$", transform=ax.transAxes,
            ha="right", fontsize=8, color="#55554f")
    ax.set_ylim(1e-17, 1e-7)
    style_axes(ax, title="Eigenpair residuals", xlabel="sector dimension",
               ylabel=r"$\|H\psi-E\psi\|$")
    ax.legend(loc="lower left")

    ax = axes[1]
    v0 = [r for r in val_rows if r["V"] == 0.0]
    for i, (lbl, key) in enumerate(((r"$|E_{\rm ED}-\sum_{a<N}\epsilon_a|$",
                                     "free_energy_check"),
                                    (r"$\max|C_{\rm ED}-C_{\rm free}|$",
                                     "free_corr_check"))):
        ax.semilogy([r["L"] for r in v0], [max(r[key], 1e-18) for r in v0],
                    lw=0, marker=MARKERS[i], ms=6, mfc="none", mew=1.5,
                    color=CAT[i], label=lbl)
    ax.set_ylim(1e-18, 1e-10)
    ax.set_xticks(L_SMALL)
    style_axes(ax, title="$V=0$: many body vs. the Step 1 free solution",
               xlabel="$L$", ylabel="deviation")
    ax.legend(loc="upper left")

    ax = axes[2]
    herm = [max(r["hermiticity"], 1e-18) for r in val_rows]
    expct = [max(r["E_expect_minus_eigenvalue"], 1e-18) for r in val_rows]
    ax.semilogy(range(len(herm)), herm, lw=0, marker="o", ms=3.4, color=CAT[0],
                label=r"$\max|H-H^{T}|$")
    ax.semilogy(range(len(expct)), expct, lw=0, marker="s", ms=3.4, color=CAT[1],
                label=r"$|\langle\psi|H|\psi\rangle-E|$")
    ax.set_ylim(1e-18, 1e-10)
    style_axes(ax, title=f"Hermiticity and $\\langle H\\rangle$ vs. $E$\n"
                         f"({len(val_rows)} small-size points, open and periodic)",
               xlabel="validation point", ylabel="deviation")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "fig6_step2_validation")


def fig_spectrum_vs_V(grid_rows) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True)
    for col, d in enumerate(DELTAS):
        ax = axes[col]
        for i, L in enumerate(sorted({r["L"] for r in grid_rows})):
            rows = sorted([r for r in grid_rows
                           if r["L"] == L and r["delta"] == d],
                          key=lambda r: r["V"])
            Vs = [r["V"] for r in rows]
            for lvl, style in ((1, "-"), (2, "--")):
                ax.plot(Vs, [r[f"E{lvl}"] - r["E0"] for r in rows], ls=style,
                        lw=1.6, marker=MARKERS[i], ms=4, color=CAT[i],
                        label=f"L = {L}, $E_{lvl}-E_0$" if col == 0 else None)
        ax.set_yscale("log")
        kind = "weak end bonds" if d < 0 else "strong end bonds"
        style_axes(ax, title=rf"$\delta={d:+.2f}$  ({kind})",
                   xlabel="$V$", ylabel=r"$E_n-E_0$" if col == 0 else None)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Low-lying many-body spectrum of the open chain\n"
                 "solid: first excited (the edge doublet for weak end bonds)   ·   "
                 "dashed: second excited",
                 x=0.012, y=0.995, ha="left", va="top", fontsize=11)
    fig.tight_layout(rect=[0, 0.09, 1, 0.90])
    save(fig, "fig7_spectrum_vs_V")


def fig_gaps_vs_V(grid_rows) -> None:
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True)
    for col, d in enumerate(DELTAS):
        ax = axes[col]
        for i, L in enumerate(sorted({r["L"] for r in grid_rows})):
            rows = sorted([r for r in grid_rows
                           if r["L"] == L and r["delta"] == d],
                          key=lambda r: r["V"])
            Vs = [r["V"] for r in rows]
            ax.plot(Vs, [r["addition_gap"] for r in rows], ls="-", lw=1.7,
                    marker=MARKERS[i], ms=4.5, color=CAT[i],
                    label=f"L = {L}, open" if col == 0 else None)
            ax.plot(Vs, [r["addition_gap_ring"] for r in rows], ls="--", lw=1.4,
                    marker=MARKERS[i], ms=4.5, mfc="none", color=CAT[i],
                    label=f"L = {L}, ring" if col == 0 else None)
        ax.axhline(2 * abs(d), color=REFERENCE, lw=1.1, ls=":")
        ax.text(0.98, 2 * abs(d), rf"free bulk $2|\delta|={2*abs(d):.1f}$  ",
                transform=ax.get_yaxis_transform(), fontsize=8,
                color="#55554f", va="bottom", ha="right")
        kind = "weak end bonds" if d < 0 else "strong end bonds"
        style_axes(ax, title=rf"$\delta={d:+.2f}$  ({kind})", xlabel="$V$",
                   ylabel=r"$E_0^{N+1}-E_0^{N}$" if col == 0 else None)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Addition gap: open chain (solid) against the ring at the same "
                 "parameters (dashed)\n"
                 "the ring value is the bulk threshold the Step 3 edge window uses",
                 x=0.012, y=0.995, ha="left", va="top", fontsize=11)
    fig.tight_layout(rect=[0, 0.09, 1, 0.90])
    save(fig, "fig8_gaps_vs_V")


def fig_state_selection(degen_rows) -> None:
    """The guide's warning, made quantitative."""
    use_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))

    # (a) a resolved eigenstate vs. an arbitrary vector from the same pair
    ax = axes[0]
    L, N, d, V = 16, 8, -0.30, 1.0
    sec = mb.Sector.make(L, N)
    H = mb.build_hamiltonian(sec, d, V)
    vals, vecs = mb.lowest_states(H, k=2)
    mix = (vecs[:, 0] + vecs[:, 1]) / np.sqrt(2)
    j = np.arange(1, L)
    series = (
        (r"$E_0$, parity $+$ (selected)", vecs[:, 0], CAT[0], "o", 4.6, "-"),
        (r"$E_1$, parity $-$ (lies on top of $E_0$)", vecs[:, 1], CAT[1], "s", 7.5, ""),
        (r"arbitrary mixture $(\psi_0+\psi_1)/\sqrt{2}$", mix, CAT[2], "^", 4.6, "-."),
    )
    spurious = {}
    for lbl, psi, col, mk, ms, ls in series:
        C = np.abs(mb.correlation_matrix(psi, L, N)[0, 1:])
        ax.plot(j, np.maximum(C, 1e-17), lw=1.5 if ls else 0, ls=ls or "-",
                marker=mk, ms=ms, color=col, label=lbl,
                mfc="none" if mk == "s" else None, mew=1.3 if mk == "s" else None)
        spurious[lbl] = C[1::2].max()          # even sites j = 2, 4, ...
    ax.set_yscale("log")
    ax.set_ylim(1e-17, 3e4)
    ax.text(0.5, 0.985,
            "on even sites an eigenstate gives "
            rf"$\leq {spurious[series[0][0]]:.0e}$," "\n"
            rf"the mixture gives ${spurious[series[2][0]]:.2f}$ — "
            "a spurious edge correlation",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.2,
            color="#55554f")
    style_axes(ax, title=rf"Edge correlator $|\langle c^\dagger_0 c_j\rangle|$  "
                         rf"($L={L}$, $\delta={d:+.2f}$, $V={V:g}$)",
               xlabel="site $j$", ylabel=r"$|C_{0j}|$")
    ax.legend(loc="lower left", fontsize=7.8)

    # (b) when does the splitting stop being resolvable?
    ax = axes[1]
    xi = sf.xi_analytic(-0.30)
    for i, V in enumerate(sorted({r["V"] for r in degen_rows})):
        rows = sorted([r for r in degen_rows if r["V"] == V], key=lambda r: r["L"])
        Ls = np.array([r["L"] for r in rows], float)
        gaps = np.array([r["neutral_gap"] for r in rows])
        ax.semilogy(Ls, gaps, lw=1.6, ls=LINESTYLES[i], marker=MARKERS[i], ms=4.5,
                    color=CAT[i], label=f"computed, $V={V:g}$")
        Lx = np.linspace(Ls[-1], 80, 50)
        ax.semilogy(Lx, gaps[-1] * np.exp(-(Lx - Ls[-1]) / xi), lw=1.0,
                    ls=(0, (1, 2)), color=CAT[i], alpha=0.8)
    ax.axhline(mb.DEGENERACY_TOL, color=REFERENCE, lw=1.2, ls="--")
    L_cross = 18 + xi * np.log(
        min(r["neutral_gap"] for r in degen_rows if r["L"] == 18) / mb.DEGENERACY_TOL)
    ax.axvline(L_cross, color=REFERENCE, lw=0.9, ls=":")
    ax.annotate(rf"$L\approx{L_cross:.0f}$", xy=(L_cross, 1e-7), fontsize=8.5,
                color="#55554f", ha="left", va="center",
                xytext=(L_cross + 2, 1e-7))
    ax.text(0.97, 0.30,
            r"below $10^{-9}$ the pair is unresolved and" "\n"
            r"the solver returns an arbitrary vector",
            transform=ax.transAxes, ha="right", fontsize=8, color="#55554f")
    ax.set_ylim(1e-13, 3)
    ax.set_xlim(6, 82)
    style_axes(ax, title=r"Edge-doublet splitting vs. $L$  ($\delta=-0.30$)"
                         "\n" r"dotted: extrapolation at the Step 1 rate $e^{-L/\xi}$",
               xlabel="$L$", ylabel=r"$E_1-E_0$")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    save(fig, "fig9_state_selection")


def save(fig, name: str) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(
                f"{r[c]:.12g}" if isinstance(r[c], (int, float, np.floating))
                and not isinstance(r[c], bool) else str(r[c]) for c in cols) + "\n")


# ---------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-l18", action="store_true")
    args = ap.parse_args(argv)

    L_list = [8, 12] if args.quick else L_GRID
    DATA.mkdir(parents=True, exist_ok=True)

    print("Project A -- Step 2: many-body ED\n" + "=" * 72)
    print("1. validation at L = 6, 8, 10 (full diagonalisation)")
    val = validation_suite()
    write_csv(DATA / "step2_validation.csv", val)
    worst_res = max(r["residual"] for r in val)
    worst_herm = max(r["hermiticity"] for r in val)
    worst_exp = max(r["E_expect_minus_eigenvalue"] for r in val)
    worst_free = max(r["free_energy_check"] for r in val if not np.isnan(r["free_energy_check"]))
    worst_corr = max(r["free_corr_check"] for r in val if not np.isnan(r["free_corr_check"]))
    print(f"   {len(val)} points | residual <= {worst_res:.1e} "
          f"(need {RES_TOL_SMALL:.0e}) | Hermiticity <= {worst_herm:.1e}")
    print(f"   <H> vs E <= {worst_exp:.1e} | V=0 energy vs Step 1 <= {worst_free:.1e}"
          f" | V=0 correlator <= {worst_corr:.1e}")

    print("\n2. first interacting grid (L = 8 everywhere, then enlarge)")
    grid: list[dict] = []
    for L in L_list:
        for d in DELTAS:
            base = []
            for V in V_BASE:
                r = grid_point(L, d, V)
                base.append(r)
                grid.append(r)
            # refine any interval where something moved
            for a, b in zip(base, base[1:]):
                need, why = needs_refinement(a, b)
                if need:
                    extra = np.arange(a["V"] + REFINE_STEP, b["V"], REFINE_STEP)
                    print(f"   refining L={L} delta={d:+.2f} V in "
                          f"({a['V']:g}, {b['V']:g}): {why}")
                    for V in extra:
                        grid.append(grid_point(L, d, float(V)))
            print(f"   L = {L:2d}, delta = {d:+.2f}: {len(base)} base points, "
                  f"dim {base[0]['dim_N']}, max residual "
                  f"{max(x['residual_E0'] for x in base):.1e}")

    if not args.quick and not args.no_l18:
        print("\n3. L = 18 at representative points")
        for d in DELTAS:
            for V in V_BIG:
                r = grid_point(18, d, V)
                grid.append(r)
                print(f"   L = 18, delta = {d:+.2f}, V = {V:g}: "
                      f"dim {r['dim_N']}, E0 = {r['E0']:.9f}, "
                      f"residual {r['residual_E0']:.1e}")

    grid.sort(key=lambda r: (r["L"], r["delta"], r["V"]))
    write_csv(DATA / "step2_grid.csv", grid)
    print(f"\n   wrote data/step2_grid.csv ({len(grid)} rows); "
          f"max production residual {max(r['residual_E0'] for r in grid):.1e}")

    print("\n4. degeneracy scan: where does the state-selection convention matter?")
    degen = []
    for V in (0.0, 1.0, 3.0):
        for L in ([8, 10, 12, 14] if args.quick else [8, 10, 12, 14, 16, 18]):
            sec = mb.Sector.make(L, L // 2)
            H = mb.build_hamiltonian(sec, -0.30, V)
            g = mb.select_ground_state(H, sec, k=4, convention=CONVENTION)
            degen.append({"L": L, "delta": -0.30, "V": V,
                          "neutral_gap": g.gap_to_first_excited,
                          "cluster_size": g.cluster_size,
                          "parity_bare": g.parity_bare,
                          "rotated_in_cluster": g.rotated,
                          "residual": g.residual})
    write_csv(DATA / "step2_degeneracy.csv", degen)
    xi = sf.xi_analytic(-0.30)
    L_crit = xi * np.log(degen[0]["neutral_gap"] / mb.DEGENERACY_TOL) + degen[0]["L"]
    print(f"   gap ~ e^(-L/xi) with xi = {xi:.3f} from Step 1; it crosses the "
          f"1e-9 tolerance near L ~ {L_crit:.0f}")

    # a case that IS degenerate now: the uniform ring at L divisible by four
    sec = mb.Sector.make(12, 6)
    Hu = mb.build_hamiltonian(sec, 0.0, 0.0, periodic=True)
    gu = mb.select_ground_state(Hu, sec, k=4, convention=CONVENTION, periodic=True)
    print(f"   uniform ring L = 12, V = 0: cluster size {gu.cluster_size}, "
          f"spread {gu.cluster_spread:.1e} -- genuinely degenerate, as the guide "
          f"notes for L divisible by four")

    print("\n5. figures")
    fig_validation(val, grid)
    fig_spectrum_vs_V([r for r in grid if r["V"] in V_BASE])
    fig_gaps_vs_V([r for r in grid if r["V"] in V_BASE])
    fig_state_selection(degen)

    meta = {
        "project": "A -- can a topological edge survive strong interactions?",
        "step": 2,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "numpy": np.__version__,
        "scipy": scipy.__version__, "platform": platform.platform(),
        "command": " ".join([Path(sys.argv[0]).name] + sys.argv[1:]),
        "conventions": {
            "inherited_from": "Step 1 / ssh_free.py",
            "J_j": "1 + (-1)^j * delta", "H": "Eq. (A.1)",
            "basis": "s = sum_j n_j 2^j, fixed N",
            "fermion_sign": "c_j, c^dag_j give (-1)^{sum_{l<j} n_l}; for c^dag_j c_k "
                            "annihilate at k first, second sign on the intermediate state",
            "periodic_closing_bond": "(L-1, 0) with J_{L-1} = 1 + (-1)^{L-1} delta, "
                                     "parity factor from the same rule",
            "state_selection": CONVENTION,
            "state_selection_rule":
                "collect the cluster within 1e-9 of E0; if every member already has "
                "|<R>| = 1 the symmetry is resolved and the lowest state is taken; "
                "otherwise R is diagonalised inside the cluster and the + eigenvector "
                "chosen. The quoted parity is the bare bit-reversal eigenvalue, which "
                "does not flip with L.",
            "random_seed": None,
        },
        "grid": {"L": L_list, "L_big": [] if (args.quick or args.no_l18) else L_BIG,
                 "delta": DELTAS, "V_base": V_BASE, "refine_step": REFINE_STEP,
                 "watched_observables": list(WATCH)},
        "tolerances": {"residual_small": RES_TOL_SMALL,
                       "residual_production": RES_TOL_PROD,
                       "degeneracy": mb.DEGENERACY_TOL},
        "worst": {"residual_validation": worst_res,
                  "residual_production": max(r["residual_E0"] for r in grid),
                  "hermiticity": worst_herm,
                  "expectation_vs_eigenvalue": worst_exp,
                  "V0_energy_vs_step1": worst_free,
                  "V0_correlator_vs_step1": worst_corr},
        "note_on_V2": "V = 2 is the uniform-chain transition. No interval was "
                      "refined merely because it contains V = 2; refinement is "
                      "triggered only by a measured change in the watched observables.",
    }
    (DATA / "step2_metadata.json").write_text(json.dumps(meta, indent=2))
    print("  wrote data/step2_metadata.json")
    print("\nStep 2 complete.  Run `python test_step2.py` and `python jw_check.py`.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
