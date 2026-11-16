#!/usr/bin/env python3
"""
Project A -- Step 1 gate checks.

    python test_step1.py          # prints a table, exits 1 if anything fails
    pytest test_step1.py -q       # same checks under pytest

The three gates the guide sets before Step 2 are checks 6, 8 and 9:
  6  the fitted length agrees with xi = 2/ln(J_s/J_w)
  8  the end-state splitting decreases with L
  9  reversing delta removes the low-energy end pair for this termination
The others are the convention and consistency checks that make those meaningful.
"""

from __future__ import annotations

import sys

import numpy as np

import ssh_free as sf

L_GRID = [40, 80, 120]
DELTA_NEG = [-0.15, -0.30, -0.50]
DELTA_POS = [0.15, 0.30, 0.50]
XI_BENCHMARK = 3.230813126


# ------------------------------------------------------------------ checks
def check_bond_convention():
    """J_j = 1+(-1)^j delta, so the first bond is J_0 = 1+delta = J_w when delta<0."""
    msgs = []
    for d in DELTA_NEG + DELTA_POS:
        J = sf.bond_strengths(40, d)
        Jw, Js = sf.weak_strong(d)
        assert np.isclose(J[0], 1.0 + d), f"J_0 != 1+delta at delta={d}"
        assert np.isclose(J[0], Jw), f"J_0 is not J_w at delta={d}"
        assert np.allclose(J[0::2], 1.0 + d) and np.allclose(J[1::2], 1.0 - d)
        if d < 0:
            assert J[0] < J[1], f"delta={d} should give a WEAK first bond"
        else:
            assert J[0] > J[1], f"delta={d} should give a STRONG first bond"
        msgs.append(f"d={d:+.2f}: J0={J[0]:.2f} J1={J[1]:.2f}")
    return "  ".join(msgs)


def check_hamiltonian_structure():
    """Real, symmetric, nearest-neighbour only, and h[j,j+1] = -J_j."""
    L, d = 40, -0.30
    h = sf.hopping_matrix(L, d)
    J = sf.bond_strengths(L, d)
    assert np.allclose(h, h.T), "not symmetric"
    assert np.allclose(np.diag(h), 0.0), "nonzero on-site energies"
    for j in range(L - 1):
        assert np.isclose(h[j, j + 1], -J[j]), f"h[{j},{j+1}] != -J_{j}"
    mask = np.ones((L, L), bool)
    for j in range(L - 1):
        mask[j, j + 1] = mask[j + 1, j] = False
    np.fill_diagonal(mask, False)
    assert np.allclose(h[mask], 0.0), "hopping beyond nearest neighbours"
    return f"L={L}: symmetric, tri-diagonal, max|h|={np.abs(h).max():.2f}"


def check_chiral_symmetry():
    """Bipartite lattice => the spectrum is symmetric under E -> -E."""
    worst = 0.0
    for L in L_GRID:
        for d in DELTA_NEG + DELTA_POS:
            E, _ = sf.single_particle_spectrum(L, d)
            worst = max(worst, float(np.max(np.abs(E + E[::-1]))))
    assert worst < 1e-12, f"spectrum not chiral-symmetric, worst {worst:.2e}"
    return f"max |E_n + E_(L-1-n)| = {worst:.2e}"


def check_edge_state_sublattice():
    """The left zero mode lives on the even sublattice (A sites) only."""
    out = []
    for L in L_GRID:
        for d in DELTA_NEG:
            psi = sf.localised_edge_states(L, d)["psi_left"]
            w_even = float(np.sum(psi[0::2] ** 2))
            assert w_even > 1 - 1e-9, f"L={L} d={d}: even-site weight {w_even}"
            out.append(w_even)
    return f"even-sublattice weight >= {min(out):.12f}"


def check_analytic_profile_overlap():
    """The numerical left state matches psi_2l = sqrt(1-r^2)(-r)^l.

    The finite chain differs from the semi-infinite solution by the admixture of
    the *other* edge, so the tolerance is not a fixed number: it is the finite-size
    parameter exp(-2L/xi) (times a margin), floored at double precision.
    """
    worst, worst_margin = 1.0, 0.0
    for L in L_GRID:
        for d in DELTA_NEG:
            psi = sf.localised_edge_states(L, d)["psi_left"]
            ov = abs(float(psi @ sf.analytic_edge_profile(L, d)))
            tol = max(10.0 * np.exp(-2.0 * L / sf.xi_analytic(d)), 1e-13)
            assert abs(1.0 - ov) < tol, (
                f"L={L} delta={d}: 1-overlap={abs(1-ov):.2e} exceeds "
                f"finite-size bound {tol:.2e}")
            worst = min(worst, ov)
            worst_margin = max(worst_margin, abs(1 - ov) / tol)
    return (f"min overlap = {worst:.12f};  worst deviation is "
            f"{worst_margin:.2f}x its exp(-2L/xi) bound")


def check_benchmark_xi():
    """GATE: at delta = -0.30 the fitted amplitude length is 3.230813126."""
    L, d = 120, -0.30
    psi = sf.localised_edge_states(L, d)["psi_left"]
    fit = sf.fit_localisation_length(psi, n_fit=8)
    xi_an = sf.xi_analytic(d)
    assert abs(xi_an - XI_BENCHMARK) < 5e-10, f"formula gives {xi_an}"
    assert abs(fit["xi_amp"] - XI_BENCHMARK) < 1e-8, f"fit gives {fit['xi_amp']}"
    return (f"xi_formula={xi_an:.9f}  xi_fit={fit['xi_amp']:.9f}  "
            f"guide={XI_BENCHMARK:.9f}  |diff|={abs(fit['xi_amp']-XI_BENCHMARK):.2e}")


def check_fit_matches_formula():
    """GATE: the fitted length agrees with the formula across the whole grid."""
    worst, where, worst_margin = 0.0, None, 0.0
    for L in L_GRID:
        for d in DELTA_NEG:
            psi = sf.localised_edge_states(L, d)["psi_left"]
            xi_an = sf.xi_analytic(d)
            tol = max(10.0 * np.exp(-L / xi_an), 1e-13)
            for nf in (6, 8, 10):
                fit = sf.fit_localisation_length(psi, n_fit=nf)
                rel = abs(fit["xi_amp"] - xi_an) / xi_an
                assert rel < tol, (
                    f"L={L} delta={d} n_fit={nf}: relative error {rel:.2e} "
                    f"exceeds the finite-size bound {tol:.2e}")
                worst_margin = max(worst_margin, rel / tol)
                if rel > worst:
                    worst, where = rel, (L, d, nf)
    return (f"worst relative error {worst:.2e} at (L,delta,n_fit)={where}; "
            f"every point within {worst_margin:.2f}x its exp(-L/xi) bound")


def check_probability_factor_two():
    """Probability decays at exactly twice the amplitude rate (the factor of 2 in A.3)."""
    worst = 0.0
    for L in L_GRID:
        for d in DELTA_NEG:
            psi = sf.localised_edge_states(L, d)["psi_left"]
            fit = sf.fit_localisation_length(psi, n_fit=8)
            worst = max(worst, abs(fit["slope_prob"] - 2 * fit["slope_amp"]))
    assert worst < 1e-10, f"slope_prob != 2*slope_amp, worst {worst:.2e}"
    return f"max |slope_prob - 2 slope_amp| = {worst:.2e}"


def check_splitting_decreases():
    """GATE: the end-state splitting falls monotonically with L (until it hits 1e-13)."""
    notes = []
    for d in DELTA_NEG:
        Ls = list(range(16, 161, 4))
        sp = [sf.localised_edge_states(L, d)["splitting"] for L in Ls]
        resolved = [(L, s) for L, s in zip(Ls, sp) if s > 1e-13]
        vals = [s for _, s in resolved]
        assert len(vals) >= 5, f"too few resolved points at delta={d}"
        assert all(b < a for a, b in zip(vals, vals[1:])), \
            f"splitting not monotonically decreasing at delta={d}"
        # the decay rate must be the analytic 1/xi to a few per cent
        Lv = np.array([L for L, _ in resolved], float)
        rate = -np.polyfit(Lv, np.log(vals), 1)[0]
        xi_fit, xi_an = 1.0 / rate, sf.xi_analytic(d)
        assert abs(xi_fit - xi_an) / xi_an < 0.05, \
            f"decay length {xi_fit:.3f} vs xi={xi_an:.3f} at delta={d}"
        notes.append(f"d={d:+.2f}: {vals[0]:.2e}->{vals[-1]:.2e}, "
                     f"1/rate={xi_fit:.3f} vs xi={xi_an:.3f}")
    return " | ".join(notes)


def check_reversed_delta_removes_edge_pair():
    """GATE: delta > 0 (strong end bonds) has no low-energy end pair."""
    notes = []
    for L in L_GRID:
        for mag in (0.15, 0.30, 0.50):
            E_neg, _ = sf.single_particle_spectrum(L, -mag)
            E_pos, _ = sf.single_particle_spectrum(L, +mag)
            i0n, _ = sf.edge_pair_indices(E_neg)
            i0p, _ = sf.edge_pair_indices(E_pos)
            lo_neg, lo_pos = abs(E_neg[i0n]), abs(E_pos[i0p])
            half_gap = 2.0 * mag
            xi = sf.xi_analytic(-mag)
            assert lo_neg < 0.05 * half_gap, \
                f"no deep in-gap state at L={L}, delta={-mag}: |E|={lo_neg:.2e}"
            assert lo_neg < max(10.0 * np.exp(-L / xi), 1e-13), \
                f"in-gap energy {lo_neg:.2e} above the hybridisation scale at L={L}"
            assert lo_pos > 0.95 * half_gap, \
                f"spurious in-gap state at L={L}, delta={+mag}: |E|={lo_pos:.4f}"
            if L == L_GRID[-1]:
                notes.append(f"|d|={mag:.2f}: |E|min = {lo_neg:.1e} (d<0) vs "
                             f"{lo_pos:.4f} (d>0), half gap {half_gap:.2f}")
    return " | ".join(notes)


def check_ring_gap_laws():
    """Ring band gap = 4|delta| and bulk addition cost = 2|delta|."""
    worst_gap = worst_add = 0.0
    for L in L_GRID:
        for d in DELTA_NEG + DELTA_POS:
            r = sf.ring_gaps(L, d)
            worst_gap = max(worst_gap, abs(r["band_gap_numeric"] - 4 * abs(d)))
            worst_add = max(worst_add, abs(r["addition_cost_numeric"] - 2 * abs(d)))
    assert worst_gap < 1e-10 and worst_add < 1e-10, \
        f"gap err {worst_gap:.2e}, addition err {worst_add:.2e}"
    return f"max |gap-4|d|| = {worst_gap:.2e},  max |E_add-2|d|| = {worst_add:.2e}"


def check_bloch_matches_ed():
    """E_pm(k) at the allowed momenta reproduces the exact ring spectrum."""
    worst = 0.0
    for L in L_GRID:
        for d in DELTA_NEG + DELTA_POS:
            worst = max(worst, sf.ring_gaps(L, d)["bloch_max_abs_error"])
    assert worst < 1e-12, f"worst Bloch-vs-ED error {worst:.2e}"
    return f"max |E_ED - E_Bloch| = {worst:.2e}"


def check_free_ground_state():
    """At V=0 and half filling the ground state is the sum of all negative levels."""
    worst = 0.0
    for L in L_GRID:
        for d in DELTA_NEG + DELTA_POS:
            E, _ = sf.single_particle_spectrum(L, d)
            worst = max(worst, abs(sf.ground_state_energy(L, d) - np.sum(E[E < 0])))
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max |E_GS - sum(E<0)| = {worst:.2e}  (hook for Step 2, check 4)"


def check_position_rotation_invariance():
    """|psi_L|^2+|psi_R|^2 equals |psi_1|^2+|psi_2|^2: the sum is rotation invariant."""
    worst = 0.0
    for L in L_GRID:
        for d in DELTA_NEG:
            r = sf.localised_edge_states(L, d)
            lhs = r["psi_left"] ** 2 + r["psi_right"] ** 2
            worst = max(worst, float(np.max(np.abs(lhs - r["prob_sum"]))))
    assert worst < 1e-12, f"worst {worst:.2e}"
    return f"max deviation {worst:.2e}"


CHECKS = [
    ("1  bond convention  J_0 = J_w for delta<0", check_bond_convention),
    ("2  Hamiltonian structure", check_hamiltonian_structure),
    ("3  chiral symmetry of the spectrum", check_chiral_symmetry),
    ("4  edge state on the even sublattice", check_edge_state_sublattice),
    ("5  overlap with the analytic zero mode", check_analytic_profile_overlap),
    ("6  GATE  benchmark xi at delta = -0.30", check_benchmark_xi),
    ("7  GATE  fitted xi = formula, whole grid", check_fit_matches_formula),
    ("8  probability decays twice as fast", check_probability_factor_two),
    ("9  GATE  splitting decreases with L", check_splitting_decreases),
    ("10 GATE  reversed delta kills the end pair", check_reversed_delta_removes_edge_pair),
    ("11 ring gap = 4|delta|, addition = 2|delta|", check_ring_gap_laws),
    ("12 Bloch bands reproduce the ring ED", check_bloch_matches_ed),
    ("13 free half-filled ground-state energy", check_free_ground_state),
    ("14 doublet sum is rotation invariant", check_position_rotation_invariance),
]

# pytest entry points
for _name, _fn in CHECKS:
    globals()["test_" + _fn.__name__[6:]] = _fn


def main() -> int:
    print("Project A -- Step 1 gate checks\n" + "=" * 78)
    failed = 0
    for label, fn in CHECKS:
        try:
            detail = fn()
            print(f"[PASS] {label}\n       {detail}")
        except AssertionError as exc:
            failed += 1
            print(f"[FAIL] {label}\n       {exc}")
    print("=" * 78)
    if failed:
        print(f"{failed} of {len(CHECKS)} checks FAILED -- do not go on to Step 2.")
    else:
        print(f"all {len(CHECKS)} checks passed -- Step 1 gates are clear for Step 2.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
