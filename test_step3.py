#!/usr/bin/env python3
"""
Project A -- Step 3 gate checks.

    python test_step3.py          # prints a table, exits 1 if anything fails
    pytest test_step3.py -q       # same checks under pytest

The guide's "before moving on" for Step 3 is that xi_W(V) and z_0(V) are reported
side by side, and that missing weight, an unstable window or a poor exponential
fit stay visible.  Checks marked GATE cover the quantities those statements rest
on: the sum rule, the V = 0 closed forms, rotation invariance of a cluster
weight, and the Lanczos measure against full ED.
"""

from __future__ import annotations

import sys

import numpy as np

import manybody as mb
import spectral as spx
import ssh_free as sf

DELTAS = [-0.30, +0.30]
VS = [0.0, 0.5, 1.0, 2.0, 3.0]


# ------------------------------------------------------------------ checks
def check_sum_rule():
    """GATE: sum_m w_jm = 1 - <n_j> at every site, both methods."""
    worst, where = 0.0, None
    for L in (8, 10, 12):
        for d in DELTAS:
            for V in VS:
                s = spx.addition_spectrum_full(L, d, V)
                for jj, j in enumerate(s.sites):
                    err = abs(s.weights[jj].sum() - (1.0 - s.n_expect[j]))
                    if err > worst:
                        worst, where = err, (L, d, V, j, "full")
    for L in (16,):
        for d in DELTAS:
            s = spx.addition_spectrum_lanczos(L, d, 1.0, 200, sites=(0, 2, L // 2))
            for jj, j in enumerate(s.sites):
                err = abs(s.weights[jj].sum() - (1.0 - s.n_expect[j]))
                if err > worst:
                    worst, where = err, (L, d, 1.0, j, "lanczos")
    assert worst < 1e-12, f"worst sum-rule error {worst:.2e} at {where}"
    return f"max |sum_m w_jm - (1 - <n_j>)| = {worst:.1e}  (all sites, both methods)"


def check_half_filling_normalisation():
    """1 - <n_j> = 1/2 exactly, so z_jm = 2 w_jm (Step 2's particle-hole result)."""
    worst = 0.0
    for L in (8, 12):
        for d in DELTAS:
            for V in (0.0, 1.0, 3.0):
                s = spx.addition_spectrum_full(L, d, V)
                worst = max(worst, float(np.max(np.abs(s.n_expect - 0.5))))
    assert worst < 1e-9, f"worst |<n_j> - 1/2| = {worst:.2e}"
    return f"max |<n_j> - 1/2| = {worst:.1e}, so the normalisation is exactly 1/2"


def check_free_limit_z0():
    """GATE: at V = 0, z_0 = 2|phi_+(0)|^2 = |psi_L(0)|^2 from Step 1."""
    worst, notes = 0.0, []
    for L in (8, 10, 12):
        d = -0.30
        s = spx.addition_spectrum_full(L, d, 0.0)
        ref = spx.free_reference(L, d)
        psiL0 = sf.localised_edge_states(L, d)["psi_left"][0] ** 2
        worst = max(worst, abs(s.z(0) - ref["z0"]), abs(s.z(0) - psiL0))
        notes.append(f"L={L}: {s.z(0):.9f}")
    assert worst < 1e-9, f"worst {worst:.2e}"
    r = (1 - 0.30) / (1 + 0.30)       # J_s/J_w inverted: r = J_w/J_s
    r = (1 + (-0.30)) / (1 - (-0.30))
    return (f"max deviation {worst:.1e}  [" + ", ".join(notes) +
            f"] -> semi-infinite 1-r^2 = {1 - r**2:.6f}")


def check_free_limit_omega():
    """GATE: at V = 0 the edge line sits at omega = eps, half Step 1's splitting."""
    worst, notes = 0.0, []
    for L in (8, 10, 12):
        d = -0.30
        s = spx.addition_spectrum_full(L, d, 0.0)
        om = np.min(s.omega_of(0))
        eps = 0.5 * sf.localised_edge_states(L, d)["splitting"]
        worst = max(worst, abs(om - eps))
        notes.append(f"L={L}: {om:.6e} vs {eps:.6e}")
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max deviation {worst:.1e}  [" + "; ".join(notes) + "]"


def check_free_limit_xi_W():
    """GATE: at V = 0, xi_W on the even sublattice reproduces Step 1's xi."""
    worst, notes = 0.0, []
    for L, m in ((12, None), (16, 200), (18, 200)):
        d = -0.30
        s = (spx.addition_spectrum_full(L, d, 0.0) if m is None
             else spx.addition_spectrum_lanczos(L, d, 0.0, m,
                                                sites=tuple(range(L // 2 + 1))))
        fit = spx.fit_weight_envelope(s, 0, n_fit=8)
        xi = sf.xi_analytic(d)
        rel = abs(fit["xi_W"] - xi) / xi
        worst = max(worst, rel)
        notes.append(f"L={L}: {fit['xi_W']:.5f} (r2={fit['r2']:.5f})")
        assert fit["r2"] > 0.999, f"L={L}: poor exponential fit, r2={fit['r2']}"
    assert worst < 0.02, f"worst relative deviation {worst:.2e}"
    return (f"xi = {sf.xi_analytic(-0.30):.5f};  " + "; ".join(notes) +
            f";  worst rel. dev. {worst:.1e}")


def check_factor_of_two():
    """GATE: the code fits W_j ~ exp(-2j/xi_W), not exp(-j/xi_W).

    A synthetic profile with a known xi is pushed through the real fitting
    routine; if the factor of two were dropped the answer would be out by 2x.
    """
    L, xi_true = 40, 5.0
    sites = tuple(range(L))
    W = np.exp(-2.0 * np.arange(L) / xi_true)
    W[1::2] = 0.0                                   # one sublattice only
    spec = spx.AdditionSpectrum(
        L=L, N=L // 2, delta=-0.30, V=0.0, method="synthetic", sites=sites,
        omegas=np.array([0.0]), weights=W[:, None],
        n_expect=np.full(L, 0.5), delta_plus_ring=1.0)
    fit = spx.fit_weight_envelope(spec, 0, n_fit=8)
    assert abs(fit["xi_W"] - xi_true) < 1e-9, \
        f"synthetic xi={xi_true} came back as {fit['xi_W']}"
    return (f"synthetic exp(-2j/{xi_true}) recovered as xi_W = {fit['xi_W']:.12f} "
            f"(a dropped factor of two would give {xi_true/2:g})")


def check_cluster_weight_rotation_invariance():
    """GATE: the cluster weight is invariant under rotations in a degenerate subspace.

    Two (near-)degenerate addition eigenstates are mixed by an arbitrary angle.
    The individual w_jm change; their sum must not.
    """
    L, d, V = 12, -0.30, 1.0
    N = L // 2
    secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
    g = mb.select_ground_state(mb.build_hamiltonian(secN, d, V), secN)
    vals, vecs = np.linalg.eigh(mb.build_hamiltonian(secP, d, V).toarray())
    phi = mb.cdag_matrix(L, N, 0) @ g.psi
    w_all = (vecs.T @ phi) ** 2
    gaps = np.diff(vals)
    # the closest pair that actually carries weight -- a pair with no weight would
    # make the test vacuous
    cand = [i for i in range(len(gaps)) if w_all[i] + w_all[i + 1] > 1e-4]
    p = min(cand, key=lambda i: gaps[i])
    base = float(((vecs[:, [p, p + 1]].T @ phi) ** 2).sum())
    worst_sum, max_individual_change = 0.0, 0.0
    w_ref = (vecs[:, [p, p + 1]].T @ phi) ** 2
    for theta in (0.3, 1.0, 2.2):
        c, s_ = np.cos(theta), np.sin(theta)
        rot = vecs[:, [p, p + 1]] @ np.array([[c, -s_], [s_, c]])
        w = (rot.T @ phi) ** 2
        worst_sum = max(worst_sum, abs(w.sum() - base))
        max_individual_change = max(max_individual_change,
                                    float(np.max(np.abs(w - w_ref))))
    assert worst_sum < 1e-12, f"cluster weight changed by {worst_sum:.2e}"
    assert max_individual_change > 1e-6, \
        "the rotation did nothing -- the test is not exercising anything"
    return (f"pair gap {gaps[p]:.2e}: cluster weight stable to {worst_sum:.1e} "
            f"while individual weights move by up to {max_individual_change:.2e}")


def check_lanczos_vs_full():
    """GATE: the Lanczos measure reproduces full ED at L = 12."""
    worst_z, worst_line = 0.0, 0.0
    for V in (0.0, 1.0, 3.0):
        full = spx.addition_spectrum_full(12, -0.30, V, sites=(0, 2, 4))
        for m in (100, 200, 400):
            lz = spx.addition_spectrum_lanczos(12, -0.30, V, m, sites=(0, 2, 4))
            for j in (0, 2, 4):
                worst_z = max(worst_z, abs(lz.z(j) - full.z(j)))
            worst_line = max(worst_line, abs(np.min(lz.omega_of(0))
                                             - np.min(full.omega_of(0))))
    assert worst_z < 1e-10 and worst_line < 1e-10
    return (f"max |z_lanczos - z_full| = {worst_z:.1e};  "
            f"max line-position error = {worst_line:.1e}  (m = 100, 200, 400)")


def check_lanczos_moments():
    """The Lanczos measure reproduces the first 2m-1 moments <v0|H^k|v0> exactly."""
    L, d, V, m = 10, -0.30, 1.0, 12
    N = L // 2
    secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
    g = mb.select_ground_state(mb.build_hamiltonian(secN, d, V), secN)
    HP = mb.build_hamiltonian(secP, d, V)
    v0 = mb.cdag_matrix(L, N, 0) @ g.psi
    theta, w, used = spx.lanczos_spectral_measure(HP, v0, m)
    worst, kmax = 0.0, 2 * used - 1
    x = v0.copy()
    for k in range(0, min(kmax, 15) + 1):
        exact = float(v0 @ x)
        approx = float(np.sum(w * theta**k))
        scale = max(abs(exact), 1.0)
        worst = max(worst, abs(exact - approx) / scale)
        x = HP @ x
    assert worst < 1e-10, f"worst relative moment error {worst:.2e}"
    return (f"moments k = 0..{min(kmax, 15)} reproduced to {worst:.1e} "
            f"(relative), with {used} Lanczos vectors")


def check_broadening_is_presentation_only():
    """The Lorentzian conserves total weight and changes no reported number."""
    s = spx.addition_spectrum_full(12, -0.30, 1.0, sites=(0,))
    om, w = s.omega_of(0), s.w_of(0)
    grid = np.linspace(om.min() - 60, om.max() + 60, 400001)
    worst_int, vals = 0.0, []
    for eta in spx.ETAS:
        curve = spx.lorentzian(grid, om, w, eta)
        worst_int = max(worst_int, abs(np.trapezoid(curve, grid) - w.sum()))
        vals.append((s.z(0), s.n_eff(0), spx.fit_weight_envelope(s, 0)["xi_W"]))
    assert worst_int < 1e-3, f"Lorentzian loses weight: {worst_int:.2e}"
    assert len(set(vals)) == 1, "a reported number depends on the broadening"
    return (f"integral matches the line weight to {worst_int:.1e}; "
            f"z_0, n_eff and xi_W identical for eta = {spx.ETAS}")


def check_n_eff_calibration():
    """n_eff = 1 for one line and 2 for two equal lines, by construction."""
    sites = (0,)
    for weights, expect in (([0.5], 1.0), ([0.25, 0.25], 2.0),
                            ([0.4, 0.1], (0.5**2) / (0.4**2 + 0.1**2))):
        W = np.array([weights])
        spec = spx.AdditionSpectrum(
            L=4, N=2, delta=-0.30, V=0.0, method="synthetic", sites=sites,
            omegas=np.linspace(0.01, 0.02, len(weights)), weights=W,
            n_expect=np.full(4, 0.5), delta_plus_ring=1.0)
        got = spec.n_eff(0)
        assert abs(got - expect) < 1e-12, f"{weights}: n_eff {got} != {expect}"
    return "n_eff = 1 (one line), 2 (two equal lines), 1.47 (0.4 / 0.1)"


def check_window_from_ring_threshold():
    """The window uses Delta+ measured on a ring, and at V=0 that is 2|delta|."""
    worst, notes = 0.0, []
    for L in (8, 12, 16):
        for d in DELTAS:
            dp = spx.ring_addition_threshold(L, d, 0.0)
            if (L // 2) % 2 == 0:            # allowed momenta include k = pi
                worst = max(worst, abs(dp - 2 * abs(d)))
            if L == 12:
                notes.append(f"d={d:+.2f}: {dp:.6f} vs 2|d|={2*abs(d):.2f}")
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max |Delta+ - 2|delta|| = {worst:.1e}  [L=12: " + "; ".join(notes) + "]"


def check_control_strong_end_bonds():
    """GATE: for strong end bonds there is no in-window line at all."""
    notes = []
    for L in (8, 12):
        for V in (0.0, 1.0, 3.0):
            s = spx.addition_spectrum_full(L, +0.30, V)
            assert s.n_lines(0) == 0, f"L={L} V={V}: unexpected in-window line"
            assert s.z(0) == 0.0
            lo = float(np.min(s.omega_of(0)))
            assert lo > 0.8 * s.delta_plus_ring
            if L == 12:
                notes.append(f"V={V:g}: lowest omega {lo:.4f} > "
                             f"{0.8*s.delta_plus_ring:.4f}")
    return "window empty at every point  [L=12: " + "; ".join(notes) + "]"


def check_edge_versus_centre():
    """For weak end bonds the window weight is an edge property, not a bulk one."""
    notes = []
    for V in (0.0, 1.0, 3.0):
        L = 12
        s = spx.addition_spectrum_full(L, -0.30, V)
        z_e, z_c = s.z(0), s.z(L // 2)
        assert z_e > 20 * max(z_c, 1e-12), \
            f"V={V}: edge {z_e:.4f} not dominant over centre {z_c:.4f}"
        notes.append(f"V={V:g}: {z_e:.4f} vs {z_c:.4f}")
    return "edge vs centre z  [" + "; ".join(notes) + "]"


def check_consistency_with_step2():
    """The lowest addition line is Step 2's open-chain addition gap."""
    worst = 0.0
    for L in (8, 12):
        for d in DELTAS:
            for V in (0.0, 1.0, 3.0):
                N = L // 2
                secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
                e0 = mb.select_ground_state(
                    mb.build_hamiltonian(secN, d, V), secN).energy
                e1 = mb.lowest_states(mb.build_hamiltonian(secP, d, V), k=1)[0][0]
                s = spx.addition_spectrum_full(L, d, V, sites=(0,))
                worst = max(worst, abs(float(np.min(s.omega_of(0))) - (e1 - e0)))
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max |omega_min - (E_0^(N+1) - E_0^N)| = {worst:.1e}"


def check_window_status_reporting():
    """GATE: an unstable window is labelled, not quietly averaged over."""
    s = spx.addition_spectrum_lanczos(18, -0.30, 3.0, 200,
                                      sites=tuple(range(10)))
    sc = spx.window_scan(s, 0)
    f_safe = spx.safe_factor(s, 0)
    assert sc["status"] == "unstable", f"expected 'unstable', got {sc['status']}"
    assert f_safe < 0.8, f"safe factor should be below 0.8, got {f_safe}"
    assert s.n_lines(0, f_safe) == 1 and s.n_lines(0, 0.9) > 1
    z_safe, z_default = s.z(0, f_safe), s.z(0, 0.8)
    assert abs(z_safe - z_default) > 1e-3, "the two windows should disagree here"
    return (f"L=18, V=3: status '{sc['status']}', safe factor {f_safe:g} "
            f"(1 line, z={z_safe:.4f}) vs 0.8 ({s.n_lines(0, 0.8)} lines, "
            f"z={z_default:.4f})")


def check_reproducible():
    """Two identical runs give identical numbers (no hidden randomness)."""
    a = spx.addition_spectrum_full(10, -0.30, 1.0)
    b = spx.addition_spectrum_full(10, -0.30, 1.0)
    assert a.z(0) == b.z(0) and a.sum_rule_error == b.sum_rule_error
    c = spx.addition_spectrum_lanczos(12, -0.30, 1.0, 150, sites=(0,))
    d_ = spx.addition_spectrum_lanczos(12, -0.30, 1.0, 150, sites=(0,))
    assert c.z(0) == d_.z(0)
    return "full ED and Lanczos both reproduce bit-for-bit"


CHECKS = [
    ("1  GATE  sum rule sum_m w_jm = 1-<n_j>", check_sum_rule),
    ("2  half filling: normalisation is exactly 1/2", check_half_filling_normalisation),
    ("3  GATE  V=0: z_0 = |psi_L(0)|^2 from Step 1", check_free_limit_z0),
    ("4  GATE  V=0: edge line sits at omega = eps", check_free_limit_omega),
    ("5  GATE  V=0: xi_W reproduces Step 1's xi", check_free_limit_xi_W),
    ("6  GATE  the fit carries the factor of two", check_factor_of_two),
    ("7  GATE  cluster weight is rotation invariant",
     check_cluster_weight_rotation_invariance),
    ("8  GATE  Lanczos measure vs full ED at L=12", check_lanczos_vs_full),
    ("9  Lanczos reproduces the moments exactly", check_lanczos_moments),
    ("10 broadening is presentation only", check_broadening_is_presentation_only),
    ("11 n_eff calibration", check_n_eff_calibration),
    ("12 window uses the ring threshold", check_window_from_ring_threshold),
    ("13 GATE  control: strong end bonds, empty window",
     check_control_strong_end_bonds),
    ("14 edge weight dominates the centre", check_edge_versus_centre),
    ("15 lowest line = Step 2's addition gap", check_consistency_with_step2),
    ("16 GATE  an unstable window is reported", check_window_status_reporting),
    ("17 runs are reproducible", check_reproducible),
]

for _label, _fn in CHECKS:
    globals()["test_" + _fn.__name__[6:]] = _fn


def main() -> int:
    print("Project A -- Step 3 gate checks\n" + "=" * 78)
    failed = 0
    for label, fn in CHECKS:
        try:
            print(f"[PASS] {label}\n       {fn()}")
        except AssertionError as exc:
            failed += 1
            print(f"[FAIL] {label}\n       {exc}")
    print("=" * 78)
    if failed:
        print(f"{failed} of {len(CHECKS)} checks FAILED -- do not go on to Step 4.")
    else:
        print(f"all {len(CHECKS)} checks passed -- Step 3 gates are clear for Step 4.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
