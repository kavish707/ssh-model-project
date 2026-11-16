#!/usr/bin/env python3
"""
Project A -- Step 4 gate checks.

    python test_step4.py          # prints a table, exits 1 if anything fails
    pytest test_step4.py -q       # same checks under pytest

The guide's "before moving on" for Step 4: the spectral and time-evolution
return probabilities agree to the numerical tolerance on the small system, and
the interval free from visible finite-size returns is stated explicitly.  Checks
marked GATE cover those and the conservation laws they rest on.
"""

from __future__ import annotations

import sys

import numpy as np

import dynamics as dy
import manybody as mb
import ssh_free as sf

DELTAS = [-0.30, +0.30]
VS = [0.0, 0.5, 1.0, 2.0, 3.0]
BOUND = dy.TOL * dy.T_MAX          # accumulated Krylov error bound, 8e-10


def centre(L):
    return 2 * (L // 4)


_CACHE: dict = {}


def run(L, d, V, site=0, t_max=dy.T_MAX):
    """Memoised production Krylov run (several checks share the same runs)."""
    k = (L, d, V, site, t_max)
    if k not in _CACHE:
        _CACHE[k] = dy.krylov_evolve(dy.prepare(L, d, V, site), t_max=t_max)
    return _CACHE[k]


# ------------------------------------------------------------------ checks
def check_spectral_vs_time_evolution():
    """GATE: G(t) from Krylov equals sum_m z_0m exp(-i omega_m t) from Step 3."""
    worst_amp, worst_P, n = 0.0, 0.0, 0
    for L in (8, 12):
        for d in DELTAS:
            for V in VS:
                for site in (0, centre(L)):
                    k = dy.krylov_evolve(dy.prepare(L, d, V, site))
                    G, _ = dy.spectral_amplitude(L, d, V, site, k.times)
                    worst_amp = max(worst_amp, float(np.max(np.abs(k.amp - G))))
                    worst_P = max(worst_P, float(np.max(np.abs(k.P_return - abs(G)**2))))
                    n += 1
    assert worst_amp < BOUND and worst_P < BOUND, \
        f"amplitude {worst_amp:.2e}, P_return {worst_P:.2e} exceed {BOUND:.0e}"
    return (f"{n} runs at L = 8, 12: max |G_Krylov - G_spectral| = {worst_amp:.1e}, "
            f"max |dP_return| = {worst_P:.1e}  (bound {BOUND:.0e}, phase included)")


def check_krylov_vs_exact():
    """GATE: at L = 8, Krylov agrees with full diagonalisation (the guide's item 1)."""
    worst = {"amp": 0.0, "P_return": 0.0, "dn": 0.0}
    for d in DELTAS:
        for V in VS:
            for site in (0, centre(8)):
                p = dy.prepare(8, d, V, site)
                c = dy.compare(dy.krylov_evolve(p), dy.exact_evolve(p))
                for k_ in worst:
                    worst[k_] = max(worst[k_], c[k_])
    assert max(worst.values()) < BOUND, f"{worst}"
    return (f"max |dG| = {worst['amp']:.1e}, |dP| = {worst['P_return']:.1e}, "
            f"|d Delta n| = {worst['dn']:.1e}  (20 runs)")


def check_tolerance_halving():
    """GATE: halving the propagation tolerance changes nothing beyond the bound."""
    worst, notes = 0.0, []
    for L in (8, 12, 16):
        for V in (0.0, 1.0, 3.0):
            p = dy.prepare(L, -0.30, V)
            a = dy.krylov_evolve(p)
            b = dy.krylov_evolve(p, tol=dy.TOL / 2)
            c = dy.compare(a, b)
            worst = max(worst, c["amp"], c["dn"])
        notes.append(f"L={L}")
    assert worst < BOUND, f"halving changed results by {worst:.2e}"
    return f"max change on halving tol {dy.TOL:g} -> {dy.TOL/2:g}: {worst:.1e} ({', '.join(notes)})"


def check_tolerance_scan_has_teeth():
    """The tolerance is a real control: error falls with it and respects the bound."""
    p = dy.prepare(8, -0.30, 1.0)
    e = dy.exact_evolve(p)
    errs = []
    for tol in (1e-3, 1e-5, 1e-7, 1e-9):
        c = dy.compare(dy.krylov_evolve(p, tol=tol), e)
        errs.append(c["amp"])
        assert c["amp"] <= tol * dy.T_MAX, f"tol {tol}: error {c['amp']} above bound"
    assert all(b < a for a, b in zip(errs, errs[1:])), f"not decreasing: {errs}"
    return "error vs exact: " + ", ".join(
        f"tol {t:.0e} -> {e_:.1e}" for t, e_ in zip((1e-3, 1e-5, 1e-7, 1e-9), errs))


def check_conservation():
    """GATE: norm, energy and sum_j Delta n_j = 1 are conserved."""
    worst = {"norm": 0.0, "energy": 0.0, "sum": 0.0}
    for L in (8, 12, 16):
        for d in DELTAS:
            for V in (0.0, 3.0):
                for site in (0, centre(L)):
                    k = run(L, d, V, site)
                    worst["norm"] = max(worst["norm"], float(np.max(np.abs(k.norm - 1))))
                    worst["energy"] = max(worst["energy"],
                                          float(np.max(np.abs(k.energy - k.energy[0]))))
                    worst["sum"] = max(worst["sum"],
                                       float(np.max(np.abs(k.total_excess() - 1))))
    assert max(worst.values()) < 1e-10, f"{worst}"
    return (f"|norm-1| <= {worst['norm']:.1e}, |E(t)-E(0)| <= {worst['energy']:.1e}, "
            f"|sum Delta n - 1| <= {worst['sum']:.1e}")


def check_t0_identities():
    """Delta n_j(0) = 1-<n_j> at the injection site, -<dn dn>_c/(1-<n_j>) elsewhere."""
    worst = 0.0
    for L in (8, 12, 16):
        for d in DELTAS:
            for V in (0.0, 1.0, 3.0):
                for site in (0, centre(L)):
                    p = dy.prepare(L, d, V, site)
                    k = dy.krylov_evolve(p, t_max=0.05)
                    pred = -p.corr_ref / (1.0 - p.n_ref[site])
                    pred[site] = 1.0 - p.n_ref[site]
                    worst = max(worst, float(np.max(np.abs(k.dn[0] - pred))))
    assert worst < 1e-12, f"worst {worst:.2e}"
    return f"max deviation {worst:.1e} (the identity holds for any V)"


def check_free_fermion_route():
    """GATE: at V = 0 the many-body dynamics equals one-body propagation with
    Step 1's hopping matrix -- an implementation that never builds a Fock space."""
    worst_G, worst_dn = 0.0, 0.0
    for L in (8, 12, 16):
        for d in DELTAS:
            for site in (0, centre(L)):
                k = run(L, d, 0.0, site)
                G, dn = dy.free_fermion_trajectory(L, d, site, k.times)
                worst_G = max(worst_G, float(np.max(np.abs(k.amp - G))))
                worst_dn = max(worst_dn, float(np.max(np.abs(k.dn - dn))))
    assert worst_G < BOUND and worst_dn < BOUND
    return f"max |G_mb - G_1body| = {worst_G:.1e}, max |dn_mb - |chi|^2| = {worst_dn:.1e}"


def check_filtered_preparation():
    """GATE: the filtered excitation has preparation probability z_0(W) from
    Step 3, never decays, and carries exactly half its charge at each end."""
    notes = []
    for L in (8, 12, 16):
        for V in (0.0, 1.0, 3.0):
            p = dy.prepare(L, -0.30, V)
            f = dy.prepare_filtered(p)
            assert f is not None
            assert abs(f.prep_weight - f.window["z_step3"]) < 1e-10, \
                f"L={L} V={V}: {f.prep_weight} vs Step 3 {f.window['z_step3']}"
            k = dy.krylov_evolve(f)
            assert np.max(np.abs(k.P_return - 1.0)) < 1e-10, "filtered state decays"
            assert np.max(np.abs(k.dn - k.dn[0])) < 1e-10, "filtered density moves"
            left = float(k.dn[0, : L // 2].sum())
            assert abs(left - 0.5) < 1e-10, f"left half carries {left}, not 1/2"
            if L == 16:
                notes.append(f"V={V:g}: z0(W)={f.prep_weight:.4f}")
    for V in (0.0, 3.0):
        assert dy.prepare_filtered(dy.prepare(12, +0.30, V)) is None, \
            "strong end bonds should leave nothing to filter onto"
    return ("prep. probability = Step 3 z_0(W) to 1e-10; P_return = 1 and static "
            "density; exactly 1/2 per half  [L=16: " + "; ".join(notes) +
            "];  delta=+0.30: empty window")


def check_finite_size_criterion():
    """GATE: the window detects returns, not initial-state offsets.

    Three cases: a run compared with itself never departs; L = 8 vs 12 is flagged
    because the initial states already differ by more than eps (the far edge is
    only 2.5 xi away); L = 16 vs 18 at V = 3 has a transient blip just above eps
    at t ~ 1.85 that is the relaxation of a static offset, which must be ignored,
    and a genuine return later on, which must be found.
    """
    a = run(16, -0.30, 1.0)
    assert not dy.finite_size_time(a, a)["reached"], "a run departed from itself"
    r8 = dy.finite_size_time(run(8, -0.30, 1.0), run(12, -0.30, 1.0))
    assert not r8["initial_state_ok"] and r8["static_offset_max"] > dy_eps(), \
        f"L = 8 should be flagged: {r8}"
    a3, b3 = run(16, -0.30, 3.0), run(18, -0.30, 3.0)
    q_a, q_b = a3.Q_left(4), b3.Q_left(4)
    i = int(round(1.85 / dy.DT_OUT))
    blip = abs((q_a[i] - q_a[0]) - (q_b[i] - q_b[0]))
    r3 = dy.finite_size_time(a3, b3)
    assert blip > 0.01, "the transient this test is about has disappeared"
    assert r3["initial_state_ok"] and r3["reached"] and r3["t_fs"] > 3.0, f"{r3}"
    return (f"self: no departure;  L=8 vs 12: flagged, static offset "
            f"{r8['static_offset_max']:.3f};  L=16 vs 18, V=3: blip {blip:.4f} at "
            f"t=1.85 ignored (below eps + 2s), return found at t = {r3['t_fs']:.2f}")


def dy_eps():
    return 0.01


def check_finite_size_window_is_stated():
    """GATE: L = 16 gets an explicit certified interval against L = 18."""
    notes = []
    for V in (0.0, 1.0, 3.0):
        r = dy.finite_size_time(run(16, -0.30, V), run(18, -0.30, V))
        assert r["initial_state_ok"], f"V={V}: L = 16 initial state not converged"
        assert r["t_fs"] > 0
        notes.append(f"V={V:g}: " + (f"t < {r['t_fs']:.2f}" if r["reached"]
                                    else f"t <= {r['t_common_max']:g} (no return)"))
    return "L = 16 vs 18: " + "; ".join(notes)


def check_late_window_logic():
    """Late averages use only the second half of the certified interval."""
    import run_step4 as r4
    assert np.allclose(r4.late_window(5.3, 8.0, True), (2.65, 5.3))
    assert np.allclose(r4.late_window(np.inf, 8.0, False), (4.0, 8.0))
    a, b = r4.late_window(1.5, 8.0, True)
    assert np.isnan(a) and np.isnan(b), "a 1.5-long window is too short to average"
    return "(5.3) -> [2.65, 5.3];  uncertified-but-clean -> [4, 8];  1.5 -> none"


def check_no_splitting_beat():
    """No two-edge beat at the N-sector splitting in the selected preparation.

    A splitting-induced beat needs a second line at omega_edge +- (E_1^N - E_0^N)
    that the electron also reaches.  Testing "no pair of lines anywhere near the
    splitting" would be meaningless -- with ~50 lines some pair always lands
    near any frequency -- so the test is about the edge line's partners, plus
    the Fourier weight of Q_edge(t) at the splitting over an exact run to t=300.
    """
    notes = []
    for V in (0.0, 1.0, 3.0):
        p = dy.prepare(8, -0.30, V)
        vals, vecs = np.linalg.eigh(p.H.toarray())
        w = np.abs(vecs.T @ p.psi0) ** 2
        om = vals - p.E0_N
        ie = int(np.argmax(w))
        partner = sum(w[m] for m in range(om.size) if m != ie and
                      min(abs(om[m] - om[ie] - p.gs_gap),
                          abs(om[m] - om[ie] + p.gs_gap)) < 0.02)
        e = dy.exact_evolve(p, t_max=300.0)
        q = e.Q_left(4) - e.Q_left(4).mean()
        a_split = 2 * abs(np.mean(q * np.exp(1j * p.gs_gap * e.times)))
        freqs = np.linspace(0.05, 3.0, 1200)
        a_max = max(2 * abs(np.mean(q * np.exp(1j * f * e.times))) for f in freqs)
        assert partner < 1e-3 * w[ie], f"V={V}: partner weight {partner:.2e}"
        assert a_split < 0.05 * a_max, f"V={V}: Fourier weight at splitting {a_split}"
        notes.append(f"V={V:g}: partner weight {partner:.0e}, Fourier "
                     f"{a_split:.1e} vs max {a_max:.1e}")
    return "; ".join(notes)


def check_diagonal_ensemble_symmetric():
    """At finite L the infinite-time average puts exactly 1/2 of the excess charge
    in each half.  Degenerate N+1 eigenvalues (common at V = 0, where many-body
    levels are sums of one-body ones) are handled by projecting onto each
    eigenspace as a block."""
    worst, ndeg = 0.0, 0
    for V in (0.0, 1.0, 3.0):
        p = dy.prepare(8, -0.30, V)
        vals, vecs = np.linalg.eigh(p.H.toarray())
        occ = p.sector.occupations().astype(float)
        dens = np.zeros(p.L)
        i = 0
        while i < vals.size:
            j = i
            while j + 1 < vals.size and vals[j + 1] - vals[i] < 1e-9:
                j += 1
            ndeg += int(j > i)
            blk = vecs[:, i:j + 1]
            phi = blk @ (blk.T @ p.psi0)
            dens += (np.abs(phi) ** 2) @ occ
            i = j + 1
        worst = max(worst, abs((dens - p.n_ref)[: p.L // 2].sum() - 0.5))
    assert worst < 1e-10, f"worst {worst:.2e}"
    return (f"|Q_left-half(inf) - 1/2| <= {worst:.1e}  (L = 8, V = 0, 1, 3; "
            f"{ndeg} degenerate blocks handled)")


def check_far_edge_response():
    """The t = 0 far-edge excess approaches (1-r^2)^2/2 at V = 0, and vanishes for
    a symmetry-broken choice from the same doublet."""
    r = (1 - 0.30) / (1 + 0.30)
    lim = (1 - r**2) ** 2 / 2
    devs = []
    for L in (8, 12, 16):
        devs.append(abs(dy.far_edge_response(L, -0.30, 0.0)["selected"]["dn_far"] - lim))
    assert all(b < a for a, b in zip(devs, devs[1:])), f"not converging: {devs}"
    assert devs[-1] < 0.01
    f = dy.far_edge_response(16, -0.30, 0.0)
    assert f["left_loaded"]["dn_far"] < 0.01 and f["right_loaded"]["dn_far"] < 0.01
    return (f"|dn_far - {lim:.4f}| = " + ", ".join(f"{x:.4f}" for x in devs) +
            f" (L = 8, 12, 16);  broken: {f['left_loaded']['dn_far']:.1e}, "
            f"{f['right_loaded']['dn_far']:.1e}")


def check_late_return_tracks_z0_squared():
    """Interpretation check: inside the CERTIFIED window the return probability
    settles near z_0(W)^2 for edge injection, and near zero for the controls."""
    import run_step4 as r4
    notes = []
    for V in (0.0, 1.0, 3.0):
        k = run(16, -0.30, V)
        fs = dy.finite_size_time(k, run(18, -0.30, V))
        t_a, t_b = r4.late_window(fs["t_fs"], dy.T_MAX, fs["reached"])
        m = (k.times >= t_a) & (k.times <= t_b + 1e-9)
        z = dy.prepare_filtered(k.prep).prep_weight
        Pl = float(k.P_return[m].mean())
        assert abs(Pl - z**2) < 0.05, f"V={V}: <P>_late {Pl:.3f} vs z0^2 {z*z:.3f}"
        t = run(16, +0.30, V)
        assert t.P_return[m].mean() < 0.05
        notes.append(f"V={V:g} [{t_a:.2f},{t_b:.2f}]: {Pl:.3f} vs {z*z:.3f}")
    return "edge <P>_late vs z0^2  " + "; ".join(notes) + ";  delta=+0.30 < 0.05"


def check_reproducible():
    p = dy.prepare(12, -0.30, 1.0)
    a, b = dy.krylov_evolve(p), dy.krylov_evolve(p)
    assert np.array_equal(a.amp, b.amp) and np.array_equal(a.dn, b.dn)
    return "bit-for-bit identical on rerun"


CHECKS = [
    ("1  GATE  spectral = time-evolution amplitude", check_spectral_vs_time_evolution),
    ("2  GATE  Krylov = full diagonalisation, L=8", check_krylov_vs_exact),
    ("3  GATE  halving the tolerance", check_tolerance_halving),
    ("4  the tolerance scan has teeth", check_tolerance_scan_has_teeth),
    ("5  GATE  norm, energy, sum Delta n = 1", check_conservation),
    ("6  exact t = 0 identities", check_t0_identities),
    ("7  GATE  independent one-body route at V=0", check_free_fermion_route),
    ("8  GATE  filtered excitation", check_filtered_preparation),
    ("9  GATE  finite-size criterion sees returns", check_finite_size_criterion),
    ("10 GATE  finite-size window is stated", check_finite_size_window_is_stated),
    ("11 late averages stay in the window", check_late_window_logic),
    ("12 no splitting-induced beat", check_no_splitting_beat),
    ("13 diagonal ensemble splits charge 1/2 : 1/2", check_diagonal_ensemble_symmetric),
    ("14 far-edge response at t = 0", check_far_edge_response),
    ("15 late return tracks z_0^2", check_late_return_tracks_z0_squared),
    ("16 reproducible", check_reproducible),
]

for _label, _fn in CHECKS:
    globals()["test_" + _fn.__name__[6:]] = _fn


def main() -> int:
    print("Project A -- Step 4 gate checks\n" + "=" * 78, flush=True)
    failed = 0
    for label, fn in CHECKS:
        try:
            print(f"[PASS] {label}\n       {fn()}", flush=True)
        except AssertionError as exc:
            failed += 1
            print(f"[FAIL] {label}\n       {exc}", flush=True)
    print("=" * 78)
    if failed:
        print(f"{failed} of {len(CHECKS)} checks FAILED -- do not go on to Step 5.")
    else:
        print(f"all {len(CHECKS)} checks passed -- Step 4 gates are clear for Step 5.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
