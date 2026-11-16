#!/usr/bin/env python3
"""
Project A -- Step 2 gate checks.

    python test_step2.py          # prints a table, exits 1 if anything fails
    pytest test_step2.py -q       # same checks under pytest

The guide's "before moving on" for Step 2 is: record the lowest few energies,
residuals, boundary condition and state symmetry, and have an explicit
state-selection convention for a near-degenerate ground-state pair.  Checks
marked GATE cover those; the rest are the construction checks of items 1-5.
"""

from __future__ import annotations

import sys

import numpy as np
import scipy.sparse as sp

import jw_check as jw
import manybody as mb
import ssh_free as sf

L_SMALL = [6, 8, 10]
L_PROD = [12, 16]
DELTAS = [-0.30, +0.30]
VS = [0.0, 0.5, 1.0, 2.0, 3.0]
RES_TOL_SMALL = 1e-10
RES_TOL_PROD = 1e-9


# ------------------------------------------------------------------ checks
def check_dimensions():
    """The sector dimension is C(L,N): 12,870 at L=16 and 48,620 at L=18."""
    from math import comb
    for L in (6, 8, 10, 12, 16, 18):
        for N in (L // 2 - 1, L // 2, L // 2 + 1):
            assert mb.Sector.make(L, N).dim == comb(L, N), (L, N)
    d16, d18 = mb.Sector.make(16, 8).dim, mb.Sector.make(18, 9).dim
    assert (d16, d18) == (12870, 48620), (d16, d18)
    return f"L=16: {d16};  L=18: {d18}  (both match the guide)"


def check_basis_number_conservation():
    """Every basis state has exactly N bits, and every hop lands back in the sector."""
    for L in (6, 8, 10, 12):
        for N in (L // 2, L // 2 + 1):
            sec = mb.Sector.make(L, N)
            assert np.all(np.bitwise_count(sec.states) == N)
            assert np.all(np.diff(sec.states) > 0), "basis not strictly ascending"
            for (j, k, _) in mb.bond_list(L, -0.30, periodic=True):
                for (p, q) in ((j, k), (k, j)):
                    rows, cols, sgn = mb._hop_entries(sec, p, q)
                    if rows.size:
                        assert np.all(np.bitwise_count(sec.states[rows]) == N)
                        assert np.all(np.abs(sgn) == 1.0)
    return "all basis states and all hop targets carry exactly N particles"


def check_hermiticity():
    """GATE: H is symmetric in every sector, open and periodic."""
    worst = 0.0
    for L in L_SMALL + L_PROD:
        for N in (L // 2 - 1, L // 2, L // 2 + 1):
            for d in DELTAS:
                for V in (0.0, 1.0, 3.0):
                    for per in (False, True):
                        H = mb.build_hamiltonian(mb.Sector.make(L, N), d, V,
                                                 periodic=per)
                        worst = max(worst, float(abs(H - H.T).max()))
    assert worst == 0.0, f"max |H - H^T| = {worst:.2e}"
    return f"max |H - H^T| = {worst:.1e} over all sectors, both boundary conditions"


def check_interaction_is_diagonal():
    """The V term is diagonal, and the hopping part does not depend on V."""
    worst_off = worst_hop = 0.0
    for L in (6, 8, 10):
        sec = mb.Sector.make(L, L // 2)
        for d in DELTAS:
            for per in (False, True):
                H0 = mb.build_hamiltonian(sec, d, 0.0, periodic=per)
                for V in (0.5, 2.0):
                    HV = mb.build_hamiltonian(sec, d, V, periodic=per)
                    D = (HV - H0).tocoo()
                    off = D.data[D.row != D.col]
                    worst_off = max(worst_off, float(np.abs(off).max()) if off.size else 0.0)
                    hop = (HV - sp.diags(HV.diagonal())) - (H0 - sp.diags(H0.diagonal()))
                    worst_hop = max(worst_hop, float(abs(hop).max()))
    assert worst_off < 1e-14 and worst_hop < 1e-14
    return (f"max off-diagonal of H(V)-H(0) = {worst_off:.1e};  "
            f"hopping part V-independent to {worst_hop:.1e}")


def check_hopping_amplitudes():
    """Off-diagonal elements are exactly -J_b times a sign, and only between
    configurations differing by one particle moved across a bond."""
    worst = 0.0
    for L in (6, 8):
        sec = mb.Sector.make(L, L // 2)
        for d in DELTAS:
            for per in (False, True):
                bonds = mb.bond_list(L, d, per)
                allowed = {round(abs(J), 12) for (_, _, J) in bonds}
                H = mb.build_hamiltonian(sec, d, 0.0, periodic=per).tocoo()
                off = H.data[H.row != H.col]
                for v in off:
                    best = min(abs(abs(v) - a) for a in allowed)
                    worst = max(worst, best)
                # every off-diagonal pair differs by exactly two flipped bits
                diff = sec.states[H.row[H.row != H.col]] ^ sec.states[H.col[H.row != H.col]]
                assert np.all(np.bitwise_count(diff) == 2)
    assert worst < 1e-14, f"an off-diagonal element is not a bond amplitude ({worst:.2e})"
    return f"every off-diagonal element equals a bond J to {worst:.1e}; all differ by one hop"


def check_anticommutators():
    """{c_i, c^dag_j} = delta_ij on the N sector -- the sign rule, independently."""
    worst = 0.0
    for L in (6, 8):
        N = L // 2
        dim = mb.Sector.make(L, N).dim
        for i in range(L):
            for j in range(L):
                A = mb.c_matrix(L, N + 1, i) @ mb.cdag_matrix(L, N, j)
                B = mb.cdag_matrix(L, N - 1, j) @ mb.c_matrix(L, N, i)
                M = (A + B).toarray()
                target = np.eye(dim) if i == j else np.zeros((dim, dim))
                worst = max(worst, float(np.max(np.abs(M - target))))
    assert worst < 1e-13, f"anticommutator violated by {worst:.2e}"
    return f"max |{{c_i, c^dag_j}} - delta_ij| = {worst:.1e}  (L = 6, 8, all i, j)"


def check_residuals_small():
    """GATE: ||H psi - E psi|| < 1e-10 at L = 6, 8, 10 by full diagonalisation."""
    worst, where = 0.0, None
    for L in L_SMALL:
        for d in DELTAS:
            for V in VS:
                for per in (False, True):
                    sec = mb.Sector.make(L, L // 2)
                    H = mb.build_hamiltonian(sec, d, V, periodic=per)
                    vals, vecs = np.linalg.eigh(H.toarray())
                    for n in range(4):
                        r = mb.residual(H, vals[n], vecs[:, n])
                        if r > worst:
                            worst, where = r, (L, d, V, "periodic" if per else "open", n)
    assert worst < RES_TOL_SMALL, f"worst residual {worst:.2e} at {where}"
    return f"worst residual over 4 lowest states = {worst:.1e}  (need < {RES_TOL_SMALL:.0e})"


def check_residuals_production():
    """GATE: ||H psi - E psi|| < 1e-9 for the production (sparse-solver) states."""
    worst, where = 0.0, None
    for L in L_PROD + [18]:
        for d in DELTAS:
            for V in (0.0, 1.0, 3.0):
                for N in (L // 2 - 1, L // 2, L // 2 + 1):
                    sec = mb.Sector.make(L, N)
                    H = mb.build_hamiltonian(sec, d, V)
                    vals, vecs = mb.lowest_states(H, k=2)
                    for n in range(2):
                        r = mb.residual(H, vals[n], vecs[:, n])
                        if r > worst:
                            worst, where = r, (L, N, d, V, n)
    assert worst < RES_TOL_PROD, f"worst residual {worst:.2e} at {where}"
    return (f"worst residual = {worst:.1e} (need < {RES_TOL_PROD:.0e}); "
            f"largest sector dim {mb.Sector.make(18, 9).dim}")


def check_expectation_equals_eigenvalue():
    """GATE: an independently evaluated <psi|H|psi> reproduces the eigenvalue."""
    worst = 0.0
    for L in L_SMALL + L_PROD:
        for d in DELTAS:
            for V in (0.0, 1.0, 3.0):
                sec = mb.Sector.make(L, L // 2)
                H = mb.build_hamiltonian(sec, d, V)
                vals, vecs = mb.lowest_states(H, k=2)
                psi = vecs[:, 0] / np.linalg.norm(vecs[:, 0])
                worst = max(worst, abs(float(psi @ (H @ psi)) - vals[0]))
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max |<psi|H|psi> - E| = {worst:.1e}"


def check_free_limit_energy():
    """GATE: at V = 0 the many-body energy is the sum of the lowest N free levels."""
    worst, where = 0.0, None
    for L in L_SMALL + L_PROD:
        for d in DELTAS:
            for per in (False, True):
                sec = mb.Sector.make(L, L // 2)
                H = mb.build_hamiltonian(sec, d, 0.0, periodic=per)
                E0 = mb.lowest_states(H, k=1)[0][0]
                free = sf.ground_state_energy(L, d, periodic=per)
                if abs(E0 - free) > worst:
                    worst, where = abs(E0 - free), (L, d, per)
    assert worst < 1e-10, f"worst {worst:.2e} at {where}"
    return f"max |E_ED - sum of lowest N one-particle levels| = {worst:.1e}"


def check_free_limit_correlator():
    """GATE: at V = 0, <c^dag_i c_j> equals the filled-orbital correlation matrix."""
    worst = 0.0
    for L in L_SMALL + [12]:
        for d in DELTAS:
            for per in (False, True):
                N = L // 2
                sec = mb.Sector.make(L, N)
                H = mb.build_hamiltonian(sec, d, 0.0, periodic=per)
                vals, vecs = mb.lowest_states(H, k=2)
                if vals[1] - vals[0] < 1e-8:
                    continue           # degenerate: the correlator is basis dependent
                C = mb.correlation_matrix(vecs[:, 0], L, N)
                Cf = mb.free_correlation_matrix(L, d, N, periodic=per)
                worst = max(worst, float(np.max(np.abs(C - Cf))))
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max |C_ED - C_free| = {worst:.1e}  (degenerate ring points skipped)"


def check_jordan_wigner():
    """GATE: Eq. (A.4) reproduces the occupation-basis matrix at L = 6."""
    L, N = 6, 3
    worst = 0.0
    for d in DELTAS:
        for V in VS:
            H_occ = np.asarray(mb.build_hamiltonian(mb.Sector.make(L, N), d, V).todense())
            H_jw = jw.fixed_n_block(jw.jw_hamiltonian(L, d, V), L, N)
            assert np.max(np.abs(H_jw.imag)) < 1e-13
            worst = max(worst, float(np.max(np.abs(H_jw.real - H_occ))))
    assert worst < 1e-12, f"worst {worst:.2e}"
    return f"max |H_JW - H_occ| = {worst:.1e} at L = 6 over {2*len(VS)} points"


def check_number_conservation_full_space():
    """[H, N] = 0, evaluated on the unrestricted 2^L space (open and periodic)."""
    worst = 0.0
    for L in (6, 8):
        Nop = jw.number_operator_full(L)
        for d in DELTAS:
            for V in (0.0, 2.0):
                for per in (False, True):
                    H = jw.fock_hamiltonian(L, d, V, periodic=per)
                    worst = max(worst, float(abs(H @ Nop - Nop @ H).max()))
    assert worst < 1e-13, f"worst {worst:.2e}"
    return f"max |[H, N]| = {worst:.1e}"


def check_ring_parity_factor():
    """GATE: the closing bond needs its parity factor; a naive spin bond is wrong.

    Also checks the correct ring against Step 1's ring_hopping_matrix at V = 0, so
    that the many-body and one-particle rings are demonstrably the same model.
    """
    notes = []
    for L in (8, 12):
        sec = mb.Sector.make(L, L // 2)
        E_right = mb.lowest_states(
            mb.build_hamiltonian(sec, -0.30, 0.0, periodic=True), k=1)[0][0]
        free = sf.ground_state_energy(L, -0.30, periodic=True)
        assert abs(E_right - free) < 1e-10, \
            f"L={L}: correct ring disagrees with Step 1 by {abs(E_right-free):.2e}"
        E_wrong = mb.lowest_states(jw._ring_without_parity(sec, -0.30, 0.0), k=1)[0][0]
        assert abs(E_wrong - free) > 1e-3, \
            f"L={L}: dropping the parity factor made no difference -- suspicious"
        notes.append(f"L={L}: correct {E_right:.9f} == Step 1 ring; "
                     f"no-parity {E_wrong:.9f}")
    return " | ".join(notes)


def check_gap_matches_step1():
    """At V = 0 the many-body neutral gap IS the one-particle end-state splitting."""
    worst, notes = 0.0, []
    for L in (8, 12, 16):
        sec = mb.Sector.make(L, L // 2)
        H = mb.build_hamiltonian(sec, -0.30, 0.0)
        vals, _ = mb.lowest_states(H, k=2)
        gap = float(vals[1] - vals[0])
        split = sf.localised_edge_states(L, -0.30)["splitting"]
        worst = max(worst, abs(gap - split))
        notes.append(f"L={L}: {gap:.6e} vs {split:.6e}")
    assert worst < 1e-10, f"worst {worst:.2e}"
    return f"max deviation {worst:.1e}  [" + "; ".join(notes) + "]"


def check_state_selection_is_deterministic():
    """GATE: the selection rule is reproducible and reports what it did."""
    notes = []
    for L in (8, 12, 16):
        sec = mb.Sector.make(L, L // 2)
        for V in (0.0, 3.0):
            H = mb.build_hamiltonian(sec, -0.30, V)
            a = mb.select_ground_state(H, sec, convention="parity+")
            b = mb.select_ground_state(H, sec, convention="parity+")
            assert abs(a.energy - b.energy) < 1e-13
            assert abs(abs(float(a.psi @ b.psi)) - 1.0) < 1e-10, "not reproducible"
            assert abs(abs(a.parity_bare) - 1.0) < 1e-6, \
                f"parity not resolved: {a.parity_bare}"
            assert a.parity_bare > 0, "convention says the + eigenvector"
            assert a.residual < RES_TOL_PROD
            if V == 0.0:
                notes.append(f"L={L}: cluster {a.cluster_size}, parity {a.parity_bare:+.0f}, "
                             f"gap {a.gap_to_first_excited:.2e}, rotated={a.rotated}")
    return " | ".join(notes)


def check_degenerate_case_is_handled():
    """GATE: on a genuinely degenerate state the rule rotates and says so.

    The uniform (delta = 0) periodic chain at L divisible by four has the free
    half-filled zero-mode degeneracy the guide mentions in Step 6; it is the one
    place in reach where the eigensolver really does return an arbitrary vector.
    """
    notes = []
    for L in (8, 12):
        sec = mb.Sector.make(L, L // 2)
        H = mb.build_hamiltonian(sec, 0.0, 0.0, periodic=True)
        g = mb.select_ground_state(H, sec, k=4, convention="parity+", periodic=True)
        assert g.cluster_size > 1, f"L={L}: expected a degenerate cluster"
        assert g.cluster_spread < mb.DEGENERACY_TOL
        notes.append(f"L={L}: cluster {g.cluster_size}, spread {g.cluster_spread:.1e}")
    return " | ".join(notes) + "  (detected, not silently accepted)"


def check_half_filling_density():
    """Particle-hole symmetry at half filling: <n_j> = 1/2 on every site.

    This is why the ground-state density profile is NOT an edge diagnostic; the
    edge shows up in the addition spectrum of Step 3 instead.
    """
    worst = 0.0
    for L in (8, 12, 16):
        sec = mb.Sector.make(L, L // 2)
        for d in DELTAS:
            for V in (0.0, 1.0, 3.0):
                H = mb.build_hamiltonian(sec, d, V)
                g = mb.select_ground_state(H, sec)
                n = mb.density_profile(g.psi, sec)
                worst = max(worst, float(np.max(np.abs(n - 0.5))))
    assert worst < 1e-9, f"worst |<n_j> - 1/2| = {worst:.2e}"
    return f"max |<n_j> - 1/2| = {worst:.1e} over the grid"


def check_sector_independence():
    """The N+1 and N-1 matrices are built independently and are consistent with
    the N sector through c^dag: <N+1,0| c^dag_j |N,0> is nonzero somewhere."""
    L, N, d, V = 12, 6, -0.30, 1.0
    secN = mb.Sector.make(L, N)
    HN = mb.build_hamiltonian(secN, d, V)
    gN = mb.select_ground_state(HN, secN)
    secP = mb.Sector.make(L, N + 1)
    HP = mb.build_hamiltonian(secP, d, V)
    valsP, vecsP = mb.lowest_states(HP, k=2)
    weights = [float(vecsP[:, 0] @ (mb.cdag_matrix(L, N, j) @ gN.psi)) ** 2
               for j in range(L)]
    total = sum(weights)
    assert total > 1e-6, "c^dag |Omega> has no overlap with the N+1 ground state"
    assert secP.dim != secN.dim
    return (f"L={L}: dim(N)={secN.dim}, dim(N+1)={secP.dim}; "
            f"sum_j |<N+1,0|c^dag_j|Omega>|^2 = {total:.4f}  (Step 3 uses these)")


CHECKS = [
    ("1  sector dimensions C(L,N)", check_dimensions),
    ("2  basis and hops conserve N", check_basis_number_conservation),
    ("3  GATE  Hermiticity, all sectors", check_hermiticity),
    ("4  interaction is diagonal, hopping V-independent", check_interaction_is_diagonal),
    ("5  off-diagonals are bond amplitudes", check_hopping_amplitudes),
    ("6  anticommutators {c_i, c^dag_j} = delta_ij", check_anticommutators),
    ("7  GATE  residuals < 1e-10 at L = 6, 8, 10", check_residuals_small),
    ("8  GATE  residuals < 1e-9 for production states", check_residuals_production),
    ("9  GATE  <psi|H|psi> equals the eigenvalue", check_expectation_equals_eigenvalue),
    ("10 GATE  V=0 energy = sum of free levels", check_free_limit_energy),
    ("11 GATE  V=0 correlator = filled-orbital matrix", check_free_limit_correlator),
    ("12 GATE  Jordan-Wigner agrees at L = 6", check_jordan_wigner),
    ("13 [H, N] = 0 on the full 2^L space", check_number_conservation_full_space),
    ("14 GATE  ring closing-bond parity factor", check_ring_parity_factor),
    ("15 V=0 neutral gap = Step 1 end splitting", check_gap_matches_step1),
    ("16 GATE  state selection is deterministic", check_state_selection_is_deterministic),
    ("17 GATE  genuine degeneracy is detected", check_degenerate_case_is_handled),
    ("18 half filling: <n_j> = 1/2 everywhere", check_half_filling_density),
    ("19 N+1 sector consistent with c^dag on |Omega>", check_sector_independence),
]

for _label, _fn in CHECKS:
    globals()["test_" + _fn.__name__[6:]] = _fn


def main() -> int:
    print("Project A -- Step 2 gate checks\n" + "=" * 78)
    failed = 0
    for label, fn in CHECKS:
        try:
            print(f"[PASS] {label}\n       {fn()}")
        except AssertionError as exc:
            failed += 1
            print(f"[FAIL] {label}\n       {exc}")
    print("=" * 78)
    if failed:
        print(f"{failed} of {len(CHECKS)} checks FAILED -- do not go on to Step 3.")
    else:
        print(f"all {len(CHECKS)} checks passed -- Step 2 gates are clear for Step 3.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
