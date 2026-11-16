"""
Project A -- Step 2, item 5: the Jordan-Wigner second implementation.

Nothing here reuses the occupation-basis matrix builder.  Pauli operators are
constructed directly over the full 2^L Fock space (basis state = the integer
s = sum_j n_j 2^j, so site j is bit j, exactly as in manybody.py), and the
Hamiltonian is assembled from Eq. (A.4):

    n_j = (1 - Z_j)/2,
    H  = -(1/2) sum_j J_j (X_j X_{j+1} + Y_j Y_{j+1}) + (V/4) sum_j Z_j Z_{j+1}

for OPEN boundaries, where the Jordan-Wigner strings of adjacent sites cancel.
Its fixed-number block is then compared with the occupation-basis matrix.

The same full-Fock-space machinery gives a second, independent check that the
occupation-basis Hamiltonian conserves particle number: we build H without any
sector restriction and evaluate [H, N].

    python jw_check.py
"""

from __future__ import annotations

import sys

import numpy as np
import scipy.sparse as sp

import manybody as mb

__all__ = ["pauli_x", "pauli_y", "pauli_z", "jw_hamiltonian",
           "fock_hamiltonian", "number_operator_full", "fixed_n_block"]


def _basis(L: int) -> np.ndarray:
    return np.arange(1 << L, dtype=np.int64)


def pauli_x(L: int, j: int) -> sp.csr_matrix:
    """X_j: |s> -> |s XOR 2^j>."""
    s = _basis(L)
    return sp.coo_matrix((np.ones(s.size), (s ^ (1 << j), s)),
                         shape=(s.size, s.size)).tocsr()


def pauli_y(L: int, j: int) -> sp.csr_matrix:
    """Y_j: Y|0> = i|1>, Y|1> = -i|0>."""
    s = _basis(L)
    bit = (s >> j) & 1
    vals = np.where(bit == 0, 1j, -1j)
    return sp.coo_matrix((vals, (s ^ (1 << j), s)),
                         shape=(s.size, s.size)).tocsr()


def pauli_z(L: int, j: int) -> sp.csr_matrix:
    """Z_j: +1 on n_j = 0, -1 on n_j = 1, consistent with n_j = (1 - Z_j)/2."""
    s = _basis(L)
    return sp.diags(1.0 - 2.0 * ((s >> j) & 1), format="csr")


def jw_hamiltonian(L: int, delta: float, V: float) -> sp.csr_matrix:
    """Eq. (A.4) over the full 2^L space, open boundaries."""
    bonds = mb.bond_list(L, delta, periodic=False)
    dim = 1 << L
    H = sp.csr_matrix((dim, dim), dtype=complex)
    for (j, k, J) in bonds:
        H = H - 0.5 * J * (pauli_x(L, j) @ pauli_x(L, k)
                           + pauli_y(L, j) @ pauli_y(L, k))
        if V != 0.0:
            H = H + 0.25 * V * (pauli_z(L, j) @ pauli_z(L, k))
    return H


def fock_hamiltonian(L: int, delta: float, V: float,
                     periodic: bool = False) -> sp.csr_matrix:
    """The occupation-basis H built over the FULL 2^L space (no sector restriction).

    Uses the same fermionic sign rule as manybody.build_hamiltonian, but with an
    unrestricted basis, so that [H, N] can be evaluated directly.
    """
    s = _basis(L)
    dim = s.size
    bonds = mb.bond_list(L, delta, periodic)
    occ = ((s[:, None] >> np.arange(L)[None, :]) & 1).astype(float) - 0.5

    diag = np.zeros(dim)
    rows, cols, vals = [], [], []
    for (j, k, J) in bonds:
        if V != 0.0:
            diag += V * occ[:, j] * occ[:, k]
        for (p, q) in ((j, k), (k, j)):
            ok = (((s >> q) & 1) == 1) & (((s >> p) & 1) == 0)
            c = np.flatnonzero(ok)
            s_in = s[c]
            sign = 1.0 - 2.0 * (np.bitwise_count(s_in & ((1 << q) - 1)) & 1)
            s_mid = s_in ^ (np.int64(1) << q)
            sign = sign * (1.0 - 2.0 * (np.bitwise_count(s_mid & ((1 << p) - 1)) & 1))
            rows.append(s_mid | (np.int64(1) << p))
            cols.append(c)
            vals.append(-J * sign)
    H = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(dim, dim)).tocsr()
    return H + sp.diags(diag, format="csr")


def number_operator_full(L: int) -> sp.csr_matrix:
    s = _basis(L)
    return sp.diags(np.bitwise_count(s).astype(float), format="csr")


def fixed_n_block(H_full: sp.csr_matrix, L: int, N: int) -> np.ndarray:
    """Restrict a full-Fock-space operator to the fixed-N block, in the same
    basis ordering manybody.Sector uses (ascending integers)."""
    idx = mb.fock_basis(L, N)
    return np.asarray(H_full.tocsr()[idx][:, idx].todense())


def main() -> int:
    print("Project A -- Step 2: Jordan-Wigner cross-check\n" + "=" * 78)
    fails = 0

    # --- 1. JW vs occupation basis at L = 6, as the guide asks
    L, N = 6, 3
    worst = 0.0
    for delta in (-0.30, +0.30):
        for V in (0.0, 0.5, 1.0, 2.0, 3.0):
            sector = mb.Sector.make(L, N)
            H_occ = np.asarray(mb.build_hamiltonian(sector, delta, V).todense())
            H_jw = fixed_n_block(jw_hamiltonian(L, delta, V), L, N)
            assert np.max(np.abs(H_jw.imag)) < 1e-13, "JW block is not real"
            worst = max(worst, float(np.max(np.abs(H_jw.real - H_occ))))
    ok = worst < 1e-12
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] L = 6 fixed-N block: Eq. (A.4) vs the "
          f"occupation basis\n       max |H_JW - H_occ| = {worst:.2e}  "
          f"(10 parameter points, dim 20)")

    # --- 2. same comparison on the whole spectrum, L = 6 and 8
    worst_spec = 0.0
    for L, N in ((6, 3), (8, 4)):
        for delta in (-0.30, +0.30):
            for V in (0.0, 1.0, 3.0):
                sector = mb.Sector.make(L, N)
                e_occ = np.linalg.eigvalsh(
                    np.asarray(mb.build_hamiltonian(sector, delta, V).todense()))
                e_jw = np.linalg.eigvalsh(fixed_n_block(
                    jw_hamiltonian(L, delta, V), L, N))
                worst_spec = max(worst_spec, float(np.max(np.abs(e_occ - e_jw))))
    ok = worst_spec < 1e-11
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] full fixed-N spectra agree at L = 6 and 8"
          f"\n       max |E_JW - E_occ| = {worst_spec:.2e}")

    # --- 3. number conservation, evaluated without assuming a sector
    worst_comm = 0.0
    for L in (6, 8):
        Nop = number_operator_full(L)
        for delta in (-0.30, +0.30):
            for V in (0.0, 2.0):
                for per in (False, True):
                    H = fock_hamiltonian(L, delta, V, periodic=per)
                    worst_comm = max(worst_comm,
                                     float(abs(H @ Nop - Nop @ H).max()))
    ok = worst_comm < 1e-13
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] [H, N] = 0 on the full 2^L space "
          f"(open and periodic)\n       max |[H, N]| = {worst_comm:.2e}")

    # --- 4. the sector matrix really is the fixed-N block of the full matrix
    worst_blk = 0.0
    for L, N in ((6, 3), (8, 4)):
        for delta in (-0.30, +0.30):
            for V in (0.0, 1.5):
                for per in (False, True):
                    sector = mb.Sector.make(L, N)
                    a = np.asarray(mb.build_hamiltonian(
                        sector, delta, V, periodic=per).todense())
                    b = fixed_n_block(fock_hamiltonian(L, delta, V, periodic=per), L, N)
                    worst_blk = max(worst_blk, float(np.max(np.abs(a - b))))
    ok = worst_blk < 1e-13
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] sector matrix == fixed-N block of the "
          f"unrestricted matrix\n       max deviation {worst_blk:.2e}")

    # --- 5. the ring's closing bond needs its parity factor
    L, N, delta, V = 8, 4, -0.30, 1.0
    sector = mb.Sector.make(L, N)
    H_right = mb.build_hamiltonian(sector, delta, V, periodic=True)
    H_wrong = _ring_without_parity(sector, delta, V)
    dE = abs(float(mb.lowest_states(H_right, 1)[0][0]
                   - mb.lowest_states(H_wrong, 1)[0][0]))
    ok = dE > 1e-3
    fails += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] dropping the parity factor on the closing "
          f"bond changes the ring\n       |E0(correct) - E0(naive spin bond)| = "
          f"{dE:.6f}  (it must not be zero)")

    print("=" * 78)
    print("Jordan-Wigner cross-checks all agree." if not fails
          else f"{fails} check(s) FAILED.")
    return 1 if fails else 0


def _ring_without_parity(sector, delta, V):
    """A deliberately WRONG ring: the closing bond treated as an open-chain spin
    bond, i.e. with the Jordan-Wigner string dropped.  Used only to show that the
    parity factor is not optional."""
    L = sector.L
    H = mb.build_hamiltonian(sector, delta, V, periodic=False)
    J = float(1.0 + ((-1.0) ** (L - 1)) * delta)
    s = sector.states
    rows, cols, vals = [], [], []
    for (p, q) in ((L - 1, 0), (0, L - 1)):
        ok = (((s >> q) & 1) == 1) & (((s >> p) & 1) == 0)
        c = np.flatnonzero(ok)
        s_out = (s[c] ^ (np.int64(1) << q)) | (np.int64(1) << p)
        rows.append(sector.index(s_out))
        cols.append(c)
        vals.append(np.full(c.size, -J))            # sign dropped on purpose
    occ = sector.occupations() - 0.5
    extra = sp.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(sector.dim, sector.dim)).tocsr()
    return H + extra + sp.diags(V * occ[:, L - 1] * occ[:, 0], format="csr")


if __name__ == "__main__":
    sys.exit(main())
