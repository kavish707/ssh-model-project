"""
Project A -- Step 2: the many-body occupation basis.

Conventions are inherited unchanged from Step 1 (ssh_free.py):

    sites      j = 0 … L-1,  L even,  N = L/2 at half filling
    bonds      J_j = 1 + (-1)^j delta,  so J_0 = 1 + delta
    H = - sum_b J_b (c^dag_j c_k + h.c.) + V sum_b (n_j - 1/2)(n_k - 1/2)

    open      bonds b = (j, j+1) for j = 0 … L-2
    periodic  plus the closing bond (L-1, 0) with J_{L-1} = 1 + (-1)^{L-1} delta,
              which for even L is J_s -- the same bond the ring hopping matrix of
              Step 1 adds, so the one- and many-body rings agree at V = 0.

An occupation state is the integer s = sum_j n_j 2^j, i.e. site j is bit j.

Fermionic signs (the guide's rule, implemented literally)
--------------------------------------------------------
Acting with c_j or c^dag_j contributes (-1)^{sum_{l<j} n_l}.  For c^dag_j c_k we
annihilate at k first, taking the sign from the *incoming* state, then create at
j, taking the second sign from the *intermediate* state.  Applied to the bond
(L-1, 0) of a ring this automatically produces the parity factor that an
open-chain spin bond would be missing -- see test_step2.check_ring_parity_factor.

Everything is sparse: the N sector has dimension C(L,N) = 12,870 at L = 16 and
48,620 at L = 18, and no dense matrix of that size is ever formed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

import ssh_free as sf

__all__ = [
    "fock_basis", "Sector", "bond_list", "build_hamiltonian", "number_operator",
    "cdag_matrix", "c_matrix", "density_profile", "correlation_matrix",
    "reflection_matrix", "lowest_states", "residual", "GroundState",
    "select_ground_state", "free_correlation_matrix", "DEGENERACY_TOL",
]

DEGENERACY_TOL = 1e-9       # states closer than this are treated as one cluster
SYMMETRY_TOL = 1e-6         # |<P>| must exceed 1 - this to count as resolved


# --------------------------------------------------------------------------
# 0. bit helpers
# --------------------------------------------------------------------------
def _popcount(x: NDArray[np.int64]) -> NDArray[np.int64]:
    return np.bitwise_count(x).astype(np.int64)


def _sign_below(states: NDArray[np.int64], j: int) -> NDArray[np.float64]:
    """(-1)^{sum_{l<j} n_l}, evaluated on whatever states are passed in."""
    return 1.0 - 2.0 * (_popcount(states & ((1 << j) - 1)) & 1)


def fock_basis(L: int, N: int) -> NDArray[np.int64]:
    """All integers with exactly N set bits below bit L, ascending (Gosper's hack)."""
    if N < 0 or N > L:
        return np.empty(0, dtype=np.int64)
    if N == 0:
        return np.zeros(1, dtype=np.int64)
    out, s, limit = [], (1 << N) - 1, 1 << L
    while s < limit:
        out.append(s)
        c = s & -s
        r = s + c
        s = (((r ^ s) >> 2) // c) | r
    return np.array(out, dtype=np.int64)


@dataclass
class Sector:
    """A fixed-particle-number sector: its basis and an index lookup."""
    L: int
    N: int
    states: NDArray[np.int64] = field(repr=False)

    @classmethod
    def make(cls, L: int, N: int) -> "Sector":
        return cls(L=L, N=N, states=fock_basis(L, N))

    @property
    def dim(self) -> int:
        return int(self.states.size)

    def index(self, s: NDArray[np.int64]) -> NDArray[np.int64]:
        """Positions of the given states; they are guaranteed to be present."""
        i = np.searchsorted(self.states, s)
        return i

    def occupations(self) -> NDArray[np.int64]:
        """dim x L array of n_j, used for the diagonal interaction and densities."""
        j = np.arange(self.L)
        return ((self.states[:, None] >> j[None, :]) & 1).astype(np.int64)


# --------------------------------------------------------------------------
# 1. bonds and the Hamiltonian
# --------------------------------------------------------------------------
def bond_list(L: int, delta: float, periodic: bool = False):
    """[(j, k, J_jk), ...] matching Step 1's hopping_matrix / ring_hopping_matrix."""
    J = sf.bond_strengths(L, delta)
    bonds = [(j, j + 1, float(J[j])) for j in range(L - 1)]
    if periodic:
        bonds.append((L - 1, 0, float(1.0 + ((-1.0) ** (L - 1)) * delta)))
    return bonds


def _hop_entries(sector: Sector, p: int, q: int):
    """Matrix entries of c^dag_p c_q within one sector, vectorised over the basis.

    Returns (rows, cols, signs).  The two signs are evaluated exactly as the guide
    prescribes: the first on the incoming state, the second on the state *after*
    the annihilation at q.
    """
    s = sector.states
    ok = (((s >> q) & 1) == 1) & (((s >> p) & 1) == 0)
    cols = np.flatnonzero(ok)
    if cols.size == 0:
        return cols, cols, np.empty(0)
    s_in = s[cols]
    sign = _sign_below(s_in, q)                 # annihilate at q, sign on |s>
    s_mid = s_in ^ (np.int64(1) << q)
    sign = sign * _sign_below(s_mid, p)         # create at p, sign on the intermediate
    s_out = s_mid | (np.int64(1) << p)
    rows = sector.index(s_out)
    return rows, cols, sign


def build_hamiltonian(sector: Sector, delta: float, V: float,
                      periodic: bool = False) -> sp.csr_matrix:
    """Sparse H in one fixed-N sector.  The interaction is diagonal."""
    L, dim = sector.L, sector.dim
    bonds = bond_list(L, delta, periodic)
    occ = sector.occupations() - 0.5

    diag = np.zeros(dim)
    rows_all, cols_all, vals_all = [], [], []
    for (j, k, J) in bonds:
        if V != 0.0:
            diag += V * occ[:, j] * occ[:, k]
        for (p, q) in ((j, k), (k, j)):          # c^dag_p c_q, both directions
            r, c, sgn = _hop_entries(sector, p, q)
            if r.size:
                rows_all.append(r)
                cols_all.append(c)
                vals_all.append(-J * sgn)

    H = sp.coo_matrix(
        (np.concatenate(vals_all), (np.concatenate(rows_all), np.concatenate(cols_all))),
        shape=(dim, dim),
    ).tocsr()
    H = H + sp.diags(diag, format="csr")
    H.sum_duplicates()
    return H


def number_operator(sector: Sector) -> NDArray[np.float64]:
    """Diagonal of the total number operator -- constant N by construction."""
    return sector.occupations().sum(axis=1).astype(float)


# --------------------------------------------------------------------------
# 2. operators between sectors
# --------------------------------------------------------------------------
def cdag_matrix(L: int, N: int, j: int) -> sp.csr_matrix:
    """c^dag_j as a (dim_{N+1} x dim_N) sparse matrix."""
    src, dst = Sector.make(L, N), Sector.make(L, N + 1)
    ok = ((src.states >> j) & 1) == 0
    cols = np.flatnonzero(ok)
    s_in = src.states[cols]
    sign = _sign_below(s_in, j)
    rows = dst.index(s_in | (np.int64(1) << j))
    return sp.coo_matrix((sign, (rows, cols)),
                         shape=(dst.dim, src.dim)).tocsr()


def c_matrix(L: int, N: int, j: int) -> sp.csr_matrix:
    """c_j as a (dim_{N-1} x dim_N) sparse matrix."""
    src, dst = Sector.make(L, N), Sector.make(L, N - 1)
    ok = ((src.states >> j) & 1) == 1
    cols = np.flatnonzero(ok)
    s_in = src.states[cols]
    sign = _sign_below(s_in, j)
    rows = dst.index(s_in ^ (np.int64(1) << j))
    return sp.coo_matrix((sign, (rows, cols)),
                         shape=(dst.dim, src.dim)).tocsr()


def density_profile(psi: NDArray[np.float64], sector: Sector) -> NDArray[np.float64]:
    """<n_j> for every site -- diagonal, so just a weighted sum of occupations."""
    return (sector.occupations() * (psi**2)[:, None]).sum(axis=0)


def correlation_matrix(psi: NDArray[np.float64], L: int, N: int) -> NDArray[np.float64]:
    """C[i,j] = <psi| c^dag_i c_j |psi> via phi_j = c_j psi in the N-1 sector."""
    phi = np.column_stack([c_matrix(L, N, j) @ psi for j in range(L)])
    return phi.T @ phi


def free_correlation_matrix(L: int, delta: float, N: int,
                            periodic: bool = False) -> NDArray[np.float64]:
    """V = 0 reference: C = sum over the N filled orbitals of phi_a(i) phi_a(j)."""
    _, vecs = sf.single_particle_spectrum(L, delta, periodic=periodic)
    occ = vecs[:, :N]
    return occ @ occ.T


# --------------------------------------------------------------------------
# 3. reflection symmetry (used to fix the state-selection convention)
# --------------------------------------------------------------------------
def _reverse_bits(states: NDArray[np.int64], L: int) -> NDArray[np.int64]:
    out = np.zeros_like(states)
    for j in range(L):
        out |= ((states >> j) & 1) << (L - 1 - j)
    return out


def reflection_phase(N: int) -> float:
    """(-1)^{N(N-1)/2}: the cost of restoring canonical operator order after
    reversing the chain.  Constant inside a fixed-N sector."""
    return -1.0 if ((N * (N - 1) // 2) & 1) else 1.0


def reflection_matrix(sector: Sector, fermionic: bool = True) -> sp.csr_matrix:
    """Reflection about the chain centre, j -> L-1-j.

    `fermionic=True` gives the true many-body operator R: reversing the
    occupation pattern leaves the creation operators in decreasing order, so
    restoring canonical order costs (-1)^{N(N-1)/2}.  `fermionic=False` gives bare
    bit reversal, i.e. R without that phase.

    Both commute with H on an *open* chain of even L, because J_j = J_{L-2-j};
    they differ only by a constant, so they have the same eigenvectors.  The bare
    eigenvalue is the one quoted as the state-selection label, because the
    fermionic phase flips with L (N = L/2) and would otherwise relabel the same
    physical state as L changes.
    """
    L, N = sector.L, sector.N
    phase = reflection_phase(N) if fermionic else 1.0
    rows = sector.index(_reverse_bits(sector.states, L))
    cols = np.arange(sector.dim)
    return sp.coo_matrix((np.full(sector.dim, phase), (rows, cols)),
                         shape=(sector.dim, sector.dim)).tocsr()


# --------------------------------------------------------------------------
# 4. eigensolving, residuals and the ground-state convention
# --------------------------------------------------------------------------
def lowest_states(H: sp.csr_matrix, k: int = 4, dense_max: int = 2000):
    """Lowest k eigenpairs, ascending.  Dense below dense_max, Lanczos above.

    The guide asks for the whole matrix at L = 6, 8, 10; with dense_max = 2000
    that is exactly what happens there (dim 20, 70, 252), and larger sizes use a
    sparse lowest-eigenpair solver.
    """
    dim = H.shape[0]
    if dim <= dense_max:
        vals, vecs = np.linalg.eigh(H.toarray())
        return vals[:k], vecs[:, :k]
    k = min(k, dim - 2)
    vals, vecs = spla.eigsh(H, k=k, which="SA", tol=0.0,
                            ncv=min(dim - 1, max(4 * k, 40)))
    order = np.argsort(vals)
    return vals[order], vecs[:, order]


def residual(H: sp.csr_matrix, E: float, psi: NDArray[np.float64]) -> float:
    """||H psi - E psi|| with psi normalised to one."""
    psi = psi / np.linalg.norm(psi)
    return float(np.linalg.norm(H @ psi - E * psi))


@dataclass
class GroundState:
    """The selected half-filled ground state and everything needed to reproduce it."""
    energy: float
    psi: NDArray[np.float64] = field(repr=False)
    energies: NDArray[np.float64] = field(repr=False)
    residual: float
    cluster_size: int
    cluster_spread: float
    gap_to_first_excited: float
    parity: float                 # <R>, the fermionic reflection eigenvalue
    parity_bare: float            # bit-reversal eigenvalue, the L-stable label
    convention: str
    rotated: bool                 # was a rotation inside a degenerate cluster needed?


def select_ground_state(H: sp.csr_matrix, sector: Sector, k: int = 4,
                        convention: str = "parity+", periodic: bool = False,
                        degeneracy_tol: float = DEGENERACY_TOL) -> GroundState:
    """Pick the half-filled ground state by an explicit, reproducible rule.

    The guide warns that silently accepting an arbitrary eigensolver vector from a
    near-degenerate pair can change the edge observables.  The rule used here:

    1.  Diagonalise and collect the cluster of states within `degeneracy_tol` of
        the lowest energy.
    2.  Since [H, R] = 0 on an open even-L chain, a *resolved* eigenvector is
        already an R eigenvector.  If every state in the cluster has |<R>| close
        to 1, nothing is ambiguous: take the lowest one and record its parity.
    3.  Only if the solver returned a mixture (|<R>| not close to 1) do we
        re-diagonalise R inside the cluster and pick the requested parity.  The
        rotated vector is an eigenvector of H only to within the cluster spread,
        which is recorded and reported.

    conventions: "parity+" / "parity-" pick that reflection eigenvalue;
                 "lowest" takes whatever the solver returned first -- kept so the
                 convention dependence can be measured rather than assumed.
    """
    vals, vecs = lowest_states(H, k=k)
    cluster = np.flatnonzero(vals - vals[0] <= degeneracy_tol)
    spread = float(vals[cluster].max() - vals[cluster].min())
    gap = float(vals[len(cluster)] - vals[0]) if len(cluster) < len(vals) else np.nan

    parity = np.nan
    rotated = False
    psi = vecs[:, 0]
    energy = float(vals[0])
    phase = reflection_phase(sector.N)

    if periodic or convention == "lowest":
        # the ring is reflection symmetric only about particular axes; we do not
        # use R there, and "lowest" deliberately takes whatever came back first
        Rb = reflection_matrix(sector, fermionic=False)
        parity_bare = float(psi @ (Rb @ psi))
        parity = phase * parity_bare
    else:
        Rb = reflection_matrix(sector, fermionic=False)
        sub = vecs[:, cluster]
        p_exp = np.array([v @ (Rb @ v) for v in sub.T])
        if np.all(np.abs(np.abs(p_exp) - 1.0) < SYMMETRY_TOL):
            parity_bare = float(p_exp[0])                  # already resolved
        else:
            rotated = True
            Rsub = sub.T @ (Rb @ sub)
            w, U = np.linalg.eigh(0.5 * (Rsub + Rsub.T))
            want = +1.0 if convention == "parity+" else -1.0
            pick = int(np.argmin(np.abs(w - want)))
            psi = sub @ U[:, pick]
            psi /= np.linalg.norm(psi)
            parity_bare = float(w[pick])
            energy = float(psi @ (H @ psi))
        parity = phase * parity_bare

    return GroundState(
        energy=energy,
        psi=psi,
        energies=vals,
        residual=residual(H, energy, psi),
        cluster_size=int(len(cluster)),
        cluster_spread=spread,
        gap_to_first_excited=gap,
        parity=parity,
        parity_bare=parity_bare,
        convention=convention,
        rotated=rotated,
    )
