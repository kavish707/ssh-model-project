"""
Project A -- Step 3: localisation and coherence, measured separately.

Conventions inherited unchanged from Steps 1 and 2 (ssh_free.py, manybody.py).

The addition spectrum at site j (Eq. A.5):

    A+_j(omega) = sum_m w_jm delta(omega - omega_m),
    w_jm = |<m, N+1| c^dag_j |Omega_N>|^2,
    omega_m = E_m^{N+1} - E_0^N,

with the exact sum rule  sum_m w_jm = 1 - <n_j>, which holds operator-identically
because sum_m |m><m| = 1 on the N+1 sector and c_j c^dag_j = 1 - n_j.  At half
filling Step 2 established <n_j> = 1/2 on every site, so 1 - <n_j> = 1/2 and the
normalised weight is simply z_jm = 2 w_jm.  The code never assumes this: it
divides by the measured 1 - <n_j>.

Two quantities, two different questions:

    xi_W   how far the edge excitation extends  -- LOCALISATION
           from W_j = sum_{m in W} w_jm  fitted to  W_j ~ exp(-2j/xi_W)
    z_0    how much of one added electron goes into it -- COHERENCE
           z_0(W) = sum_{m in W} w_0m / (1 - <n_0>)

A third number is needed to say whether the window weight sits in ONE line or is
spread over several -- the question "does the electron operator efficiently create
one eigen-excitation".  We use the effective line count

    n_eff = (sum_m w_jm)^2 / sum_m w_jm^2        (= 1 for a single line)

taken over the window.  Like the cluster weight it is unchanged by rotations
inside an exactly degenerate subspace only when the whole subspace lies in the
window, which is why the window is always reported together with it.

At V = 0 everything here is checkable in closed form: the window holds the single
+eps orbital, and z_0 = |psi_L(0)|^2 -> 1 - r^2, xi_W = xi = 2/ln(J_s/J_w).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray
from scipy.linalg import eigh_tridiagonal

import manybody as mb
import ssh_free as sf

__all__ = [
    "AdditionSpectrum", "ring_addition_threshold", "addition_spectrum_full",
    "lanczos_spectral_measure", "addition_spectrum_lanczos", "lorentzian",
    "window_report", "window_scan", "fit_weight_envelope", "WEIGHT_FLOOR",
]

WEIGHT_FLOOR = 1e-12      # lines below this carry no usable weight
DEFAULT_FACTOR = 0.8      # the guide's starting edge window, 0 <= omega < 0.8 Delta+
FACTORS = (0.6, 0.7, 0.8, 0.9)
ETAS = (0.02, 0.01, 0.005)


# --------------------------------------------------------------------------
# 1. the spectrum object
# --------------------------------------------------------------------------
@dataclass
class AdditionSpectrum:
    L: int
    N: int
    delta: float
    V: float
    method: str                                   # "full" or "lanczos(m)"
    sites: tuple[int, ...]
    omegas: NDArray[np.float64] = field(repr=False)        # (M,) or (nsites, M)
    weights: NDArray[np.float64] = field(repr=False)       # (nsites, M)
    n_expect: NDArray[np.float64] = field(repr=False)
    E0_N: float = 0.0
    delta_plus_ring: float = np.nan
    addition_gap_open: float = np.nan
    sum_rule_error: float = np.nan
    gs_parity: float = np.nan
    gs_residual: float = np.nan
    gs_gap: float = np.nan
    gs_cluster: int = 1

    def _row(self, j: int) -> int:
        return self.sites.index(j)

    def omega_of(self, j: int) -> NDArray[np.float64]:
        """Line positions for site j (per-site for Lanczos, shared for full ED)."""
        return self.omegas[self._row(j)] if self.omegas.ndim == 2 else self.omegas

    def w_of(self, j: int) -> NDArray[np.float64]:
        return self.weights[self._row(j)]

    def window_mask(self, j: int, factor: float = DEFAULT_FACTOR):
        om = self.omega_of(j)
        return (om >= 0.0) & (om < factor * self.delta_plus_ring)

    def cluster_weight(self, j: int, factor: float = DEFAULT_FACTOR) -> float:
        """W_j = sum_{m in W} w_jm, the unnormalised window weight."""
        return float(self.w_of(j)[self.window_mask(j, factor)].sum())

    def z(self, j: int, factor: float = DEFAULT_FACTOR) -> float:
        """z_j(W) = W_j / (1 - <n_j>)."""
        return self.cluster_weight(j, factor) / (1.0 - self.n_expect[j])

    def n_eff(self, j: int, factor: float = DEFAULT_FACTOR) -> float:
        """Effective number of lines sharing the window weight at site j."""
        w = self.w_of(j)[self.window_mask(j, factor)]
        w = w[w > WEIGHT_FLOOR]
        if w.size == 0:
            return np.nan
        return float(w.sum() ** 2 / np.sum(w**2))

    def n_lines(self, j: int, factor: float = DEFAULT_FACTOR) -> int:
        w = self.w_of(j)[self.window_mask(j, factor)]
        return int(np.sum(w > WEIGHT_FLOOR))

    def isolation_gap(self, j: int = 0, factor: float = DEFAULT_FACTOR) -> float:
        """Distance in omega between the last line inside the window and the first
        line with usable weight outside it.  Large means the window is isolated."""
        om, w = self.omega_of(j), self.w_of(j)
        inside = self.window_mask(j, factor) & (w > WEIGHT_FLOOR)
        outside = (~self.window_mask(j, factor)) & (om >= 0) & (w > WEIGHT_FLOOR)
        if not inside.any() or not outside.any():
            return np.nan
        return float(om[outside].min() - om[inside].max())

    def weight_profile(self, factor: float = DEFAULT_FACTOR) -> NDArray[np.float64]:
        """W_j for every site in `sites`, in site order."""
        return np.array([self.cluster_weight(j, factor) for j in self.sites])


# --------------------------------------------------------------------------
# 2. the bulk threshold and the full spectrum
# --------------------------------------------------------------------------
def ring_addition_threshold(L: int, delta: float, V: float) -> float:
    """Delta+ = E_0,ring^{N+1} - E_0,ring^N, measured on a ring at the SAME parameters.

    This is the bulk addition threshold the edge window is defined against; at
    V = 0 it reduces to the 2|delta| of Step 1 whenever the ring's allowed momenta
    include k = pi (i.e. L/2 even).
    """
    N = L // 2
    secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
    e0 = mb.lowest_states(mb.build_hamiltonian(secN, delta, V, periodic=True), k=1)[0][0]
    e1 = mb.lowest_states(mb.build_hamiltonian(secP, delta, V, periodic=True), k=1)[0][0]
    return float(e1 - e0)


def addition_spectrum_full(L: int, delta: float, V: float,
                           sites: tuple[int, ...] | None = None,
                           convention: str = "parity+") -> AdditionSpectrum:
    """Full spectrum in the N+1 sector (the guide's item 1, for L <= 12)."""
    N = L // 2
    secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
    HN = mb.build_hamiltonian(secN, delta, V)
    g = mb.select_ground_state(HN, secN, convention=convention)
    HP = mb.build_hamiltonian(secP, delta, V)
    vals, vecs = np.linalg.eigh(HP.toarray())

    n = mb.density_profile(g.psi, secN)
    sites = tuple(range(L)) if sites is None else tuple(sites)
    W = np.empty((len(sites), vals.size))
    sum_err = 0.0
    for row, j in enumerate(sites):
        phi = mb.cdag_matrix(L, N, j) @ g.psi
        W[row] = (vecs.T @ phi) ** 2
        sum_err = max(sum_err, abs(W[row].sum() - (1.0 - n[j])))

    return AdditionSpectrum(
        L=L, N=N, delta=delta, V=V, method="full", sites=sites,
        omegas=vals - g.energy, weights=W, n_expect=n, E0_N=g.energy,
        delta_plus_ring=ring_addition_threshold(L, delta, V),
        addition_gap_open=float(vals[0] - g.energy),
        sum_rule_error=sum_err, gs_parity=g.parity_bare, gs_residual=g.residual,
        gs_gap=g.gap_to_first_excited, gs_cluster=g.cluster_size,
    )


# --------------------------------------------------------------------------
# 3. the Lanczos spectral measure (item 5)
# --------------------------------------------------------------------------
def lanczos_spectral_measure(H: sp.csr_matrix, v0: NDArray[np.float64], m: int,
                             reorth: bool = True, beta_floor: float = 1e-12):
    """Approximate lines and weights of <v0| delta(omega - H) |v0>.

    Builds the Krylov basis from v0, with FULL reorthogonalisation against every
    stored vector (twice, which is what makes the small weights trustworthy), then
    diagonalises the tridiagonal matrix.  The eigenvalues are the approximate
    lines; the squared first components of its eigenvectors, times ||v0||^2, are
    the approximate weights.

    Exactness: the first 2m-1 moments <v0|H^k|v0> are reproduced exactly, and the
    total weight equals ||v0||^2 by construction.  Neither fact means an
    individual line has converged -- only comparison against full ED, or stability
    in m, can say that.
    """
    nrm = float(np.linalg.norm(v0))
    if nrm == 0.0:
        return np.empty(0), np.empty(0), 0
    dim = v0.size
    m = min(m, dim)
    Q = np.empty((m, dim))                       # preallocated: reorth stays O(m^2 dim)
    Q[0] = v0 / nrm
    alpha: list[float] = []
    beta: list[float] = []
    for i in range(m):
        w = H @ Q[i]
        a = float(Q[i] @ w)
        alpha.append(a)
        w = w - a * Q[i] - (beta[-1] * Q[i - 1] if i > 0 else 0.0)
        if reorth:
            Qi = Q[: i + 1]
            for _ in range(2):
                w -= Qi.T @ (Qi @ w)
        b = float(np.linalg.norm(w))
        if b < beta_floor or i == m - 1:
            break
        beta.append(b)
        Q[i + 1] = w / b
    theta, U = eigh_tridiagonal(np.array(alpha), np.array(beta))
    return theta, nrm**2 * U[0, :] ** 2, len(alpha)


def addition_spectrum_lanczos(L: int, delta: float, V: float, m: int,
                              sites: tuple[int, ...] | None = None,
                              convention: str = "parity+") -> AdditionSpectrum:
    """Item 5: the Lanczos spectral measure with initial vector c^dag_j |Omega_N>."""
    N = L // 2
    secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
    HN = mb.build_hamiltonian(secN, delta, V)
    g = mb.select_ground_state(HN, secN, convention=convention)
    HP = mb.build_hamiltonian(secP, delta, V)
    n = mb.density_profile(g.psi, secN)
    sites = tuple(range(L)) if sites is None else tuple(sites)

    oms, wts, sum_err, used = [], [], 0.0, []
    for j in sites:
        phi = mb.cdag_matrix(L, N, j) @ g.psi
        theta, w, k = lanczos_spectral_measure(HP, phi, m)
        oms.append(theta - g.energy)
        wts.append(w)
        used.append(k)
        sum_err = max(sum_err, abs(w.sum() - (1.0 - n[j])))

    width = max(len(o) for o in oms)
    O = np.full((len(sites), width), np.inf)
    W = np.zeros((len(sites), width))
    for i, (o, w) in enumerate(zip(oms, wts)):
        O[i, : o.size], W[i, : w.size] = o, w

    return AdditionSpectrum(
        L=L, N=N, delta=delta, V=V, method=f"lanczos({min(used)}-{max(used)})",
        sites=sites, omegas=O, weights=W, n_expect=n, E0_N=g.energy,
        delta_plus_ring=ring_addition_threshold(L, delta, V),
        addition_gap_open=float(np.min([o.min() for o in oms])),
        sum_rule_error=sum_err, gs_parity=g.parity_bare, gs_residual=g.residual,
        gs_gap=g.gap_to_first_excited, gs_cluster=g.cluster_size,
    )


# --------------------------------------------------------------------------
# 4. presentation and window analysis
# --------------------------------------------------------------------------
def lorentzian(grid: NDArray[np.float64], omegas: NDArray[np.float64],
               weights: NDArray[np.float64], eta: float) -> NDArray[np.float64]:
    """Lines convolved with a Lorentzian of HWHM eta.

    This is PRESENTATION ONLY.  eta is a drawing choice, not a lifetime: a finite
    system has discrete lines and nothing here broadens them physically.  Every
    number reported in Step 3 is computed from the lines, never from this curve.
    """
    ok = np.isfinite(omegas) & (weights > 0.0)
    om, w = omegas[ok], weights[ok]
    out = np.zeros_like(grid, dtype=float)
    # chunked so that a fine grid times many lines cannot blow up memory
    chunk = max(1, int(4_000_000 // max(grid.size, 1)))
    for s in range(0, om.size, chunk):
        d = grid[:, None] - om[None, s:s + chunk]
        out += (w[None, s:s + chunk] * eta / np.pi / (d**2 + eta**2)).sum(axis=1)
    return out


def window_report(spec: AdditionSpectrum, j: int = 0,
                  factor: float = DEFAULT_FACTOR) -> dict:
    """Everything about the edge window at one site."""
    return {
        "site": j, "factor": factor,
        "window_hi": factor * spec.delta_plus_ring,
        "delta_plus_ring": spec.delta_plus_ring,
        "n_lines": spec.n_lines(j, factor),
        "W_j": spec.cluster_weight(j, factor),
        "z": spec.z(j, factor),
        "n_eff": spec.n_eff(j, factor),
        "isolation_gap": spec.isolation_gap(j, factor),
        "omega_lowest": float(np.min(spec.omega_of(j))),
    }


def window_scan(spec: AdditionSpectrum, j: int = 0,
                factors: tuple[float, ...] = FACTORS) -> dict:
    """Vary the window factor over 0.6-0.9 and report whether anything moves.

    The window counts as isolated when the cluster weight and the line count are
    unchanged across the range -- i.e. the boundary sits in a spectral gap wherever
    it is put.  If it is not, that is reported, not smoothed over.
    """
    zs = np.array([spec.z(j, f) for f in factors])
    ns = np.array([spec.n_lines(j, f) for f in factors])
    spread = float(zs.max() - zs.min())
    rel = spread / max(abs(zs.mean()), 1e-15)
    stable = bool(ns.min() == ns.max())
    if ns.max() == 0:
        status = "empty"          # no line below the bulk threshold at all
    elif stable and rel < 1e-6:
        status = "isolated"       # the boundary sits in a gap wherever it is put
    else:
        status = "unstable"       # the answer depends on where the boundary is
    return {
        "factors": factors, "z": zs, "n_lines": ns,
        "z_spread": spread, "z_rel_spread": rel,
        "line_count_stable": stable,
        "status": status,
        "isolated": status == "isolated",
    }


def safe_factor(spec: AdditionSpectrum, j: int = 0,
                factors: tuple[float, ...] = FACTORS) -> float:
    """Largest window factor at which the line count still equals its value at the
    smallest factor -- i.e. the widest window that has not yet swallowed another line.

    The guide fixes 0.8 as the starting factor and asks that it be varied over
    0.6-0.9.  Where 0.8 turns out not to sit in a spectral gap, this is the factor
    at which a defensible number can still be quoted; both are reported, and the
    difference between them is the systematic uncertainty from the window choice.
    """
    base = spec.n_lines(j, factors[0])
    best = factors[0]
    for f in factors:
        if spec.n_lines(j, f) == base:
            best = f
        else:
            break
    return best


def fit_weight_envelope(spec: AdditionSpectrum, sublattice: int = 0,
                        factor: float = DEFAULT_FACTOR, n_fit: int = 6,
                        max_site_frac: float = 0.5, floor: float = 1e-14,
                        from_right: bool = False) -> dict:
    """Fit the edge envelope of W_j to W_j ~ exp(-2 j / xi_W).

    Follows Step 1's fit discipline: one sublattice at a time, the window kept to
    j < L/2 so the far edge cannot contaminate it, and the fit range and the
    omitted weight both reported rather than hidden.

    Note the factor of two: W_j is a weight (a probability), so the slope of
    ln W_j against j is -2/xi_W, not -1/xi_W.

    With `from_right`, distance is measured from the far end (j -> L-1-j).  This
    is what the *other* sublattice needs: at delta < 0 the window weight sits on
    even sites near j = 0 and on odd sites near j = L-1, because the orbital being
    filled is the antibonding combination of the two end states.  Fitting the odd
    sublattice from the left correctly reports a non-decaying envelope.
    """
    L = spec.L
    W = spec.weight_profile(factor)
    sites = np.array(spec.sites)
    dist = (L - 1 - sites) if from_right else sites
    sel = (sites % 2 == sublattice) & (dist < max_site_frac * L) & (W > floor)
    order = np.argsort(dist[sel])
    js, ws = dist[sel][order], W[sel][order]
    sub_total = float(W[(sites % 2 == sublattice)].sum())
    out = {"sublattice": sublattice, "factor": factor, "from_right": from_right,
           "n_available": int(js.size),
           "weight_this_sublattice": sub_total,
           "weight_all_sites": float(W.sum())}
    if js.size < 3:
        out.update({"xi_W": np.nan, "r2": np.nan, "ok": False,
                    "reason": f"only {js.size} usable sites on this sublattice"})
        return out
    n_use = min(n_fit, js.size)
    js, ws = js[:n_use], ws[:n_use]
    slope, intercept = np.polyfit(js, np.log(ws), 1)
    resid = np.log(ws) - (slope * js + intercept)
    ss_tot = float(np.sum((np.log(ws) - np.log(ws).mean()) ** 2))
    out.update({
        "xi_W": float(-2.0 / slope) if slope < 0 else np.nan,
        "slope": float(slope), "intercept": float(intercept),
        "fit_site_min": int(js[0]), "fit_site_max": int(js[-1]), "n_fit": int(n_use),
        "r2": float(1.0 - np.sum(resid**2) / ss_tot) if ss_tot > 0 else np.nan,
        "max_resid": float(np.max(np.abs(resid))),
        "weight_in_fit": float(ws.sum()),
        # omitted relative to this sublattice, and relative to the whole window
        "weight_omitted_frac_sublattice":
            float(1.0 - ws.sum() / sub_total) if sub_total > 0 else np.nan,
        "weight_omitted_frac_total":
            float(1.0 - ws.sum() / W.sum()) if W.sum() > 0 else np.nan,
        "ok": bool(slope < 0),
        "reason": "" if slope < 0 else "envelope does not decay on this sublattice",
    })
    return out


def free_reference(L: int, delta: float) -> dict:
    """The V = 0 closed-form expectations, for cross-checking.

    The window holds the single +eps orbital, so
        omega = eps,  z_0 = |psi_L(0)|^2 -> 1 - r^2,  xi_W = xi.
    """
    ev, evec = sf.single_particle_spectrum(L, delta)
    i0, i1 = sf.edge_pair_indices(ev)
    r = (1.0 + delta) / (1.0 - delta)
    return {
        "omega_edge": float(ev[i1]),
        "z0": float(2.0 * evec[0, i1] ** 2),
        "z0_semi_infinite": float(1.0 - r**2),
        "xi": sf.xi_analytic(delta),
    }
