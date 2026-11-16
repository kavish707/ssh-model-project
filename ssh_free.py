"""
Project A -- Step 1: the free (V = 0) SSH edge problem.

Conventions (fixed by the guide, Eq. A.1):
    sites          j = 0, ..., L-1          (L even, N = L/2 particles)
    bonds          J_j = 1 + (-1)^j * delta,   j = 0, ..., L-2
    hopping matrix h[j, j+1] = h[j+1, j] = -J_j
    hbar = 1, lattice spacing = 1, J = 1.

delta < 0  =>  J_0 = 1 + delta is the WEAK end bond, so
               J_w = 1 + delta,  J_s = 1 - delta,  r = J_w/J_s < 1.

Semi-infinite zero mode on the even sublattice (Eq. A.2):
    J_w psi_{2l} + J_s psi_{2l+2} = 0   =>   psi_{2l} = sqrt(1-r^2) (-r)^l
Amplitude decay  ~ exp(-j/xi),  probability decay ~ exp(-2j/xi), with (Eq. A.3)
    xi = 2 / ln(J_s / J_w) ~= 1/|delta|.

Ring (two sites per cell):
    E_pm(k) = +- sqrt(J_w^2 + J_s^2 + 2 J_w J_s cos k)
    band gap            = 4|delta|
    bulk addition cost  = 2|delta|   (half filling, V = 0)

Everything here is pure NumPy/SciPy and deterministic: no random seeds are used.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "bond_strengths",
    "hopping_matrix",
    "weak_strong",
    "xi_analytic",
    "single_particle_spectrum",
    "edge_pair_indices",
    "localised_edge_states",
    "sublattice_profile",
    "fit_localisation_length",
    "analytic_edge_profile",
    "bloch_bands",
    "ring_hopping_matrix",
    "ring_gaps",
    "ground_state_energy",
]

TOL_DEGENERATE = 1e-12


# --------------------------------------------------------------------------
# 1. Bonds and the one-particle Hamiltonian
# --------------------------------------------------------------------------
def bond_strengths(L: int, delta: float) -> NDArray[np.float64]:
    """J_j = 1 + (-1)^j delta for the L-1 open-chain bonds j = 0..L-2."""
    if L % 2 != 0:
        raise ValueError("L must be even (half filling N = L/2).")
    j = np.arange(L - 1)
    return 1.0 + ((-1.0) ** j) * delta


def weak_strong(delta: float) -> tuple[float, float]:
    """(J_w, J_s) with the guide's convention for delta < 0: J_w = 1+delta, J_s = 1-delta.

    For delta > 0 the same expressions are returned, so that J_0 = 1 + delta is
    always the *first* bond; it is then the strong one and there is no edge mode
    for this termination.
    """
    return 1.0 + delta, 1.0 - delta


def hopping_matrix(L: int, delta: float) -> NDArray[np.float64]:
    """Real symmetric L x L matrix with h[j,j+1] = h[j+1,j] = -J_j."""
    J = bond_strengths(L, delta)
    h = np.zeros((L, L), dtype=float)
    idx = np.arange(L - 1)
    h[idx, idx + 1] = -J
    h[idx + 1, idx] = -J
    return h


def ring_hopping_matrix(L: int, delta: float) -> NDArray[np.float64]:
    """Periodic chain: the open matrix plus the closing bond J_{L-1} = 1 + (-1)^{L-1} delta.

    For even L the closing bond is J_{L-1} = 1 - delta = J_s, so the ring has an
    unbroken ... J_w J_s J_w J_s ... pattern and no edges.  (At V = 0 and one
    particle the fermionic sign on the closing bond is irrelevant; it matters
    only in the many-body problem of Step 2.)
    """
    h = hopping_matrix(L, delta)
    J_close = 1.0 + ((-1.0) ** (L - 1)) * delta
    h[L - 1, 0] = h[0, L - 1] = -J_close
    return h


def xi_analytic(delta: float) -> float:
    """xi = 2 / ln(J_s / J_w) with J_w = 1+delta, J_s = 1-delta.

    Positive (a genuine decay length) only for delta < 0, i.e. weak end bonds.
    At delta = -0.30 this returns 3.230813126 -- the value quoted in the guide.
    """
    Jw, Js = weak_strong(delta)
    if Jw <= 0.0 or Js <= 0.0:
        raise ValueError("|delta| < 1 required.")
    return 2.0 / np.log(Js / Jw)


# --------------------------------------------------------------------------
# 2. Spectrum and the edge doublet
# --------------------------------------------------------------------------
def single_particle_spectrum(L: int, delta: float, periodic: bool = False):
    """Return (energies ascending, eigenvectors as columns)."""
    h = ring_hopping_matrix(L, delta) if periodic else hopping_matrix(L, delta)
    evals, evecs = np.linalg.eigh(h)
    return evals, evecs


def edge_pair_indices(evals: NDArray[np.float64]) -> tuple[int, int]:
    """Indices of the two eigenvalues closest to zero, returned in ascending energy."""
    order = np.argsort(np.abs(evals))[:2]
    return tuple(sorted(order.tolist()))  # type: ignore[return-value]


def localised_edge_states(L: int, delta: float):
    """Disentangle the hybridised near-zero doublet into left/right localised states.

    On a finite chain the two ends hybridise, so the eigensolver returns the
    bonding/antibonding combinations at energies -eps and +eps.  Any unitary
    rotation inside that 2-d subspace is still an eigenbasis of H to within eps,
    so we fix the rotation by diagonalising the position operator X = diag(j)
    restricted to the subspace.  The eigenvector with the smaller <x> is the
    left-localised state.

    Returns a dict with:
        energies      the two eigenvalues (ascending)
        splitting     E_+ - E_-  (= 2 eps)
        psi_left      left-localised state, L-vector, unit norm
        psi_right     right-localised state
        x_left/right  <x> of each
        prob_sum      |psi_1|^2 + |psi_2|^2 of the raw doublet (the guide's
                      alternative, rotation invariant, construction)
    """
    evals, evecs = single_particle_spectrum(L, delta)
    i0, i1 = edge_pair_indices(evals)
    e_pair = evals[[i0, i1]]
    V = evecs[:, [i0, i1]]                      # L x 2

    x = np.arange(L, dtype=float)
    Xsub = V.T @ (x[:, None] * V)               # 2 x 2, symmetric
    xvals, U = np.linalg.eigh(Xsub)
    W = V @ U                                   # columns ordered by <x> ascending

    psi_left = W[:, 0] / np.linalg.norm(W[:, 0])
    psi_right = W[:, 1] / np.linalg.norm(W[:, 1])
    # fix an overall sign: make the largest-magnitude component positive
    for p in (psi_left, psi_right):
        p *= np.sign(p[np.argmax(np.abs(p))])

    return {
        "energies": e_pair,
        "splitting": float(e_pair[1] - e_pair[0]),
        "eps": float(0.5 * (e_pair[1] - e_pair[0])),
        "psi_left": psi_left,
        "psi_right": psi_right,
        "x_left": float(xvals[0]),
        "x_right": float(xvals[1]),
        "prob_sum": (V[:, 0] ** 2 + V[:, 1] ** 2),
        "all_energies": evals,
    }


def sublattice_profile(psi: NDArray[np.float64], sublattice: int = 0):
    """(sites, amplitudes) on the even (0) or odd (1) sublattice."""
    sites = np.arange(sublattice, psi.size, 2)
    return sites, psi[sites]


# --------------------------------------------------------------------------
# 3. Localisation-length fit
# --------------------------------------------------------------------------
def fit_localisation_length(
    psi: NDArray[np.float64],
    n_fit: int = 8,
    sublattice: int = 0,
    skip: int = 0,
    max_site_frac: float = 0.5,
    floor: float = 1e-13,
):
    """Fit the first `n_fit` nonzero sublattice weights to an exponential.

    The guide asks for the first 6-10 nonzero even-site weights, *excluding the
    opposite edge*: we therefore keep only sites with j < max_site_frac * L and
    amplitudes above `floor`.

    Two fits are performed and both are returned, to make the factor of two in
    A.3 explicit and checkable:
        amplitude     ln|psi_j|   vs j   has slope -1/xi
        probability   ln psi_j^2  vs j   has slope -2/xi

    Returns a dict with xi_amp, xi_prob, the fit window, residuals and R^2.
    """
    L = psi.size
    sites, amps = sublattice_profile(psi, sublattice)
    keep = (np.abs(amps) > floor) & (sites < max_site_frac * L)
    sites, amps = sites[keep], amps[keep]
    sites, amps = sites[skip:], amps[skip:]
    if sites.size < n_fit:
        raise ValueError(f"only {sites.size} usable points, need {n_fit}")
    sites, amps = sites[:n_fit], amps[:n_fit]

    out: dict[str, object] = {
        "fit_sites": sites,
        "fit_amplitudes": amps,
        "n_fit": int(n_fit),
        "sublattice": int(sublattice),
    }
    for label, y, factor in (
        ("amp", np.log(np.abs(amps)), 1.0),
        ("prob", np.log(amps**2), 2.0),
    ):
        slope, intercept = np.polyfit(sites, y, 1)
        resid = y - (slope * sites + intercept)
        ss_tot = np.sum((y - y.mean()) ** 2)
        out[f"xi_{label}"] = float(-factor / slope)
        out[f"slope_{label}"] = float(slope)
        out[f"intercept_{label}"] = float(intercept)
        out[f"max_resid_{label}"] = float(np.max(np.abs(resid)))
        out[f"r2_{label}"] = float(1.0 - np.sum(resid**2) / ss_tot) if ss_tot > 0 else np.nan
    return out


def analytic_edge_profile(L: int, delta: float) -> NDArray[np.float64]:
    """Semi-infinite zero mode psi_{2l} = sqrt(1-r^2)(-r)^l on even sites, zero on odd."""
    Jw, Js = weak_strong(delta)
    r = Jw / Js
    psi = np.zeros(L)
    ell = np.arange((L + 1) // 2)
    psi[2 * ell] = np.sqrt(max(1.0 - r**2, 0.0)) * (-r) ** ell
    return psi


# --------------------------------------------------------------------------
# 4. Ring / Bloch bands
# --------------------------------------------------------------------------
def bloch_bands(delta: float, k: NDArray[np.float64]):
    """E_pm(k) = +- sqrt(J_w^2 + J_s^2 + 2 J_w J_s cos k)."""
    Jw, Js = weak_strong(delta)
    E = np.sqrt(Jw**2 + Js**2 + 2.0 * Jw * Js * np.cos(k))
    return -E, E


def ring_gaps(L: int, delta: float):
    """Numerical ring gaps compared with the analytic predictions 4|delta|, 2|delta|.

    At V = 0 and half filling the N = L/2 lowest levels are filled, so
        addition cost   = E_{N}          (lowest empty level, index N)
        removal cost    = -E_{N-1}       (highest filled level)
        band gap        = E_{N} - E_{N-1}
    """
    evals, _ = single_particle_spectrum(L, delta, periodic=True)
    N = L // 2
    e_homo, e_lumo = evals[N - 1], evals[N]
    kgrid = np.linspace(-np.pi, np.pi, 20001)
    Em, Ep = bloch_bands(delta, kgrid)
    return {
        "L": L,
        "delta": delta,
        "band_gap_numeric": float(e_lumo - e_homo),
        "band_gap_analytic": 4.0 * abs(delta),
        "addition_cost_numeric": float(e_lumo),
        "addition_cost_analytic": 2.0 * abs(delta),
        "removal_cost_numeric": float(-e_homo),
        "band_gap_from_dispersion": float(Ep.min() - Em.max()),
        "bloch_max_abs_error": float(
            _bloch_vs_exact_error(L, delta, evals)
        ),
    }


def _bloch_vs_exact_error(L, delta, evals):
    """Max |exact ring level - nearest allowed Bloch level|."""
    ncell = L // 2
    k = 2.0 * np.pi * np.arange(ncell) / ncell
    Em, Ep = bloch_bands(delta, k)
    predicted = np.sort(np.concatenate([Em, Ep]))
    return np.max(np.abs(np.sort(evals) - predicted))


def ground_state_energy(L: int, delta: float, periodic: bool = False) -> float:
    """Sum of the lowest N = L/2 one-particle energies (the V = 0 half-filled GS)."""
    evals, _ = single_particle_spectrum(L, delta, periodic=periodic)
    return float(np.sum(evals[: L // 2]))
