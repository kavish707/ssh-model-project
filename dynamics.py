"""
Project A -- Step 4: watch an edge particle evolve.

Conventions inherited unchanged from Steps 1-3 (ssh_free.py, manybody.py,
spectral.py).  The preparation and evolution are Eq. (A.6):

    |Psi(0)> = c^dag_j |Omega_N> / sqrt(1 - <n_j>),
    |Psi(t)> = exp(-i (H_{N+1} - E_0^N) t) |Psi(0)>,

with |Omega_N> chosen by Step 2's state-selection rule ("parity+").  The guide's
primary calculation uses j = 0; j = L/2 is the centre-injection control.

Observables (Eqs. A.7-A.9):

    Delta n_j(t) = <Psi(t)|n_j|Psi(t)> - <Omega|n_j|Omega>
    Q_edge(t)    = sum_{j < l_e} Delta n_j(t),   l_e = 2, 4
    G(t)         = <Psi(0)|Psi(t)> = sum_m z_jm exp(-i omega_m t)
    P_return(t)  = |G(t)|^2

The E_0^N shift is kept in the propagator, so G(t) is compared with Step 3's
spectral sum as a COMPLEX amplitude, phase included -- a stronger test than
comparing |G|^2 alone.

Exact identities used as checks (all hold for any V):

    sum_j Delta n_j(t) = 1                         (one particle added)
    Delta n_j(0) = 1 - <n_j>                       for the injection site j
    Delta n_k(0) = -<dn_j dn_k>_c / (1 - <n_j>)    for k != j
    norm and <H_{N+1}> are conserved

Krylov propagation
------------------
Each output interval dt is covered by one short-iterative-Lanczos step whose
Krylov dimension grows until the a posteriori error estimate

    err = ||psi|| * beta_m * |e_m^T exp(-i T_m dt) e_1|

falls below tol * dt.  The tolerance is therefore an error per unit time, so the
accumulated bound after time t is tol * t.  If the dimension cap is hit the step
is split in two, recursively.  Halving `tol` is the guide's item 1 check.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray
from scipy.linalg import eigh_tridiagonal

import manybody as mb
import spectral as spx
import ssh_free as sf

__all__ = [
    "Preparation", "prepare", "prepare_filtered", "Trajectory",
    "krylov_evolve", "exact_evolve", "spectral_amplitude", "compare",
    "finite_size_time", "return_velocity", "free_fermion_trajectory",
    "far_edge_response", "free_front_velocity", "subthreshold_lines",
    "round_trip",
    "TOL", "DT_OUT", "T_MAX", "ELL_E",
]

TOL = 1e-10          # Krylov error per unit time (production)
DT_OUT = 0.05        # output spacing, as the guide asks
T_MAX = 8.0          # initial time range
ELL_E = (2, 4)       # edge-region sizes for Q_edge
M_MAX = 40           # Krylov dimension cap before a step is split


# --------------------------------------------------------------------------
# 1. preparation
# --------------------------------------------------------------------------
@dataclass
class Preparation:
    L: int
    N: int
    delta: float
    V: float
    site: int
    psi0: NDArray[np.complex128] = field(repr=False)     # normalised, N+1 sector
    H: sp.csr_matrix = field(repr=False)                 # H_{N+1}
    sector: mb.Sector = field(repr=False)                # the N+1 sector
    E0_N: float = 0.0
    n_ref: NDArray[np.float64] = field(repr=False, default=None)  # <Omega|n_j|Omega>
    corr_ref: NDArray[np.float64] = field(repr=False, default=None)  # <dn_site dn_k>_c
    prep_weight: float = 1.0      # 1 - <n_site> unfiltered;  z_site(W) filtered
    filtered: bool = False
    window: dict | None = None
    gs_parity: float = np.nan
    gs_residual: float = np.nan
    gs_gap: float = np.nan

    @property
    def label(self) -> str:
        kind = "filtered" if self.filtered else "electron"
        where = "edge" if self.site == 0 else ("centre" if self.site == self.L // 2
                                                else f"j={self.site}")
        return f"L={self.L} d={self.delta:+.2f} V={self.V:g} {where} {kind}"


def prepare(L: int, delta: float, V: float, site: int = 0,
            convention: str = "parity+") -> Preparation:
    """|Psi(0)> = c^dag_site |Omega_N> / sqrt(1 - <n_site>)  (Eq. A.6)."""
    N = L // 2
    secN, secP = mb.Sector.make(L, N), mb.Sector.make(L, N + 1)
    HN = mb.build_hamiltonian(secN, delta, V)
    g = mb.select_ground_state(HN, secN, convention=convention)
    n_ref = mb.density_profile(g.psi, secN)

    # connected density-density correlations of Omega with the injection site,
    # needed for the exact t = 0 identity
    occ = secN.occupations().astype(float)
    p = g.psi**2
    nn = (occ * (occ[:, site] * p)[:, None]).sum(axis=0)
    corr = nn - n_ref[site] * n_ref

    phi = mb.cdag_matrix(L, N, site) @ g.psi
    w = float(phi @ phi)                    # = 1 - <n_site>, the sum rule
    HP = mb.build_hamiltonian(secP, delta, V)
    return Preparation(
        L=L, N=N, delta=delta, V=V, site=site,
        psi0=(phi / np.sqrt(w)).astype(complex), H=HP, sector=secP,
        E0_N=g.energy, n_ref=n_ref, corr_ref=corr, prep_weight=w,
        gs_parity=g.parity_bare, gs_residual=g.residual,
        gs_gap=g.gap_to_first_excited,
    )


def prepare_filtered(prep: Preparation, factor: float | None = None,
                     lanczos_m: int = 200) -> Preparation | None:
    """Item 5: project |Psi(0)> onto the edge spectral window and renormalise.

    The window is Step 3's: 0 <= omega < f Delta_+, with f the widest isolated
    factor in 0.6-0.9 unless one is given.  Returns None when the window is
    empty (strong end bonds): there is then nothing to filter onto.

    The preparation probability ||P_W Psi(0)||^2 equals z_site(W) of Step 3 --
    this is checked, not assumed.
    """
    L, d, V, j = prep.L, prep.delta, prep.V, prep.site
    spec = (spx.addition_spectrum_full(L, d, V, sites=(j,)) if L <= 12 else
            spx.addition_spectrum_lanczos(L, d, V, lanczos_m, sites=(j,)))
    f = spx.safe_factor(spec, j) if factor is None else factor
    n_w = spec.n_lines(j, f)
    hi = f * spec.delta_plus_ring
    if n_w == 0:
        return None

    dim = prep.H.shape[0]
    if dim <= 2000:
        vals, vecs = np.linalg.eigh(prep.H.toarray())
    else:
        vals, vecs = spla.eigsh(prep.H, k=min(n_w + 6, dim - 2), which="SA", tol=0.0)
        o = np.argsort(vals)
        vals, vecs = vals[o], vecs[:, o]
    om = vals - prep.E0_N
    inw = (om >= 0.0) & (om < hi)
    sub = vecs[:, inw]
    amp = sub.T @ prep.psi0                      # <m|Psi(0)>
    keep = np.abs(amp) ** 2 > spx.WEIGHT_FLOOR   # lines with usable weight only
    proj = sub[:, keep] @ amp[keep]
    pw = float(np.vdot(proj, proj).real)

    return Preparation(
        L=L, N=prep.N, delta=d, V=V, site=j, psi0=proj / np.sqrt(pw), H=prep.H,
        sector=prep.sector, E0_N=prep.E0_N, n_ref=prep.n_ref, corr_ref=prep.corr_ref,
        prep_weight=pw, filtered=True,
        window={"factor": f, "omega_hi": hi, "n_lines": int(keep.sum()),
                "n_lines_step3": n_w, "z_step3": spec.z(j, f),
                "delta_plus_ring": spec.delta_plus_ring},
        gs_parity=prep.gs_parity, gs_residual=prep.gs_residual, gs_gap=prep.gs_gap,
    )


# --------------------------------------------------------------------------
# 2. trajectories
# --------------------------------------------------------------------------
@dataclass
class Trajectory:
    prep: Preparation = field(repr=False)
    method: str
    tol: float
    times: NDArray[np.float64] = field(repr=False)
    amp: NDArray[np.complex128] = field(repr=False)       # G(t) = <Psi(0)|Psi(t)>
    dn: NDArray[np.float64] = field(repr=False)           # (T, L) Delta n_j(t)
    norm: NDArray[np.float64] = field(repr=False)
    energy: NDArray[np.float64] = field(repr=False)       # <H_{N+1}> - E_0^N
    krylov_dim_max: int = 0
    krylov_vectors_total: int = 0      # differs between two tolerances only if
    n_splits: int = 0                  # the tighter one did more work
    err_estimate: float = 0.0

    @property
    def P_return(self) -> NDArray[np.float64]:
        return np.abs(self.amp) ** 2

    def Q_left(self, ell: int) -> NDArray[np.float64]:
        """Q_edge of Eq. (A.8): excess charge on sites 0 .. ell-1."""
        return self.dn[:, :ell].sum(axis=1)

    def Q_right(self, ell: int) -> NDArray[np.float64]:
        return self.dn[:, -ell:].sum(axis=1)

    def Q_around(self, centre: int, ell: int) -> NDArray[np.float64]:
        """Excess charge in ell sites around `centre` -- the centre-injection
        analogue of Q_edge."""
        lo = max(0, centre - ell // 2)
        return self.dn[:, lo:lo + ell].sum(axis=1)

    def Q_local(self, ell: int) -> NDArray[np.float64]:
        """Charge remaining near wherever the particle was injected."""
        return self.Q_left(ell) if self.prep.site == 0 else \
            self.Q_around(self.prep.site, ell)

    def total_excess(self) -> NDArray[np.float64]:
        return self.dn.sum(axis=1)


def _observe(psi, prep: Preparation, occ):
    prob = np.abs(psi) ** 2
    dn = occ.T @ prob - prep.n_ref
    e = float(np.vdot(psi, prep.H @ psi).real) - prep.E0_N
    return complex(np.vdot(prep.psi0, psi)), dn, float(prob.sum()), e


def _krylov_step(H, psi, dt, tol, m_max, shift):
    """One short-iterative-Lanczos step exp(-i (H - shift) dt) psi."""
    beta0 = float(np.linalg.norm(psi))
    n = psi.size
    Vk = np.empty((m_max + 1, n), dtype=complex)
    Vk[0] = psi / beta0
    alpha: list[float] = []
    beta: list[float] = []
    err = np.inf
    for j in range(m_max):
        w = H @ Vk[j] - shift * Vk[j]
        a = float(np.vdot(Vk[j], w).real)
        alpha.append(a)
        w = w - a * Vk[j]
        if j > 0:
            w = w - beta[-1] * Vk[j - 1]
        Vj = Vk[: j + 1]
        for _ in range(2):                         # full reorthogonalisation
            w -= Vj.T @ (Vj.conj() @ w)
        b = float(np.linalg.norm(w))
        if j == 0:
            theta, S = np.array(alpha), np.ones((1, 1))
        else:
            theta, S = eigh_tridiagonal(np.array(alpha), np.array(beta))
        u = S @ (np.exp(-1j * theta * dt) * S[0, :])
        err = beta0 * b * abs(u[-1])
        if b < 1e-13 or err <= tol * dt:
            return beta0 * (Vj.T @ u), err, j + 1, True
        beta.append(b)
        Vk[j + 1] = w / b
    return None, err, m_max, False


def krylov_evolve(prep: Preparation, t_max: float = T_MAX, dt: float = DT_OUT,
                  tol: float = TOL, m_max: int = M_MAX) -> Trajectory:
    """Item 1: Krylov time evolution with output spacing dt."""
    times = np.round(np.arange(0.0, t_max + 0.5 * dt, dt), 12)
    occ = prep.sector.occupations().astype(float)
    T, L = times.size, prep.L
    amp = np.empty(T, complex)
    dn = np.empty((T, L))
    nrm = np.empty(T)
    en = np.empty(T)
    stats = {"m": 0, "m_total": 0, "splits": 0, "err": 0.0}

    def advance(psi, h):
        out, err, m, ok = _krylov_step(prep.H, psi, h, tol, m_max, prep.E0_N)
        if ok:
            stats["m"] = max(stats["m"], m)
            stats["m_total"] += m
            stats["err"] += err
            return out
        stats["splits"] += 1
        return advance(advance(psi, 0.5 * h), 0.5 * h)

    psi = prep.psi0.copy()
    amp[0], dn[0], nrm[0], en[0] = _observe(psi, prep, occ)
    for k in range(1, T):
        psi = advance(psi, times[k] - times[k - 1])
        amp[k], dn[k], nrm[k], en[k] = _observe(psi, prep, occ)

    return Trajectory(prep=prep, method="krylov", tol=tol, times=times, amp=amp,
                      dn=dn, norm=nrm, energy=en, krylov_dim_max=stats["m"],
                      krylov_vectors_total=stats["m_total"],
                      n_splits=stats["splits"], err_estimate=stats["err"])


def exact_evolve(prep: Preparation, t_max: float = T_MAX,
                 dt: float = DT_OUT) -> Trajectory:
    """Full-diagonalisation reference (the guide's L = 8 comparison)."""
    times = np.round(np.arange(0.0, t_max + 0.5 * dt, dt), 12)
    vals, vecs = np.linalg.eigh(prep.H.toarray())
    c = vecs.T @ prep.psi0
    phase = np.exp(-1j * np.outer(vals - prep.E0_N, times))   # (dim, T)
    occ = prep.sector.occupations().astype(float)
    T, L = times.size, prep.L
    amp, dn = np.empty(T, complex), np.empty((T, L))
    nrm, en = np.empty(T), np.empty(T)
    for k in range(T):
        psi = vecs @ (c * phase[:, k])
        amp[k], dn[k], nrm[k], en[k] = _observe(psi, prep, occ)
    return Trajectory(prep=prep, method="exact", tol=0.0, times=times, amp=amp,
                      dn=dn, norm=nrm, energy=en)


def spectral_amplitude(L: int, delta: float, V: float, site: int,
                       times: NDArray[np.float64], lanczos_m: int | None = None):
    """G(t) = sum_m z_jm exp(-i omega_m t), built from Step 3's addition spectrum.

    Full spectrum for L <= 12; for larger L the Lanczos spectral measure, whose
    G(t) is the single-Krylov-space approximation of the same amplitude.
    """
    spec = (spx.addition_spectrum_full(L, delta, V, sites=(site,))
            if lanczos_m is None else
            spx.addition_spectrum_lanczos(L, delta, V, lanczos_m, sites=(site,)))
    om, w = spec.omega_of(site), spec.w_of(site)
    ok = np.isfinite(om)
    z = w[ok] / (1.0 - spec.n_expect[site])
    return (z[None, :] * np.exp(-1j * np.outer(times, om[ok]))).sum(axis=1), spec


# --------------------------------------------------------------------------
# 3. comparisons and the finite-size window
# --------------------------------------------------------------------------
def compare(a: Trajectory, b: Trajectory) -> dict:
    """Largest differences between two trajectories on a common time grid."""
    n = min(a.times.size, b.times.size)
    return {
        "amp": float(np.max(np.abs(a.amp[:n] - b.amp[:n]))),
        "P_return": float(np.max(np.abs(a.P_return[:n] - b.P_return[:n]))),
        "dn": float(np.max(np.abs(a.dn[:n] - b.dn[:n]))),
        "energy": float(np.max(np.abs(a.energy[:n] - b.energy[:n]))),
    }


def finite_size_time(small: Trajectory, big: Trajectory, eps: float = 0.01,
                     ells: tuple[int, ...] = ELL_E) -> dict:
    """First time at which `small` visibly departs from the larger chain `big`.

    Compared: P_return and the local charge Q (Q_edge for edge injection, the
    centred window for centre injection), for each l_e.  Before this time the
    smaller chain behaves like the larger one -- no finite-size return has reached
    the observation region.  Returns t_max of the common grid with
    `reached = False` if no departure is seen.
    """
    n = min(small.times.size, big.times.size)
    t = small.times[:n]
    series = {"P_return": (small.P_return[:n], big.P_return[:n])}
    for ell in ells:
        series[f"Q{ell}"] = (small.Q_local(ell)[:n], big.Q_local(ell)[:n])
    first, offset = {}, {}
    for key, (x, y) in series.items():
        # Two different finite-size effects are kept apart.
        #  * Initial state: the chains already differ at t = 0 by a static amount
        #    s -- the edge-edge hybridisation of |Omega>, of order exp(-L/xi) --
        #    and that difference relaxes over a short transient of comparable
        #    size.  It is reported as `static_offset`, not mistaken for a return.
        #  * Return: a change in the dynamics larger than eps + 2 s, which the
        #    initial-state difference alone cannot produce.
        offset[key] = float(abs(x[0] - y[0]))
        d = np.abs((x - x[0]) - (y - y[0]))
        idx = np.flatnonzero(d > eps + 2.0 * offset[key])
        first[key] = float(t[idx[0]]) if idx.size else np.inf
    t_fs = min(first.values())
    s_max = max(offset.values())
    return {"t_fs": float(t_fs) if np.isfinite(t_fs) else float(t[-1]),
            "reached": bool(np.isfinite(t_fs)), "by": first, "eps": eps,
            "static_offset": offset, "static_offset_max": s_max,
            # below eps the small chain's initial state is converged at the
            # visibility level, so size-independence can be claimed at all
            "initial_state_ok": bool(s_max < eps),
            "t_common_max": float(t[-1])}


def round_trip(L: int, ell: int = 4) -> int:
    """Sites travelled from the injection site to the far end and back to the
    observation region [0, ell)."""
    return 2 * (L - 1) - (ell - 1)


def return_velocity(L: int, t_fs: float, ell: int = 4) -> float:
    """Speed implied by a certified return: round-trip distance / return time.

    Uses the same visibility criterion as the window itself.  Single-site
    arrival times are not used: they are spoiled here by the sublattice
    structure, and by precursors that enter the observation region LINEARLY
    through interference with the large static edge amplitude, so a visible
    return can precede the light cone of the main front.  The number is
    therefore an effective (upper-side) speed, not a group velocity.
    """
    return float(round_trip(L, ell) / t_fs)


def free_fermion_trajectory(L: int, delta: float, site: int,
                            times: NDArray[np.float64]):
    """V = 0 dynamics from Step 1's one-body matrix alone -- no Fock space.

    With C the filled-orbital correlation matrix, c^dag_site |Omega> creates a
    particle in chi(0) = (1 - C) e_site / sqrt(1 - C_ss), because the occupied
    part of e_site is Pauli-blocked.  Then chi(t) = exp(-i h t) chi(0),
    G(t) = <chi(0)|chi(t)> and Delta n_j(t) = |chi_j(t)|^2.  This is an
    independent implementation of the whole Step 4 calculation at V = 0.
    """
    h = sf.hopping_matrix(L, delta)
    eps, U = np.linalg.eigh(h)
    occ = U[:, : L // 2]
    C = occ @ occ.T
    e = np.zeros(L)
    e[site] = 1.0
    chi0 = (e - C @ e) / np.sqrt(1.0 - C[site, site])
    c = U.T @ chi0
    chi_t = U @ (c[:, None] * np.exp(-1j * np.outer(eps, times)))   # (L, T)
    G = chi0 @ chi_t
    return G, (np.abs(chi_t) ** 2).T


def subthreshold_lines(prep: Preparation, delta_plus: float, k: int = 12) -> list[dict]:
    """Every N+1 line below the ring threshold Delta_+ that the electron reaches.

    For each: omega, z (its share of the injected electron), and the excess
    charge it carries in the first 2 and 4 sites and in the last 4.  The sum of
    z * q over these lines is the edge charge that cannot leave -- the
    diagonal-ensemble share of the bound lines.  Comparing it with the measured
    late Q_edge separates bound charge from charge still in transit.
    """
    dim = prep.H.shape[0]
    if dim <= 2000:
        vals, vecs = np.linalg.eigh(prep.H.toarray())
    else:
        vals, vecs = spla.eigsh(prep.H, k=min(k, dim - 2), which="SA", tol=0.0)
        o = np.argsort(vals)
        vals, vecs = vals[o], vecs[:, o]
    occ = prep.sector.occupations().astype(float)
    out = []
    for m in range(vals.size):
        om = vals[m] - prep.E0_N
        if om >= delta_plus:
            break
        z = float(abs(vecs[:, m] @ prep.psi0) ** 2)
        if z < 1e-8:
            continue
        dn = (occ * (vecs[:, m] ** 2)[:, None]).sum(0) - prep.n_ref
        out.append({"omega": float(om), "z": z, "q2": float(dn[:2].sum()),
                    "q4": float(dn[:4].sum()), "q_far4": float(dn[-4:].sum()),
                    "q_left_half": float(dn[: prep.L // 2].sum())})
    return out


def far_edge_response(L: int, delta: float, V: float) -> dict:
    """Delta n at the far end at t = 0, for the selected state and for the two
    symmetry-broken alternatives built from the N-sector doublet.

    Delta n_{L-1}(0) = -<dn_0 dn_{L-1}>_c / (1 - <n_0>): a ground-state
    correlation, not transport.  In the parity-selected |Omega> the two end modes
    share one fermion, so this is O(1) and L-independent; in a left- or
    right-loaded combination of the doublet it is not.  The doublet is split by
    exp(-L/xi), so which of these is 'the' ground state is exactly Step 2's
    state-selection question.
    """
    N = L // 2
    sec = mb.Sector.make(L, N)
    H = mb.build_hamiltonian(sec, delta, V)
    g = mb.select_ground_state(H, sec)
    _, vecs = mb.lowest_states(H, k=2)
    occ = sec.occupations().astype(float)
    imb = occ[:, : L // 2].sum(1) - occ[:, L // 2:].sum(1)
    sub = vecs[:, :2]
    D = sub.T @ (imb[:, None] * sub)
    w, U = np.linalg.eigh(0.5 * (D + D.T))
    right_loaded, left_loaded = sub @ U[:, 0], sub @ U[:, 1]

    def resp(psi):
        p = psi**2
        n = (occ * p[:, None]).sum(0)
        nn = float((occ[:, 0] * occ[:, -1] * p).sum())
        conn = nn - n[0] * n[-1]
        return {"dn_far": -conn / (1.0 - n[0]), "n0": float(n[0]),
                "conn_0_far": conn}
    return {"selected": resp(g.psi), "left_loaded": resp(left_loaded),
            "right_loaded": resp(right_loaded),
            "doublet_splitting": float(g.gap_to_first_excited)}


def free_front_velocity(delta: float) -> float:
    """V = 0 maximal group velocity of the upper band, in sites per unit time.

    dE_+/dk = -J_w J_s sin k / E_+(k) per cell; its maximum is min(J_w, J_s),
    and a cell is two sites.
    """
    Jw, Js = sf.weak_strong(delta)
    k = np.linspace(0, np.pi, 200001)
    E = np.sqrt(Jw**2 + Js**2 + 2 * Jw * Js * np.cos(k))
    return float(2.0 * np.max(Jw * Js * np.sin(k) / E))
