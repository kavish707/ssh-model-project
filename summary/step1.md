# Step 1 summary: the free SSH edge problem ($V = 0$)

**Result.** With no interactions, a dimerised chain with weak end bonds has one exponentially localised state at each end, and the numerics confirm the analytic form of that state.

- The localisation length is $\xi = 2/\ln(J_s/J_w)$. At $\delta = -0.30$ this gives $\xi = 3.230813126$ sites, the value quoted in the project guide.
- The two end states split by $\Delta E \simeq 2J_s(1-r^2)\,e^{-L/\xi}$. This formula was derived here, and both its slope and prefactor were confirmed.
- Reversing the sign of $\delta$ removes the end states but leaves the bulk bands unchanged.
- All 37 checks in the notebook pass. One limit: the guide's supplied benchmark script was not available, so it has not been run (see [Section 7](#7-limits-and-what-was-not-tested)).

The code is in Step 1 of [`project_A.ipynb`](../project_A.ipynb). The write-up for submission is in [`report.tex`](../report.tex).

## Contents

1. [System and task](#1-system-and-task)
2. [Reduction to a matrix problem](#2-reduction-to-a-matrix-problem)
3. [Analytic derivations](#3-analytic-derivations)
4. [Code](#4-code)
5. [Results](#5-results)
6. [Verification](#6-verification)
7. [Limits and what was not tested](#7-limits-and-what-was-not-tested)
8. [Conclusions and what carries forward](#8-conclusions-and-what-carries-forward)
9. [How to reproduce](#9-how-to-reproduce)

---

## 1. System and task

The project studies spinless fermions on an open chain of $L$ sites (with $L$ even), labelled $j = 0, \dots, L-1$, at half filling $N = L/2$:

$$
H = -\sum_{j=0}^{L-2} J_j \left(c_j^\dagger c_{j+1} + \text{h.c.}\right) + V \sum_{j=0}^{L-2} \left(n_j - \tfrac12\right)\left(n_{j+1} - \tfrac12\right), \qquad J_j = 1 + (-1)^j \delta .
$$

Bond $j$ joins sites $j$ and $j+1$, and the bonds alternate in strength (the SSH chain). Units are $\hbar = 1$ with the mean hopping set to 1.

**Sign convention.** Since $J_0 = 1 + \delta$, a negative $\delta$ means the chain starts on a weak bond. Throughout,

$$
J_w = 1 - |\delta|, \qquad J_s = 1 + |\delta|, \qquad r = \frac{J_w}{J_s} < 1, \qquad M = \frac{L}{2}.
$$

Because $L$ is even, there is an odd number of bonds and $J_{L-2} = J_0$, so both ends have the same kind of bond.

**Step 1 sets $V = 0$.** The guide asks for four things, and all four were done:

| Guide item | Where it is done |
|---|---|
| Build the hopping matrix for $L = 40, 80, 120$ and $\delta = \pm 0.15, \pm 0.30, \pm 0.50$ | Sections [2](#2-reduction-to-a-matrix-problem), [4](#4-code) |
| Plot the energies and the two states nearest zero; separate the hybridised ends | Sections [3.4](#34-separating-the-two-ends), [5.1](#51-spectrum-and-edge-states) |
| Fit the first 6–10 nonzero even-site weights; recover $\xi = 3.230813126$ at $\delta = -0.30$ | Section [5.2](#52-localisation-length) |
| Derive the ring bands, gap $4\lvert\delta\rvert$ and addition cost $2\lvert\delta\rvert$ | Sections [3.5](#35-bulk-bands-on-a-ring), [5.4](#54-ring) |

The guide's three conditions for moving on are also met: the fitted length agrees with the formula, the end splitting decreases with $L$, and reversing $\delta$ removes the low-energy end pair.

---

## 2. Reduction to a matrix problem

At $V = 0$ the Hamiltonian is quadratic:

$$
H_0 = \sum_{ij} c_i^\dagger h_{ij} c_j, \qquad h_{j,j+1} = h_{j+1,j} = -J_j, \quad \text{all other entries zero.}
$$

Let $h\phi_n = \varepsilon_n\phi_n$ with real orthonormal eigenvectors. The operators $d_n = \sum_j \phi_n(j)\,c_j$ then satisfy $\{d_n, d_m^\dagger\} = \delta_{nm}$ and $H_0 = \sum_n \varepsilon_n d_n^\dagger d_n$.

So the many-body problem reduces to an $L \times L$ matrix: the ground state fills the $N$ lowest orbitals. This is why chains as long as $L = 120$ take milliseconds here, whereas the interacting problem in Step 2 is limited to about $L = 18$.

---

## 3. Analytic derivations

### 3.1 Sublattice (chiral) symmetry

Hopping only connects even sites to odd sites. With $\Gamma = \mathrm{diag}\big((-1)^j\big)$, this means $\Gamma h \Gamma = -h$. If $h\phi = \varepsilon\phi$, then $h(\Gamma\phi) = -\varepsilon\,\Gamma\phi$. Two consequences follow:

- The spectrum is symmetric under $\varepsilon \to -\varepsilon$, and partners $\phi$ and $\Gamma\phi$ have equal probability on every site.
- A state at exactly $\varepsilon = 0$ can be chosen to live on one sublattice only.

### 3.2 The edge wavefunction

Take a semi-infinite chain $j \ge 0$ with $\delta < 0$, so $J_{2\ell} = J_w$ and $J_{2\ell+1} = J_s$. Look for $E = 0$ with $\psi_{2\ell+1} = 0$.

The rows of $h\psi = 0$ at even sites contain only odd-site amplitudes, so they hold automatically. The row at odd site $2\ell+1$ gives

$$
J_w\psi_{2\ell} + J_s\psi_{2\ell+2} = 0 \quad\Longrightarrow\quad \psi_{2\ell} = \psi_0(-r)^\ell .
$$

Normalising with $\sum_{\ell \ge 0} r^{2\ell} = 1/(1-r^2)$ gives

$$
\boxed{\psi_{2\ell} = \sqrt{1-r^2}\,(-r)^\ell, \qquad \psi_{2\ell+1} = 0.}
$$

Setting $j = 2\ell$ gives $|\psi_j| = \sqrt{1-r^2}\,r^{j/2} = \sqrt{1-r^2}\,e^{-j/\xi}$, with

$$
\boxed{\xi = \frac{2}{\ln(J_s/J_w)} = \frac{1}{\operatorname{artanh}|\delta|} \simeq \frac{1}{|\delta|}.}
$$

The factor of 2 appears because the amplitude drops by $r$ per two-site cell, not per site. The probability decays twice as fast as the amplitude, and this is the form every fit must use:

$$
|\psi_j|^2 = (1-r^2)\,e^{-2j/\xi} \qquad (j \text{ even}).
$$

| $\lvert\delta\rvert$ | $r$ | $1-r^2$ | $\xi$ (exact) | $1/\lvert\delta\rvert$ |
|---|---|---|---|---|
| 0.15 | 0.739130 | 0.453686 | 6.616363079 | 6.667 |
| 0.30 | 0.538462 | 0.710059 | 3.230813126 | 3.333 |
| 0.50 | 0.333333 | 0.888889 | 1.820478453 | 2.000 |

**Why a strong end bond has no such state.** For $\delta > 0$ the same recursion gives $\psi_{2\ell} \propto (-1/r)^\ell$, which grows and cannot be normalised. An odd-sublattice state fails too: the row at site 0 reads $-J_0\psi_1 = 0$, so $\psi_1 = 0$, and the recursion then sets every odd amplitude to zero. So a left end has a zero mode only if its first bond is weak.

**The right end.** Site $L-1$ is odd. By reflection, a $\delta < 0$ chain also has $\psi^R_{L-1-2\ell} = \sqrt{1-r^2}\,(-r)^\ell$, which lives on the odd sublattice. Call the left state $\psi^L$.

### 3.3 Splitting of the two end states

Truncate $\psi^L$ to $L$ sites. It still solves every row of $h\psi^L = 0$ except the last, where a neighbour is missing:

$$
(h\psi^L)_{L-1} = -J_w\,\psi^L_{L-2} .
$$

Now work in the two-dimensional space spanned by $\psi^L$ and $\psi^R$. These are orthogonal because they occupy different sublattices, and the diagonal elements of $h$ vanish because $h$ maps one sublattice to the other. The only remaining element is the tunnelling amplitude

$$
t = \langle\psi^R|h|\psi^L\rangle = -J_w\,\psi^R_{L-1}\,\psi^L_{L-2} = -J_w(1-r^2)(-r)^{M-1}.
$$

Diagonalising the $2\times 2$ matrix with zero diagonal and off-diagonal $t$ gives $E_\pm = \pm|t|$. Using $J_w = rJ_s$ and $r^M = e^{-L/\xi}$,

$$
\boxed{\Delta E = E_+ - E_- \simeq 2J_s(1-r^2)\,r^{L/2} = 2J_s(1-r^2)\,e^{-L/\xi}.}
$$

This makes two testable predictions: $\ln\Delta E$ against $L$ has slope $-1/\xi$, and the prefactor is $2J_s(1-r^2)$. The sign of $t$ is $(-1)^M$, which matters in [Section 3.6](#36-the-half-filled-state-and-the-values-at-v--0).

### 3.4 Separating the two ends

The two exact eigenstates nearest zero, $\phi_\pm$, are spread over both ends, and each one is a sum or difference of $\psi^L$ and $\psi^R$. Worse, once $\Delta E$ reaches round-off, an eigensolver may return any rotation of them.

So we fix the rotation by diagonalising the position operator $X = \mathrm{diag}(j)$ inside $\mathrm{span}\{\phi_+, \phi_-\}$. The eigenvector with the smaller $\langle X\rangle$ is taken as $\psi^L$.

Diagonalising $\Gamma$ in the same subspace gives the same state, for this reason:

- $\Gamma$ maps $\phi_+$ to $\phi_-$ (up to a sign), so the subspace is $\Gamma$-invariant.
- It is therefore spanned by the sublattice-polarised states $\psi^{L,R} \propto \phi_+ \pm \Gamma\phi_+$.
- $X$ is already diagonal in that basis, because $\langle\psi^L|X|\psi^R\rangle = 0$ for states on disjoint sublattices.

The code computes both routes as a cross-check.

### 3.5 Bulk bands on a ring

Group the sites into $M$ two-site cells, $a_n = c_{2n}$ and $b_n = c_{2n+1}$, with $a_M \equiv a_0$. The bond inside a cell is $J_0 = 1+\delta$ and the bond between cells is $J_1 = 1-\delta$:

$$
H_0 = -\sum_{n=0}^{M-1}\Big[J_0\big(a_n^\dagger b_n + \text{h.c.}\big) + J_1\big(b_n^\dagger a_{n+1} + \text{h.c.}\big)\Big].
$$

Substitute $a_n = M^{-1/2}\sum_k e^{ikn}a_k$ (and likewise for $b_n$), with $k = 2\pi m/M$ and $m = 0, \dots, M-1$. Using $\sum_n e^{i(k-k')n} = M\delta_{kk'}$, the bond between cells becomes $\sum_n b_n^\dagger a_{n+1} = \sum_k e^{ik}\,b_k^\dagger a_k$, so

$$
H_0 = \sum_k \begin{pmatrix} a_k^\dagger & b_k^\dagger \end{pmatrix}
\begin{pmatrix} 0 & -f(k)^* \\ -f(k) & 0 \end{pmatrix}
\begin{pmatrix} a_k \\ b_k \end{pmatrix},
\qquad f(k) = J_0 + J_1 e^{ik}.
$$

The eigenvalues are $\pm|f(k)|$, and $|f|^2 = J_0^2 + J_1^2 + 2J_0J_1\cos k$. This is symmetric under $J_0 \leftrightarrow J_1$, so for either sign of $\delta$

$$
\boxed{E_\pm(k) = \pm\sqrt{J_w^2 + J_s^2 + 2J_wJ_s\cos k}.}
$$

**Gap.** $|E|$ is smallest at $k = \pi$, where $|E| = J_s - J_w = 2|\delta|$. It is largest at $k = 0$, where $|E| = 2$. The band gap is therefore $E_+(\pi) - E_-(\pi) = 4|\delta|$.

**Addition cost.** At half filling, $N = M$ particles exactly fill the lower band. An added particle must enter the upper band, and the cheapest state is the band bottom, so $E_0^{N+1} - E_0^N = E_+(\pi) = 2|\delta|$. Removing a particle costs the same, which adds up to the full gap.

**When this holds.** Both results need $k = \pi$ to be an allowed momentum, which requires $M$ even ($L$ divisible by 4). For odd $M$ the allowed momenta miss $k = \pi$, so the gap is $2\min_k E_+(k)$ over the allowed $k$, which is larger than $4|\delta|$.

**What the $\pm\delta$ symmetry means.** The two dimerisations have identical bulk spectra, so the bulk energies alone cannot tell you whether an edge state exists.

### 3.6 The half-filled state and the values at $V = 0$

At half filling exactly one of the end pair is occupied. Choosing which one is a real issue. Once $\Delta E$ is below round-off, "fill the lowest $N$ orbitals" picks an arbitrary state, and edge observables then change by order one. So we fix the choice with an explicit rule: fill the bonding orbital

$$
\phi_\mp = \frac{\psi^L \mp s\,\psi^R}{\sqrt 2}, \qquad s = \operatorname{sgn} t = (-1)^{L/2},
$$

where the overall signs are set by $\psi^L_0 > 0$ and $\psi^R_{L-1} > 0$. Here $\phi_-$ is filled and $\phi_+$ is empty. The sign $s$ comes from the analytic result in [Section 3.3](#33-splitting-of-the-two-end-states), not from a numerical $t$, which is unreliable at the same points where the choice matters.

Two exact results follow, and later steps must reproduce them at $V = 0$.

**Uniform density.** Completeness gives $\sum_n \phi_n(j)^2 = 1$. Chiral partners carry equal weight on each site, and the rule above fills exactly one member of each pair (note that $\Gamma\phi_- = \phi_+$). Therefore

$$
\langle n_j\rangle = \tfrac12 \quad \text{on every site.}
$$

**Edge addition weight.** Adding a particle at site 0 can reach the empty orbital $\phi_+$, at energy $\omega = \Delta E/2 \approx 0$, while all bulk lines lie at $\omega \ge 2|\delta|$. Since $\psi^R_0 = 0$ exactly, $|\phi_+(0)|^2 = |\psi^L_0|^2/2$, and so

$$
z_0 = \frac{|\phi_+(0)|^2}{1 - \langle n_0\rangle} = |\psi^L_0|^2 \simeq 1 - r^2 .
$$

The addition weight also decays as $W_j \propto |\psi^L_j|^2 \propto e^{-2j/\xi}$, so the spectral localisation length at $V = 0$ is $\xi_W = \xi$. These are the starting points of the Step 3 curves $z_0(V)$ and $\xi_W(V)$.

---

## 4. Code

Everything is in the Step 1 section of [`project_A.ipynb`](../project_A.ipynb). The notebook imports `numpy`, `pandas`, `matplotlib`, `mpmath` and `IPython`.

### Parameters (the cell at the top of Step 1)

| Variable | Value | Meaning |
|---|---|---|
| `L_LIST` | `[40, 80, 120]` | Open-chain lengths |
| `DELTA_LIST` | `[-0.50, -0.30, -0.15, 0.15, 0.30, 0.50]` | Dimerisation; negative means weak end bonds |
| `N_FIT_LIST` | `[6, 8, 10]` | Number of even sites used in the fit |
| `L_SCAN` | `12, 16, ..., 120` | Lengths for the splitting plot |
| `L_EX`, `DELTA_EX` | `40`, `-0.30` | Example used in the spectrum, end-pair and ring plots |
| `L_PROFILE` | `80` | Length used in the edge-profile plot |
| `XI_GUIDE` | `3.230813126` | The guide's value of $\xi$ at $\delta = -0.30$ |
| `L_ODD_M` | `[18, 42, 62, 82]` | Open chains with odd $M$, which test $s = -1$ |
| `L_RING_ODD` | `[10, 14, 18]` | Rings with odd $M$ (the Step 6 sizes) |

### Functions (section 1.1 of the notebook)

| Function | What it does |
|---|---|
| `bonds(L, delta, periodic)` | Returns $J_j = 1 + (-1)^j\delta$ for $L-1$ bonds, or $L$ on a ring |
| `hopping_matrix(L, delta, periodic)` | Builds the matrix $h$ |
| `weak_strong(delta)`, `ratio(delta)` | Return $(J_w, J_s)$ and $r$ |
| `xi_exact(delta)` | Returns $2/\ln(J_s/J_w)$ |
| `zero_mode(L, delta)` | Returns the analytic $\psi^L$, truncated to $L$ sites |
| `edge_pair(L, delta)` | Diagonalises $h$ and builds $\psi^L$, $\psi^R$, $\phi_\mp$ |
| `fit_xi(prob, n_fit)` | Fits $\xi$ from even-site probabilities |
| `splitting_asym(L, delta)` | Returns the formula $2J_s(1-r^2)r^{L/2}$ |
| `splitting_sturm(L, delta)` | Returns the exact splitting in 80-digit arithmetic |
| `bloch_bands(k, delta)` | Returns $E_\pm(k)$ |

**The hopping matrix.** On a ring, the `% L` adds the closing bond $J_{L-1} = 1-\delta$:

```python
def hopping_matrix(L, delta, periodic=False):
    h = np.zeros((L, L))
    for j, Jj in enumerate(bonds(L, delta, periodic)):
        k = (j + 1) % L
        h[j, k] = h[k, j] = -Jj
    return h
```

**Separating the ends.** This is the core of `edge_pair`, implementing [Section 3.4](#34-separating-the-two-ends) and the rule in [Section 3.6](#36-the-half-filled-state-and-the-values-at-v--0). Here `P` holds the two eigenvectors nearest zero as columns. `P.T @ X @ P` is the $2\times2$ position matrix within that subspace, and multiplying by its eigenvectors gives the left and right states:

```python
E, U = np.linalg.eigh(h)
pair = np.sort(np.argsort(np.abs(E))[:2])
P = U[:, pair]

_, vx = np.linalg.eigh(P.T @ np.diag(np.arange(L, dtype=float)) @ P)
left, right = P @ vx[:, 0], P @ vx[:, 1]              # smaller <X> first
left, right = left * np.sign(left[0]), right * np.sign(right[-1])

s = (-1) ** (L // 2)
bonding     = (left - s * right) / np.sqrt(2)         # phi_-, filled
antibonding = (left + s * right) / np.sqrt(2)         # phi_+, empty
```

**The fit.** `prob` is $|\psi_j|^2$. The sites used are $j = 0, 2, \dots, 2(n_\text{fit}-1)$, and the slope of $\ln|\psi_j|^2$ is $-2/\xi$:

```python
def fit_xi(prob, n_fit):
    j = 2 * np.arange(n_fit)
    slope = np.polyfit(j, np.log(prob[j]), 1)[0]
    return -2.0 / slope
```

**The exact splitting.** Ordinary floating point cannot resolve a splitting of $10^{-29}$, so an independent high-precision route is needed. It works like this:

1. Eliminating the tridiagonal matrix $h - x$ row by row gives the pivots $q_0 = -x$ and $q_j = -x - J_{j-1}^2/q_{j-1}$.
2. By Sylvester's law of inertia, the number of negative pivots equals the number of eigenvalues below $x$.
3. The spectrum has $M$ negative levels. So this count is $M$ for $0 < x < \varepsilon_0$ and $M+1$ just above $\varepsilon_0$, where $\varepsilon_0 = \Delta E/2$ is the smallest positive level.
4. Bisecting on $x$ finds that step.

The search uses the geometric mean because $\varepsilon_0$ can be many orders of magnitude below 1. The starting bracket is $[10^{-70}, |\delta|]$, and an `assert` confirms it contains exactly one step. The search stops when the bracket ratio is within $10^{-30}$ of 1. The `tiny` value guards against an exactly zero pivot, which really does occur at $x = 0.5$, $\delta = -0.5$:

```python
def splitting_sturm(L, delta, dps=80):
    with mpmath.workdps(dps):
        J2 = [(1 + (-1) ** j * mpmath.mpf(delta)) ** 2 for j in range(L - 1)]
        tiny = mpmath.mpf(10) ** (-2 * dps)

        def count_below(x):
            q, n = -x, int(x > 0)
            for b in J2:
                q = -x - b / (q if q != 0 else tiny)
                n += int(q < 0)
            return n

        lo, hi = mpmath.mpf(10) ** (-dps + 10), mpmath.mpf(abs(delta))
        assert count_below(lo) == L // 2 and count_below(hi) == L // 2 + 1
        while hi / lo - 1 > mpmath.mpf(10) ** -30:
            mid = mpmath.sqrt(lo * hi)
            lo, hi = (mid, hi) if count_below(mid) == L // 2 else (lo, mid)
        return float(2 * lo)
```

### Notebook sections

| Section | What it computes | Output |
|---|---|---|
| 1.2 | Bond convention, symmetry of $h$, the $\pm E$ spectrum, $\xi$, zero-mode residual | 7 checks, table of $r$ and $\xi$ |
| 1.3 | Number of in-gap states for every $(L, \delta)$ | 3 checks, `step1_spectrum.pdf` |
| 1.4 | $\psi^L$ and fitted $\xi$ for $\delta < 0$ | 4 checks, `step1_edge_states.pdf` |
| 1.5 | Splitting by three routes, and its slope | 7 checks, `step1_splitting.pdf` |
| 1.6 | Ring spectrum against $E_\pm(k)$, gap, addition cost | 4 checks, `step1_ring_bands.pdf` |
| 1.7 | The half-filling rule, $\langle n_0\rangle$, $z_0$ | 4 checks |
| 1.8 | The same rule and the splitting formula at odd $M$; odd-$M$ rings | 8 checks |

Two eigensolvers appear. `edge_pair` uses `np.linalg.eigh`, because it needs eigenvectors, while the spectrum, splitting and ring cells use `np.linalg.eigvalsh`. They use different algorithms and do not reach the same accuracy on the tiny splitting: at $L = 80$, $\delta = -0.5$, `eigvalsh` returns round-off while the eigenvectors from `eigh` still give the correct $|t|$. The "resolved" column in the table below refers to `eigvalsh`.

---

## 5. Results

### 5.1 Spectrum and edge states

([`step1_spectrum.pdf`](../figures/step1_spectrum.pdf), [`step1_edge_states.pdf`](../figures/step1_edge_states.pdf))

- For every $\delta < 0$ there are exactly two levels with $|E| < |\delta|$.
- For every $\delta > 0$ there are none. The lowest $|E|$ is $1.001$ to $1.092$ times $2|\delta|$, which places it at the bulk band edge.
- The two eigenstates nearest zero are each spread over both ends, as expected (panel (a) of the edge-state figure).
- The separated left state matches $(1-r^2)e^{-2j/\xi}$ over about fourteen decades with no free parameters (panel (b)). Its total weight on odd sites is below $3\times10^{-29}$.

### 5.2 Localisation length

Fractional errors are $(\xi_\text{fit} - \xi)/\xi$. $\Delta E$ (exact) comes from the 80-digit calculation. "Resolved" means `eigvalsh` agrees with the exact value to $10^{-6}$.

| $\delta$ | $L$ | $\xi_\text{fit}$ ($n_\text{fit} = 10$) | Error, $n_\text{fit} = 6$ | Error, $n_\text{fit} = 8$ | Error, $n_\text{fit} = 10$ | $\Delta E$ (exact) | Resolved? |
|---|---|---|---|---|---|---|---|
| $-0.15$ | 40 | 6.613912628 | $-6.0\times10^{-5}$ | $-1.5\times10^{-4}$ | $-3.7\times10^{-4}$ | $2.4713\times10^{-3}$ | yes |
| $-0.15$ | 80 | 6.616363065 | $-3.4\times10^{-10}$ | $-8.3\times10^{-10}$ | $-2.1\times10^{-9}$ | $5.8523\times10^{-6}$ | yes |
| $-0.15$ | 120 | 6.616363079 | $-1.5\times10^{-15}$ | $-4.8\times10^{-15}$ | $-1.2\times10^{-14}$ | $1.3860\times10^{-8}$ | yes |
| $-0.30$ | 40 | 3.230812682 | $-2.3\times10^{-9}$ | $-1.7\times10^{-8}$ | $-1.4\times10^{-7}$ | $7.7511\times10^{-6}$ | yes |
| $-0.30$ | 80 | 3.230813126 | $0$ | $1.6\times10^{-15}$ | $2.6\times10^{-15}$ | $3.2543\times10^{-11}$ | yes |
| $-0.30$ | 120 | 3.230813126 | $1.4\times10^{-15}$ | $1.4\times10^{-16}$ | $-4.1\times10^{-16}$ | $1.3663\times10^{-16}$ | no |
| $-0.50$ | 40 | 1.820478453 | $-2.4\times10^{-16}$ | $-3.2\times10^{-14}$ | $-1.7\times10^{-12}$ | $7.6479\times10^{-10}$ | yes |
| $-0.50$ | 80 | 1.820478453 | $3.7\times10^{-16}$ | $2.4\times10^{-16}$ | $2.4\times10^{-16}$ | $2.1934\times10^{-19}$ | no |
| $-0.50$ | 120 | 1.820478453 | $0$ | $-2.4\times10^{-16}$ | $0$ | $6.2906\times10^{-29}$ | no |

The fit reproduces the formula to round-off on long chains. The largest error, $3.7\times10^{-4}$, is at $\delta = -0.15$, $L = 40$, where $L/\xi$ is only 6.

**Why that point deviates.** The error comes from the opposite edge, not from the fit. Two things support this. The error grows as the fit reaches further into the chain, and it falls rapidly with $L$ at fixed $\delta$.

It can also be estimated. The exact finite-chain $\psi^L$ is a combination of $\phi_\pm$, so $h\psi^L = \pm\varepsilon_0\psi^R$ rather than zero. The recursion therefore gains a source term:

$$
J_w\psi^L_{2\ell} + J_s\psi^L_{2\ell+2} = \mp\varepsilon_0\,\psi^R_{2\ell+1}.
$$

Relative to $J_s\psi^L_{2\ell+2}$, this source has size $(1-r^2)\,r^{L-2\ell-2}$. At $L = 40$, $|\delta| = 0.15$ that is about $1\times10^{-4}$ at $\ell = 5$ (the end of the $n_\text{fit} = 6$ window) and $1\times10^{-3}$ at $\ell = 9$ (the end of the $n_\text{fit} = 10$ window). The measured errors are $6\times10^{-5}$ and $3.7\times10^{-4}$. This is an order-of-magnitude estimate, and it agrees with the table at that level.

### 5.3 Splitting

([`step1_splitting.pdf`](../figures/step1_splitting.pdf))

The splitting falls monotonically with $L$. Fitting $\ln\Delta E$ against $L$ on the exact values with $L \ge 40$ gives:

| $\delta$ | Fitted slope | $-1/\xi$ | Relative difference |
|---|---|---|---|
| $-0.15$ | $-0.151141$ | $-0.151140$ | $1.5\times10^{-6}$ |
| $-0.30$ | $-0.309520$ | $-0.309520$ | $2.8\times10^{-12}$ |
| $-0.50$ | $-0.549306$ | $-0.549306$ | $2.2\times10^{-16}$ |

The full formula, prefactor included, is within $4.8\times10^{-5}$ of the exact value for all $L \ge 40$. The worst case is $\delta = -0.15$, $L = 40$. It reaches round-off for long chains and worsens on short ones: at $L = 18$ the error is $1.6\times10^{-2}$ for $\delta = -0.15$ and $8.8\times10^{-5}$ for $\delta = -0.30$.

Direct diagonalisation stops being reliable near $10^{-15}$. `eigvalsh` fails at three of the nine grid points, returning round-off of about $10^{-16}$ instead of the true value. The figure shows this floor directly. This is why the exact route was needed.

### 5.4 Ring

([`step1_ring_bands.pdf`](../figures/step1_ring_bands.pdf))

Over the whole grid (all $L$ divisible by 4), the sorted ring eigenvalues match the sorted values of $E_\pm(k)$ at the allowed $k$:

| Quantity | Largest deviation |
|---|---|
| Spectrum against $E_\pm(k)$ | $4.9\times10^{-15}$ |
| Gap against $4\lvert\delta\rvert$ | $1.6\times10^{-15}$ |
| Addition cost against $2\lvert\delta\rvert$ | $7.8\times10^{-16}$ |
| Spectrum at $+\delta$ against $-\delta$ | $3.8\times10^{-15}$ |

In the figure, ring eigenvalues are given momenta by matching sorted lists. That assignment is only for plotting; the check itself compares sorted eigenvalues.

**Odd-$M$ rings.** On the rings $L = 10, 14, 18$ the spectrum also matches $E_\pm(k)$ to $2.2\times10^{-15}$, and the gap equals $2\min_k E_+(k)$ to the same accuracy. The gap can be far above $4|\delta|$ on these small rings, and the excess is largest at small $|\delta|$:

| $\lvert\delta\rvert$ | $4\lvert\delta\rvert$ | Gap, $L = 10$ | Gap, $L = 14$ | Gap, $L = 18$ |
|---|---|---|---|---|
| 0.15 | 0.6 | 1.361 | 1.065 | 0.912 |
| 0.30 | 1.2 | 1.682 | 1.470 | 1.371 |
| 0.50 | 2.0 | 2.268 | 2.143 | 2.089 |

Does this spoil the Step 3 edge window $0 \le \omega < 0.8\Delta_+$, if $\Delta_+$ is taken from an odd-$M$ ring? At $V = 0$ it does not. A separate one-off comparison for $L = 8, 10, \dots, 18$ and the three negative $\delta$ found $0.8\Delta_+$ below the lowest open-chain bulk level of the same length in every case, because short open chains also have raised bulk levels. The closest case is $\delta = -0.50$, $L = 18$, with $0.8\Delta_+ = 0.835$ against a lowest bulk level of $1.049$. The same holds for the guide's largest factor, $0.9\Delta_+$.

This comparison uses same-length chains at $V = 0$ only, and it is not one of the notebook's checks. It says nothing about $V > 0$.

### 5.5 The half-filled state

The problem is real, not hypothetical. At $L = 120$, $\delta = -0.5$ the eigensolver returns the pure left and right states, with their tiny tails set exactly to zero. Filling the lowest $L/2$ orbitals then fills the whole left state, which gives $\langle n_0\rangle = 0.944$ instead of $\tfrac12$. That equals the prediction $r^2/2 + (1-r^2)$ for this situation. The edge addition weight also comes out as zero, which is wrong.

The numerically computed $t$ at that point is $-3.1\times10^{-29}$. Its magnitude is close to the true value, but its sign is wrong, which is why the rule uses the analytic sign.

With the rule from [Section 3.6](#36-the-half-filled-state-and-the-values-at-v--0), on all nine grid points:

- $\langle n_0\rangle = \tfrac12$ to $1.0\times10^{-15}$;
- $z_0 = |\psi^L_0|^2$ to round-off;
- $z_0 = 1 - r^2$ to within $4.0\times10^{-5}$, with the worst case again at $\delta = -0.15$, $L = 40$.

For example, $z_0 = 0.710059$ at $\delta = -0.30$.

**Odd $M$.** All the grid lengths have even $M$, where $s = +1$, so section 1.8 of the notebook tests $s = -1$ at $L = 18, 42, 62, 82$ for the three negative $\delta$. The numerical sign of $t$ is $-1$ at all 12 points. That includes the 2 points where `eigvalsh` does not resolve the splitting ($\delta = -0.5$ at $L = 62$ and $82$), but the check itself only counts the 10 resolved points. At those 10 points the bonding orbital also lies below the antibonding one. At all 12 points, $\langle n_0\rangle = \tfrac12$ and $z_0 = |\psi^L_0|^2$ to $1.0\times10^{-15}$. The splitting formula holds to $2.8\times10^{-5}$ for $L \ge 42$.

---

## 6. Verification

Each derivation is paired with a separate numerical calculation:

| Derived by hand | Checked against |
|---|---|
| Edge wavefunction and $\xi$ | Fit to eigenvectors of $h$; direct substitution into $h\psi$ |
| Splitting formula | 80-digit exact eigenvalue |
| $E_\pm(k)$, gap, addition cost | Direct diagonalisation of the ring matrix |
| $\langle n_j\rangle = \tfrac12$, $z_0 = \lvert\psi^L_0\rvert^2$ | Direct sum over filled orbitals |

The notebook runs 37 checks, printed as `[pass]` or `[FAIL]`. All pass. "Exact" means the check is a yes/no or integer comparison.

| Notebook section | Check | Largest deviation | Tolerance |
|---|---|---|---|
| 1.2 | $J_0 = J_w$ and $J_{L-2} = J_0$, for three $\delta$ (3 checks) | exact | exact |
| 1.2 | $h$ symmetric | exact | exact |
| 1.2 | Spectrum symmetric under $E \to -E$ | $2.4\times10^{-15}$ | $10^{-12}$ |
| 1.2 | $\xi(-0.30)$ against the guide's value | $4.4\times10^{-11}$ | $10^{-9}$ |
| 1.2 | Truncated zero mode: $h\psi = 0$ except the last row | $5.6\times10^{-17}$ | $10^{-12}$ |
| 1.3 | $\delta < 0$: exactly two in-gap states | exact | exact |
| 1.3 | $\delta > 0$: no in-gap states | exact | exact |
| 1.3 | $\delta > 0$: lowest $\lvert E\rvert \ge 2\lvert\delta\rvert$ | exact | $10^{-12}$ |
| 1.4 | Fitted $\xi$ against the formula, whole grid | $3.7\times10^{-4}$ | $10^{-3}$ |
| 1.4 | Fitted $\xi$ at $\delta = -0.30$ | $1.4\times10^{-7}$ | $10^{-3}$ |
| 1.4 | Left state: total weight on odd sites | $2.4\times10^{-29}$ | $10^{-12}$ |
| 1.4 | $X$ and $\Gamma$ routes give the same left state | $5.8\times10^{-15}$ | $10^{-12}$ |
| 1.5 | Splitting decreases with $L$, for three $\delta$ (3 checks) | exact | exact |
| 1.5 | Slope equals $-1/\xi$, for three $\delta$ (3 checks) | $1.5\times10^{-6}$ | $10^{-3}$ |
| 1.5 | Splitting formula against exact, $L \ge 40$ | $4.8\times10^{-5}$ | $10^{-4}$ |
| 1.6 | Ring spectrum against $E_\pm(k)$ | $4.9\times10^{-15}$ | $10^{-12}$ |
| 1.6 | Gap equals $4\lvert\delta\rvert$ | $1.6\times10^{-15}$ | $10^{-12}$ |
| 1.6 | Addition cost equals $2\lvert\delta\rvert$ | $7.8\times10^{-16}$ | $10^{-12}$ |
| 1.6 | Spectra identical for $\pm\delta$ | $3.8\times10^{-15}$ | $10^{-12}$ |
| 1.7 | $\operatorname{sgn} t = (-1)^{L/2}$ where resolved | no mismatch | exact |
| 1.7 | $\langle n_0\rangle = \tfrac12$ | $1.0\times10^{-15}$ | $10^{-12}$ |
| 1.7 | $z_0 = \lvert\psi^L_0\rvert^2$ | round-off | $10^{-12}$ |
| 1.7 | $z_0$ against $1-r^2$ | $4.0\times10^{-5}$ | $10^{-3}$ |
| 1.8 | Odd $M$: $\operatorname{sgn} t = -1$ where resolved | no mismatch | exact |
| 1.8 | Odd $M$: bonding below antibonding where resolved | exact | exact |
| 1.8 | Odd $M$: $\langle n_0\rangle = \tfrac12$ | $8.3\times10^{-16}$ | $10^{-12}$ |
| 1.8 | Odd $M$: $z_0 = \lvert\psi^L_0\rvert^2$ | $1.0\times10^{-15}$ | $10^{-12}$ |
| 1.8 | Odd $M$: splitting formula against exact, $L \ge 42$ | $2.8\times10^{-5}$ | $10^{-4}$ |
| 1.8 | Odd-$M$ rings: spectrum against $E_\pm(k)$ | $2.2\times10^{-15}$ | $10^{-12}$ |
| 1.8 | Odd-$M$ rings: gap equals $2\min_k E_+(k)$ | $2.2\times10^{-15}$ | $10^{-12}$ |
| 1.8 | Odd-$M$ rings: gap $> 4\lvert\delta\rvert$ | exact | exact |

The $10^{-3}$ tolerances are deliberately loose: they cover the known finite-size error at $L/\xi = 6$, not numerical noise.

**A failure that exposed the half-filling problem.** The checks on $\langle n_0\rangle$ and $z_0$ originally failed at $L = 120$, $\delta = -0.5$, with a deviation of 0.44 in $\langle n_0\rangle$. That failure is what revealed the issue described in [Section 5.5](#55-the-half-filled-state). The fix was to define the filled state explicitly. No tolerance was loosened to make the check pass.

---

## 7. Limits and what was not tested

- **Guide benchmark not run.** The guide mentions `checks/benchmark_checks.py` and `checks/benchmark_results.json`, but they were not in the available files. The only comparison with the guide's numbers is the quoted $\xi = 3.230813126$.
- **The results $4|\delta|$ and $2|\delta|$ apply only to rings with $L$ divisible by 4.** The Step 6 rings have odd $M$ and larger gaps ([Section 5.4](#54-ring)).
- **The sign rule rests on the analytic argument where the splitting is unresolved.** Numerically it can only be confirmed where $t$ is above round-off. It is a chosen convention at those points, not a measured result.
- **The $V = 0$ values for $\langle n_j\rangle$, $z_0$ and $\xi_W$ were not tested in a many-body code.** $\langle n_0\rangle$ and $z_0$ were computed from orbitals at site 0 only, and $\xi_W$ is a prediction. Testing them is a task for Steps 2 and 3.
- **The size estimate in [Section 5.2](#52-localisation-length) is order-of-magnitude only.**
- **The derived splitting formula is leading order.** It is accurate on the long chains tested here but less so on the short chains used in Step 2. Those short chains also have a much larger splitting: at $\delta = -0.30$, $L = 18$, $\Delta E = 7.0\times10^{-3}$. The two end states there are not close to degenerate on the scale of a converged eigensolver.

---

## 8. Conclusions and what carries forward

**What Step 1 establishes:**

1. The edge state is the analytic zero mode. On long chains the fit reproduces $\xi$ to round-off. At small $L/\xi$ the only visible correction comes from the opposite edge.
2. The two ends communicate only through an exponentially small tunnelling, $\Delta E \simeq 2J_s(1-r^2)e^{-L/\xi}$.
3. The edge state depends on how the bulk is terminated, not on the bulk energies: $+\delta$ and $-\delta$ have identical ring spectra, yet only the weak-end chain has an end pair.

**What later steps need from this:**

1. **An explicit ground-state rule.** Once the end pair is nearly degenerate, the half-filled state must be selected deliberately, or edge observables change by order one. The guide's Step 2 warns about this same problem for the interacting chain.
2. **Values to reproduce at $V = 0$:** $\langle n_j\rangle = \tfrac12$, $z_0 = 1 - r^2$ (up to finite size) and $\xi_W = \xi$. If the many-body code misses these, it is wrong.
3. **The bulk threshold:** at $V = 0$, $\Delta_+ = 2|\delta|$ for rings with $L$ divisible by 4, and half the larger gap in [Section 5.4](#54-ring) for odd $M$. This sets the edge window in Step 3, so a threshold should be quoted with its ring size.
4. **Fit windows must stay away from the far edge.** The error scales as $r^{L-2\ell}$. This matters more in Steps 2 and 3, where chains are much shorter.

---

## 9. How to reproduce

1. Use an environment with `numpy`, `pandas`, `matplotlib`, `mpmath` and `ipykernel`. These results came from Python 3.12.3, NumPy 2.5.2, pandas 3.0.5, Matplotlib 3.11.1 and mpmath 1.3.0.
2. Open [`project_A.ipynb`](../project_A.ipynb) and run all cells from the top.
3. The last Step 1 cell should print `Step 1: 37/37 checks passed`. Figures are written to [`figures/`](../figures/).
