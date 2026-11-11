# Step 3 summary: localisation and coherence

Read this after [`step2.md`](step2.md), or at least its Sections 3 and 6. Each section answers one question. Sources: the guide is `QFT_student_research_guides.pdf`, Project A, p. 8; the code is Step 3 of [`project_A.ipynb`](../project_A.ipynb); the write-up is Section 4 of [`report.tex`](../report.tex).

## Short version

Steps 1 and 2 built the edge state and the many-body machine. Step 3 asks a physical question with it: **if you add one electron at the end of the chain, how much of it becomes one single eigenstate, the edge excitation, and how far into the chain does that excitation reach?** The first is $z_0$, the second $\xi_W$.

**Done and verified, on $L = 8, 12, 16, 18$, $\delta = \pm0.30$, $V = 0, 0.5, 1, 2, 3$:**

- The addition spectrum at $L = 8, 12$ (full spectrum) and $L = 16, 18$ (Lanczos spectral measure), with 14 checks (Section 5). They include a hand-derived $L = 2$ result at $V \neq 0$, the exact $V = 0$ solution at the same $L$ (Step 1, notebook section 1.9), a separate calculation in the full Fock space at $L = 4, 6$, and the Lanczos measure tested against the full spectrum at $L = 12$.
- Weak end bonds ($\delta = -0.30$): one sharp edge line, the lowest state of the $N+1$ sector, far below the bulk threshold at every $L$ and $V$. Strong end bonds ($\delta = +0.30$): nothing in the window, at every $L$, $V$ and window factor.
- Its weight falls with repulsion and is **converged in size**: $z_{\rm edge} = 0.710, 0.702, 0.681, 0.597, 0.456$ at $V = 0, 0.5, 1, 2, 3$ ($L = 18$); $L = 16$ and $18$ agree to $0.002$ or better.
- **The window at $V = 3$ is not isolated at $L = 16, 18$.** A second line ($z = 0.104$ at $L = 16$, $0.084$ at $L = 18$) enters it for $f = 0.9$ ($L = 16$) and $f \ge 0.8$ ($L = 18$), so the window weight $z_0$ there depends on $f$ (0.457 or 0.562; 0.456 or 0.540). Everywhere else $z_0$ is the same for every $f$ from 0.6 to 0.9.
- **Localisation changes less than coherence.** At $L \ge 16$, $\xi_W$ decreases with $V$, by 3–4% on a fixed fit range $j = 0, 2, 4$ and by 9–12% on all even sites $j < L/2$; the profile stops being a single exponential (rms deviation of $\ln W_j$ up to about 0.1). $L = 12$ shows no decrease and $L = 8$ an increase: the small sizes are not converged.
- Odd sites carry the far edge; adding at the centre gives almost no window weight (at most 0.017 at $L = 12$, 0.005 at $L = 16, 18$).

**Open:** what the second line is; sizes beyond $L = 18$.

| Symbol | Meaning |
|---|---|
| $\Omega_N$ | the selected half-filled ground state (Step 2, Section 6) |
| $w_{jm}$ | probability that adding a particle at site $j$ gives exactly eigenstate $m$ of the $N+1$ sector |
| $\omega_m$ | $E_m^{N+1} - E_0^N$, the cost of that outcome |
| $z_{jm}$ | $w_{jm}/(1-\langle n_j\rangle) = 2w_{jm}$ |
| $\Delta_+$ | ring addition threshold at the same $(L,\delta,V)$ (Step 2) |
| $W$ | edge window, $0 \le \omega < f\Delta_+$, $f = 0.8$ (0.6–0.9 tested) |
| $z_0$ | total $z_{0m}$ inside $W$ |
| $z_{\rm edge}$ | $z_{0m}$ of the lowest line alone; equal to $z_0$ where the window holds one line |
| $W_j$ | total $w_{jm}$ inside $W$ at site $j$: the profile of the window line |
| $\xi_W$ | decay length from a least-squares fit $\ln W_j = a - 2j/\xi_W$ |
| $\phi_+$ | the antibonding end orbital of Step 1, `edge_pair(L, delta)["antibonding"]` |

---

## 1. What did Step 3 have to produce?

| Guide requirement (p. 8) | Where | Status |
|---|---|---|
| 1. Full addition spectrum at $L \le 12$; lines first; Lorentzian widths 0.02, 0.01, 0.005 for presentation only | notebook 3.1, 3.3, 3.9 | done at $L = 8, 12$; broadened figure `step3_broadened` |
| 2. $z_{jm}$ for isolated lines; summed weight for clusters | 3.1 (`merge_lines`), 3.2 (A3, A9) | done |
| 3. Edge window from the ring, factor 0.8, varied 0.6 to 0.9; say so if no isolated window | 3.2 (A5, A6), 3.4, 3.6 | done; not isolated at $V = 3$, $L = 16$ ($f = 0.9$) and $L = 18$ ($f \ge 0.8$), reported in Table `step3_lanczos.csv` and report Table 'window at $V=3$' |
| 4. $W_j$; fit $e^{-2j/\xi_W}$, even and odd sites separately; fit range and omitted weight; compare edge and centre spectra | 3.5, 3.7, 3.8 | done at every $L$: fit range, rms, fixed-range fit, odd-site fit, weight on fitted sites, $\sum_j W_j$; edge and centre spectra `step3_edge_centre` |
| 5. $L = 16, 18$ with the Lanczos spectral measure, 100, 200, 400 vectors with reorthogonalisation; no full sum rule from a few lines | 3.6 (A10–A14) | done; converged at 100 vectors |
| Before moving on: one figure with $\xi_W(V)$ and $z_0(V)$ side by side; missing weight, an unstable window or a poor fit must stay visible | 3.7, `figures/step3_z0_xi.pdf` | done for $L = 8$–$18$; the poor fit at $L = 18$, $V = 3$ is left in |

## 2. What is being measured, and why?

Take the half-filled chain in its ground state and add one electron exactly on site $j$. Then ask which eigenstate it is in. Quantum mechanically there is no single answer: the result $c_j^\dagger|\Omega_N\rangle$ is a superposition of the eigenstates $|m, N+1\rangle$, and

$$w_{jm} = \big|\langle m,N+1|\,c_j^\dagger\,|\Omega_N\rangle\big|^2$$

is the probability of finding it in state $m$. Plotting $w_{jm}$ at $\omega_m$ gives the addition spectrum (guide Eq. A.5).

- **One tall line** means the electron operator creates mostly one eigen-excitation: the edge excitation is **coherent**.
- **Many small lines** mean the added electron spreads over many states, even if it stays near the edge.

That is why the guide separates coherence ($z_0$, weight at one site) from localisation ($W_j$ and $\xi_W$, where along the chain the window line lives): they answer different questions, and here they behave differently.

**At $V = 0$** (Step 1, notebook sections 1.7 and 1.9): the added electron can only enter empty orbitals, with weight $\phi_n(j)^2$. The only empty orbital in the window is the antibonding end orbital $\phi_+$, so $W_j = \phi_+(j)^2$. On even sites this is $|\psi_L(j)|^2/2$, so $z_0 = |\psi^L_0|^2$ exactly at every $L$. On odd sites it is $|\psi_R(j)|^2/2$: the right edge state's tail, which grows towards the far end.

The benchmarks at $V = 0$ are the exact values **at the same $L$** (summary/step1.md §5.6), not the $L\to\infty$ limits:

| $L$ | $z_0$ | $\xi_W$ | fit sites |
|---|---|---|---|
| 8 | 0.73091 | 3.1652 | 0, 2 |
| 12 | 0.71289 | 3.2161 | 0, 2, 4 |
| 16 | 0.71040 | 3.2276 | 0, 2, 4, 6 |
| 18 | 0.71017 | 3.2285 | 0, 2, 4, 6, 8 |
| limit | $1 - r^2 = 0.71006$ | $\xi = 3.2308$ | |

## 3. Why do the weights add to ½?

Insert the complete set of $N+1$ eigenstates, then use $c_jc_j^\dagger = 1 - n_j$:

$$\sum_m w_{jm} = \sum_m \langle\Omega_N|c_j|m\rangle\langle m|c_j^\dagger|\Omega_N\rangle = \langle\Omega_N|c_jc_j^\dagger|\Omega_N\rangle = 1 - \langle n_j\rangle .$$

$\langle n_j\rangle = \tfrac12$ at every $V$ by the particle-hole symmetry $C$ (Step 2, Section 6), so the total is $\tfrac12$. Physically, you can only add an electron where the site is empty, which is half the time. Interactions can move weight between lines but cannot change the total.

## 4. How does the code compute it?

`addition_spectrum(L, delta, V, j)` in notebook 3.1 follows the definition:

1. `half_filled_state` gives $\Omega_N$; its energy is $E_0^N = \langle\Omega|H|\Omega\rangle$.
2. `c_dag(L, N, j) @ psi` gives $v = c_j^\dagger\Omega_N$, with the Jordan-Wigner sign string, in sector $N+1$.
3. `np.linalg.eigh` on the dense $N+1$ Hamiltonian gives every eigenpair $(E_m, u_m)$; at most $792\times792$ at $L = 12$.
4. `(U.T @ v) ** 2` gives every $w_{jm} = (u_m\cdot v)^2$ at once.

`edge_window` returns $f\Delta_+$ using Step 2's `ring_gap`. `merge_lines` adds up the weights of lines closer than `ETA_DEGENERATE` (Step 2), because only the total weight of a degenerate subspace is defined. `brute_force_addition` (check A9 only) builds every $c_i$ as a $2^L\times2^L$ Kronecker product and shares no code with Step 2.

**$L = 16, 18$** (notebook 3.6): the $N+1$ sector has 11 440 and 43 758 states. `lanczos_spectrum(H, v, k)` starts from $q_0 = v/|v|$, builds $k$ orthonormal vectors by applying $H$ and removing every earlier direction (twice), and diagonalises the resulting $k\times k$ tridiagonal matrix. Its eigenvalues are the lines; $|v|^2$ times the squared first components of its eigenvectors are their weights. `addition_spectra_lanczos` does this for a list of sites with the ground state and $H$ built once. These weights add to $|v|^2$ by construction, so the sum rule is not a test at these sizes.

**Fits** (notebook 3.7): `fit_line(W, sites)` fits a straight line to $(j, \ln W_j)$ and returns $\xi = -2/\text{slope}$ and the rms deviation of $\ln W_j$ from the line. It is used on the even sites $j < L/2$ (the `BENCH` sites), on the fixed range $j = 0, 2, 4$, on the odd sites $j < L/2$, and on the lowest line's own profile.

**Broadening** (notebook 3.9): `broaden` replaces each line by a Lorentzian of width $\eta$, for the figure only.

## 5. How do we know it is right?

Step 2 verified the matrices and $c_j^\dagger$. What is new is how they are put together and the Lanczos measure, so the checks target those. The worst deviations below are the ones in [`data/step3_checks.csv`](../data/step3_checks.csv).

| Id | What | Worst deviation | Independent of the new code because |
|---|---|---|---|
| A1 | sum rule $=\tfrac12$, every site, $L = 8, 12$, both $\delta$, every $V$ | $3.2\times10^{-13}$ | exact identity; fails if $\Omega$, $c^\dagger$ or the eigenbasis is wrong |
| A2 | $L = 2$: one line at $\omega = J_0 + V/2$, weight $\tfrac12$ | $8.9\times10^{-16}$ | derived by hand: $E_0 = -V/4 - J_0$, $E(\lvert 11\rangle) = V/4$; tests $V \neq 0$ |
| A3 | $V = 0$: lines at the empty orbital energies with weights $\phi_n(j)^2$, nothing elsewhere, every site | $6.1\times10^{-14}$ | uses only the $L\times L$ hopping matrix of Step 1 |
| A4 | $V = 0$: $z_0 = \lvert\psi^L_0\rvert^2$ at the same $L$ | $3.5\times10^{-14}$ | Step 1's separated left edge state |
| A5 | one line with weight at $j = 0$ in the window for $\delta = -0.30$, none for $+0.30$, every $V$ and $f$ ($L = 8, 12$) | exact | a result as well as a check: the window is isolated there |
| A6 | one eigenstate of sector $N+1$ in the window for $\delta = -0.30$, whatever its weight, none for $+0.30$ ($L = 8, 12$) | exact | so $W_j$ is the profile of one state |
| A7 | $V = 0$: $W_j = \phi_+(j)^2$ at every site, same $L$ | $6.4\times10^{-14}$ | Step 1's antibonding orbital |
| A8 | $V = 0$: $\xi_W$ = `BENCH[(L, delta)]["xi_W"]`, same fit sites (relative) | $1.9\times10^{-15}$ | Step 1's `free_benchmark` |
| A9 | $L = 4, 6$: lines and weights = full-Fock-space calculation, both $\delta$, every $V$ and site | $8.9\times10^{-15}$ | separate construction, no Step 2 code |
| A10 | $L = 12$: Lanczos with 100, 200, 400 vectors = full spectrum (window weight every $f$, lowest three lines and weights), both $\delta$, every $V$ | $1.2\times10^{-14}$ | the full spectrum of A1–A9 |
| A11 | $L = 16, 18$: window weight with 100 and 200 vectors = with 400, every point and $f$; 200-vector profile = 400 at $j = 0$ | $0$ | convergence in the number of vectors |
| A12 | $L = 16, 18$, $V = 0$: $z_0$ and the whole $W_j$ = Step 1 same-size values | $1.2\times10^{-12}$ | Step 1's `free_benchmark` |
| A13 | $L = 16, 18$: lowest and second line at $j = 0$ = Step 2's `add_gap` and `add_gap + add_level1_gap` | $6.3\times10^{-14}$ | Step 2's own eigenvalues (`data/step2_grid.csv`) |
| A14 | $L = 16, 18$, $\delta = +0.30$: window empty, every $V$ and $f$ | exact | the control |

A9 prints the $L = 4$ lines at $j = 0$, $\delta = -0.30$: weights $0.4201, 0.0799$ at $V = 0$ and $0.3723, 0.1194, 0.0067, 0.0016$ at $V = 2$.

## 6. What do the results show?

Figures: [`step3_spectra`](../figures/step3_spectra.pdf), [`step3_z0_xi`](../figures/step3_z0_xi.pdf), [`step3_profile`](../figures/step3_profile.pdf), [`step3_edge_centre`](../figures/step3_edge_centre.pdf), [`step3_broadened`](../figures/step3_broadened.pdf). Tables: [`step3_z0.csv`](../data/step3_z0.csv) ($L = 8, 12$), [`step3_lanczos.csv`](../data/step3_lanczos.csv) ($L = 16, 18$), [`step3_xi.csv`](../data/step3_xi.csv) (every size), [`step3_profile.csv`](../data/step3_profile.csv).

**Weight of the edge line**, $z_{\rm edge}$ at $j = 0$, $\delta = -0.30$ (equal to $z_0$ for every $f$ except where noted):

| $V$ | $L = 8$ | $L = 12$ | $L = 16$ | $L = 18$ |
|---|---|---|---|---|
| 0 | 0.731 | 0.713 | 0.710 | 0.710 |
| 0.5 | 0.706 | 0.700 | 0.702 | 0.702 |
| 1 | 0.673 | 0.675 | 0.680 | 0.681 |
| 2 | 0.589 | 0.588 | 0.595 | 0.597 |
| 3 | 0.503 | 0.464 | 0.457 ($z_0 = 0.562$ at $f = 0.9$) | 0.456 ($z_0 = 0.540$ at $f \ge 0.8$) |

**The window at $V = 3$:**

| $L$ | $\Delta_+$ | $0.8\Delta_+$ | $0.9\Delta_+$ | lowest line $\omega$, $z$ | second line $\omega$, $z$ |
|---|---|---|---|---|---|
| 12 | 1.566 | 1.253 | 1.410 | 0.162, 0.464 | 1.502, 0.163 |
| 16 | 1.481 | 1.185 | 1.333 | 0.064, 0.457 | 1.255, 0.104 |
| 18 | 1.556 | 1.245 | 1.400 | 0.041, 0.456 | 1.181, 0.084 |

**Decay length** ($\delta = -0.30$):

| $V$ | 0 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| $L = 8$, $j = 0, 2$ | 3.1652 | 3.4101 | 3.5844 | 3.7330 | 3.6211 |
| $L = 12$, $j = 0, 2, 4$ | 3.2161 | 3.2026 | 3.2042 | 3.2303 | 3.2201 |
| $L = 16$, $j = 0..6$ | 3.2276 | 3.0899 | 3.0086 | 2.9484 | 2.9529 |
| $L = 18$, $j = 0..8$ | 3.2285 | 3.0347 | 2.9176 | 2.8212 | 2.8296* |
| $L = 16$, $j = 0, 2, 4$ | 3.2296 | 3.1760 | 3.1443 | 3.1232 | 3.1224 |
| $L = 18$, $j = 0, 2, 4$ | 3.2305 | 3.1699 | 3.1334 | 3.1055 | 3.1057* |
| rms of $\ln W_j$, $L = 18$, $j < L/2$ | 0.001 | 0.048 | 0.079 | 0.108 | 0.297 (window) |

\* lowest line alone: at $L = 18$, $V = 3$ the window holds two lines and its profile gives $\xi_W = 4.122$ with rms 0.297, a poor fit, left in the figure.

Reading it:

- **One sharp edge line** for weak end bonds at every $L$ and $V$, the lowest $N+1$ state, at Step 2's addition gap. **None** for strong end bonds. The line is an edge effect.
- **Its weight falls with repulsion** from about 0.71 to 0.456, and **does not fragment with size**: $L = 16$ and $18$ agree to 0.002. The fall at $V = 3$ between $L = 8$ and $12$ levels off.
- **The window is not isolated at $V = 3$, $L \ge 16$.** The second line moves down with $L$ (1.502, 1.255, 1.181) while its weight at $j = 0$ falls (0.163, 0.104, 0.084), and it crosses $f\Delta_+$ for $f = 0.9$ at $L = 16$ and $f \ge 0.8$ at $L = 18$. We report $z_0$ for each $f$ and the lowest line separately. What the second line is has not been established; at $L = 12$, $V = 3$ its weight is largest at the ends, so we do not call it bulk.
- **Where the rest goes** ($L = 12$): at or above $\Delta_+$, except the second line at $V = 3$ (split 0.464 / 0.163 / 0.373).
- **Localisation.** At $V = 0$, $\xi_W$ equals the same-size value at every $L$. With repulsion, at $L = 16, 18$, $\xi_W$ decreases (3–4% on $j = 0, 2, 4$; 9–12% on $j < L/2$) and the rms deviation grows from 0.001 to about 0.1, so the profile is not a single exponential and $\xi_W$ depends on the fit range. $L = 12$ gives no decrease on the same range and $L = 8$ an increase: those sizes are not converged. Read at $L \ge 16$: the edge line becomes slightly more confined to the edge while its weight falls.
- **Odd sites and omitted weight.** Odd-site fits give negative lengths ($-2.3$ to $-4.6$): those sites carry the right edge. The fitted even sites hold 46–50% of $\sum_j W_j$; the rest is mostly the mirror image at the right end. $\sum_j W_j$ falls from 1 to 0.65–0.67 at $V = 3$.
- **Edge against centre.** Added at the centre, the window weight is at most 0.017 ($L = 12$) and 0.005 ($L = 16, 18$; 0.012 at $L = 18$, $V = 3$, where the window holds two lines); the centre spectrum lies above $\Delta_+$.

**A defensible statement** (in the guide's form): *on $L \le 18$ and up to $V = 3$, adding an electron at the weak-bond end creates one sharp edge excitation whose weight falls from 0.710 to 0.456 and is converged to 0.002 between $L = 16$ and $18$; its decay length changes much less (a decrease of 3–12% at $L \ge 16$, depending on the fit range). At $V = 3$, $L \ge 16$, a second line with weight 0.08–0.10 lies inside the guide's window for some $f$, so the window weight there depends on $f$.*

## 7. Limits

- **Sizes.** Up to $L = 18$, i.e. $L/\xi \le 5.6$. The edge weight is converged on this grid; $\xi_W$ is not converged at $L \le 12$ and depends on the fit range at all sizes.
- **Odd $M$ at $L = 18$.** $L = 8, 12, 16$ have even $M = L/2$, so $\Delta_+ = 2|\delta| = 0.600$ at $V = 0$; $L = 18$ has $M = 9$ and $\Delta_+ = 0.685$, so its window edge is higher for a reason unrelated to the edges.
- **Lanczos.** The weights add to $\frac12$ by construction, so the sum rule is not tested at $L = 16, 18$; the low lines are converged (A11) and agree with Step 2's eigenvalues (A13).
- **The second line** is not identified.
- **Broadening** is presentation only.

## 8. Files and how to reproduce

| File | Contents |
|---|---|
| [`project_A.ipynb`](../project_A.ipynb) | Step 3 sections 3.1 to 3.9; Steps 1 and 2 unchanged |
| [`report.tex`](../report.tex) | Section 4 |
| [`figures/step3_spectra.pdf`](../figures/step3_spectra.pdf) | addition spectra, $L = 12$, $j = 0$, $\delta = \pm0.30$, $V = 0, 1, 2, 3$ |
| [`figures/step3_z0_xi.pdf`](../figures/step3_z0_xi.pdf) | $z_0(V)$ and $\xi_W(V)$ side by side, $L = 8$–$18$ |
| [`figures/step3_profile.pdf`](../figures/step3_profile.pdf) | $z_j$ of the window line and of the next line along the chain, $L = 12$ |
| [`figures/step3_edge_centre.pdf`](../figures/step3_edge_centre.pdf) | spectra at $j = 0$ and $j = 6$, $L = 12$ |
| [`figures/step3_broadened.pdf`](../figures/step3_broadened.pdf) | Lorentzian widths 0.02, 0.01, 0.005, presentation only |
| [`data/step3_z0.csv`](../data/step3_z0.csv) | $L = 8, 12$: $z_0$ for each $f$, weight split at $\Delta_+$, edge and next line, $\Delta_+$ |
| [`data/step3_lanczos.csv`](../data/step3_lanczos.csv) | $L = 16, 18$: $z_0$ and number of lines in the window for each $f$, lowest two lines and weights, $\Delta_+$ |
| [`data/step3_xi.csv`](../data/step3_xi.csv) | every size: $z_0$, $z_{\rm edge}$, lines in the window, $\xi_W$ with rms, fixed-range and odd-site fits, lowest-line fit, weight on fitted sites, $\sum_j W_j$, centre weight |
| [`data/step3_profile.csv`](../data/step3_profile.csv) | $W_j$, $z_j$ of the window line and of the next line, every site, every size |
| [`data/step3_checks.csv`](../data/step3_checks.csv) | every Step 3 check with its worst deviation |

Run the whole notebook, or the Step 3 cells after Steps 1 and 2. Step 3 takes about six minutes (most of it the $L = 18$ profile) and prints `Step 3: 14/14 checks passed`.

## 9. Record of LLM assistance

The guide asks for a short record of substantive assistance and how it was verified.

- **Tool.** Claude (Anthropic), in claude.ai sessions on 2026-11-04, 2026-11-05, 2026-11-11 and 2026-11-12.
- **What it did.** Explained the observable, the sum rule and the window. Wrote `addition_spectrum`, which I typed in and ran myself after a line-by-line explanation. Wrote the check cells A1 to A14, `merge_lines`, `brute_force_addition`, `lanczos_spectrum`, `addition_spectra_lanczos`, `fit_line`, `broaden`, the profile, fit, table and figure cells, the notebook text, report Section 4 and this summary, built on my first runs. On 2026-11-05 it made the changes asked for in the group leader's review. On 2026-11-12 it added guide items 4 and 5 that were missing: the Lanczos measure at $L = 16, 18$, the odd-site and fixed-range fits with residuals and omitted weight, edge and centre spectra, and broadening.
- **How it was verified.** A1 to A4, A7, A8 and A12 compare against results that do not use the new code: the sum rule, a hand-derived $L = 2$ formula, and Step 1's exact solution at the same $L$. A9 compares with a separate construction that shares no code with Step 2. A10 tests the Lanczos measure against the full spectrum, A11 its convergence, and A13 its lines against Step 2's own eigenvalues. A9 and the Lanczos code were written by the same LLM as `addition_spectrum`, so they guard against errors in assembling the pieces, not against a shared misunderstanding of the definition. The numbers in this summary are the ones the notebook prints.
- **What I still need to do.** Be able to derive the sum rule and the $L = 2$ check by hand, explain each line of `addition_spectrum` and `lanczos_spectrum`, and decide how to interpret the second line at $V = 3$.
