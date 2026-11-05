# Step 3 summary: localisation and coherence (in progress)

Read this after [`step2.md`](step2.md), or at least its Sections 3 and 6. Each section answers one question. Sources: the guide is `QFT_student_research_guides.pdf`, Project A, p. 8; the code is Step 3 of [`project_A.ipynb`](../project_A.ipynb); the write-up is Section 4 of [`report.tex`](../report.tex).

## Short version

Steps 1 and 2 built the edge state and the many-body machine. Step 3 asks a physical question with it: **if you add one electron at the end of the chain, how much of it becomes one single eigenstate, the edge excitation, and how far into the chain does that excitation reach?** The first is $z_0$, the second $\xi_W$.

**Done and verified, for $L = 8, 12$:**

- The addition spectrum is built from Step 2's matrices and operators, and passes 9 checks (Section 5). They include a hand-derived $L = 2$ result at $V \neq 0$, the exact $V = 0$ solution at the same $L$ (Step 1, notebook section 1.9), and a separate calculation in the full Fock space at $L = 4, 6$.
- Weak end bonds ($\delta = -0.30$): exactly one eigenstate lies inside the edge window at every $V \le 3$. Strong end bonds ($\delta = +0.30$): none.
- That line's weight falls with repulsion: $z_0 = 0.713, 0.700, 0.675, 0.588, 0.464$ at $V = 0, 0.5, 1, 2, 3$ ($L = 12$).
- The rest of the weight lies at or above $\Delta_+$ at every point except $L = 12$, $V = 3$. There the split is 0.464 in the window, 0.163 in one line at $\omega = 1.502$ just below $\Delta_+ = 1.566$, and 0.373 above $\Delta_+$. What that line is has not been established (Section 6).
- $L = 8$ and $L = 12$ agree for $V \le 2$ (to $0.0063$ or better for $V > 0$) and differ at $V = 3$ ($0.503$ against $0.464$).
- $z_0$ does not depend on the window factor between 0.6 and 0.9.
- At $L = 12$, $\xi_W$ stays within 0.5% of its $V = 0$ value 3.216 up to $V = 3$. At $L = 8$ (a two-site fit) it rises to 3.73, so $\xi_W(V)$ is not converged in size.

**Not done yet:** $L = 16, 18$ with the Lanczos spectral measure (item 5), where the window must be re-tested (Section 7); the odd-site profile and the edge-versus-centre comparison (item 4); presentation broadening (item 1).

| Symbol | Meaning |
|---|---|
| $\Omega_N$ | the selected half-filled ground state (Step 2, Section 6) |
| $w_{jm}$ | probability that adding a particle at site $j$ gives exactly eigenstate $m$ of the $N+1$ sector |
| $\omega_m$ | $E_m^{N+1} - E_0^N$, the cost of that outcome |
| $z_{jm}$ | $w_{jm}/(1-\langle n_j\rangle) = 2w_{jm}$ |
| $\Delta_+$ | ring addition threshold at the same $(L,\delta,V)$ (Step 2) |
| $W$ | edge window, $0 \le \omega < f\Delta_+$, $f = 0.8$ |
| $z_0$ | total $z_{0m}$ inside $W$ |
| $W_j$ | total $w_{jm}$ inside $W$ at site $j$: the profile of the window line |
| $\xi_W$ | decay length from a fit $W_j \propto e^{-2j/\xi_W}$ on the even sites $j < L/2$ |
| $\phi_+$ | the antibonding end orbital of Step 1, `edge_pair(L, delta)["antibonding"]` |

---

## 1. What did Step 3 have to produce?

| Guide requirement (p. 8) | Where | Status |
|---|---|---|
| 1. Full addition spectrum at $L \le 12$; lines first; broadening 0.02, 0.01, 0.005 for presentation only | notebook 3.1, 3.3 | lines done at $L = 8, 12$; broadening not done |
| 2. $z_{jm}$ for isolated lines; summed weight for clusters | 3.1 (`merge_lines`), 3.2 (A3, A9) | done |
| 3. Edge window from the ring, factor 0.8, varied 0.6 to 0.9; say so if no isolated window | 3.2 (A5, A6), 3.4 | done at $L = 8, 12$: window isolated at every point |
| 4. $W_j$, fit $e^{-2j/\xi_W}$ on even and odd sites; compare edge and centre spectra | 3.2 (A7, A8), 3.5 | $W_j$ and the even-site fit done at $L = 8, 12$; odd sites and edge-versus-centre not done |
| 5. $L = 16, 18$ with Lanczos (100, 200, 400 vectors, reorthogonalised) | | not done |
| Before moving on: one figure with $\xi_W(V)$ and $z_0(V)$ side by side | 3.5, `figures/step3_z0_xi.pdf` | done at $L = 8, 12$ |

## 2. What is being measured, and why?

Take the half-filled chain in its ground state and add one electron exactly on site $j$. Then ask which eigenstate it is in. Quantum mechanically there is no single answer: the result $c_j^\dagger|\Omega_N\rangle$ is a superposition of the eigenstates $|m, N+1\rangle$, and

$$w_{jm} = \big|\langle m,N+1|\,c_j^\dagger\,|\Omega_N\rangle\big|^2$$

is the probability of finding it in state $m$. Plotting $w_{jm}$ at $\omega_m$ gives the addition spectrum (guide Eq. A.5).

- **One tall line** means the electron operator creates mostly one eigen-excitation: the edge excitation is **coherent**.
- **Many small lines** mean the added electron spreads over many states, even if it stays near the edge.

That is why the guide separates coherence ($z_0$, weight at one site) from localisation ($W_j$ and $\xi_W$, where along the chain the window line lives): they answer different questions.

**At $V = 0$** (Step 1, notebook sections 1.7 and 1.9): the added electron can only enter empty orbitals, with weight $\phi_n(j)^2$. The only empty orbital in the window is the antibonding end orbital $\phi_+$, so $W_j = \phi_+(j)^2$. On even sites this is $|\psi_L(j)|^2/2$, so $z_0 = |\psi^L_0|^2$ exactly at every $L$. On odd sites it is $|\psi_R(j)|^2/2$: the right edge state's tail, which grows towards the far end. That is why $\xi_W$ is fitted on even sites only.

The benchmarks at $V = 0$ are the exact values **at the same $L$** (summary/step1.md §5.6), not the $L\to\infty$ limits:

| $L$ | $z_0$ | $\xi_W$ | fit sites |
|---|---|---|---|
| 8 | 0.73091 | 3.1652 | 0, 2 |
| 12 | 0.71289 | 3.2161 | 0, 2, 4 |
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
4. `(U.T @ v) ** 2` gives every $w_{jm} = (u_m\cdot v)^2$ at once. The $u_m$ are orthonormal, so $u_m\cdot v$ is the component of $v$ along $u_m$. All amplitudes are real on the open chain.

`edge_window` returns $f\Delta_+$ using Step 2's `ring_gap`. `merge_lines` adds up the weights of lines closer than `ETA_DEGENERATE` (Step 2), because only the total weight of a degenerate subspace is defined.

The profile is $W_j$ = the weight inside the window at each site, from the spectra at every $j$. $\xi_W$ comes from Step 1's `fit_xi` on the even sites $j < L/2$, the same call `free_benchmark` makes.

`brute_force_addition` is used only by check A9. It builds $c_i = Z\otimes\cdots\otimes Z\otimes a\otimes 1\otimes\cdots\otimes 1$ as $2^L\times2^L$ matrices ($a$ removes a particle, $Z = \mathrm{diag}(1,-1)$), writes $H$ from them and picks the sectors by particle number. It shares no code with Step 2.

## 5. How do we know it is right?

Step 2 verified the matrices and $c_j^\dagger$. What is new is how they are put together, so the checks target that. The worst deviations below are the ones in [`data/step3_checks.csv`](../data/step3_checks.csv).

| Id | What | Worst deviation | Independent of the new code because |
|---|---|---|---|
| A1 | sum rule $=\tfrac12$, every site, $L = 8, 12$, both $\delta$, every $V$ | $1.4\times10^{-13}$ | exact identity; fails if $\Omega$, $c^\dagger$ or the eigenbasis is wrong |
| A2 | $L = 2$: one line at $\omega = J_0 + V/2$, weight $\tfrac12$ | $8.9\times10^{-16}$ | derived by hand: $E_0 = -V/4 - J_0$, $E(\lvert 11\rangle) = V/4$; tests $V \neq 0$ |
| A3 | $V = 0$: lines at the empty orbital energies with weights $\phi_n(j)^2$, nothing elsewhere, every site | $4.9\times10^{-14}$ | uses only the $L\times L$ hopping matrix of Step 1 |
| A4 | $V = 0$: $z_0 = \lvert\psi^L_0\rvert^2$ at the same $L$ | $2.8\times10^{-14}$ | Step 1's separated left edge state |
| A5 | one line with weight at $j = 0$ in the window for $\delta = -0.30$, none for $+0.30$, every $V$ and $f$ | exact | a result as well as a check: the window is isolated |
| A6 | one eigenstate of sector $N+1$ in the window for $\delta = -0.30$, whatever its weight, none for $+0.30$, every $V$ and $f$ | exact | so $W_j$ is the profile of one state |
| A7 | $V = 0$: $W_j = \phi_+(j)^2$ at every site, same $L$ | $5.2\times10^{-14}$ | Step 1's antibonding orbital, $L\times L$ matrix only |
| A8 | $V = 0$: $\xi_W$ = `BENCH[(L, delta)]["xi_W"]`, same fit sites (relative) | $8.9\times10^{-16}$ | Step 1's `free_benchmark` |
| A9 | $L = 4, 6$: lines and weights = full-Fock-space calculation, both $\delta$, every $V$ and site | $7.1\times10^{-15}$ | separate construction, no Step 2 code; ground state unique there (smallest gap 0.303) |

A9 prints the $L = 4$ lines at $j = 0$, $\delta = -0.30$: weights $0.4201, 0.0799$ at $V = 0$ and $0.3723, 0.1194, 0.0067, 0.0016$ at $V = 2$.

## 6. What do the results show?

Look at [`figures/step3_spectra.pdf`](../figures/step3_spectra.pdf), [`figures/step3_z0_xi.pdf`](../figures/step3_z0_xi.pdf) and [`figures/step3_profile.pdf`](../figures/step3_profile.pdf). The full tables are [`data/step3_z0.csv`](../data/step3_z0.csv), [`data/step3_xi.csv`](../data/step3_xi.csv) and [`data/step3_profile.csv`](../data/step3_profile.csv).

**Weight at the edge** ($\delta = -0.30$). The last four columns are $L = 12$. $z_{\rm mid}$ is the weight between $0.8\Delta_+$ and $\Delta_+$, $z_{\rm above}$ the weight at or above $\Delta_+$; with $z_0$ they add to one.

| $V$ | $z_0$, $L = 8$ | $z_0$, $L = 12$ | $z_{\rm mid}$ | $z_{\rm above}$ | $\omega_{\rm edge}$ | next line | $0.8\Delta_+$ | $\Delta_+$ |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.731 | 0.713 | 0 | 0.287 | 0.023 | 0.828 | 0.480 | 0.600 |
| 0.5 | 0.706 | 0.700 | 0 | 0.300 | 0.024 | 0.982 | 0.609 | 0.762 |
| 1 | 0.673 | 0.675 | 0 | 0.325 | 0.029 | 1.116 | 0.735 | 0.919 |
| 2 | 0.589 | 0.588 | 0 | 0.412 | 0.061 | 1.342 | 0.983 | 1.228 |
| 3 | 0.503 | 0.464 | 0.163 | 0.373 | 0.162 | 1.502 | 1.253 | 1.566 |

**Decay length** ($\delta = -0.30$):

| $V$ | 0 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| $\xi_W$, $L = 8$ ($j = 0, 2$) | 3.1652 | 3.4101 | 3.5844 | 3.7330 | 3.6211 |
| $\xi_W$, $L = 12$ ($j = 0, 2, 4$) | 3.2161 | 3.2026 | 3.2042 | 3.2303 | 3.2201 |

Reading it:

- **One isolated edge line at every $V$** for weak end bonds; **none** for strong end bonds. The line is an edge effect.
- **Its weight falls** from 0.713 to 0.464 at $L = 12$.
- **Where the rest goes.** At every point except $L = 12$, $V = 3$, all of it lies at or above $\Delta_+$. At $L = 12$, $V = 3$, 0.163 sits in one line at $\omega = 1.502$, below $\Delta_+ = 1.566$. That line is the second state of the $N+1$ sector (Step 2's second addition level). At $V = 0$ it is the lowest bulk orbital of the open chain ($\omega = 0.828$): largest in the middle ($z_j = 0.33$ at $j = 5, 6$), 0.083 at the ends. By $V = 3$ its weight in the middle has fallen to 0.063 and at the ends risen to 0.163, the largest along the chain (profile figure, right). So it is not a plain bulk line at $V = 3$, and what it is has not been established. We do not call it bulk.
- **Sizes:** $L = 8$ and $12$ agree for $V \le 2$. At $V = 3$ the longer chain is lower by 0.039, so the $V = 3$ value is not converged in size.
- **$V = 0$:** 0.731 and 0.713 are the exact same-size values (check A4), not approximations to $1 - r^2 = 0.710$.
- **Profile:** at $L = 12$ the window weight falls with $V$ while its decay along the even sites barely changes ($\xi_W$ within 0.5% of 3.216). At $L = 8$ the fit uses only two sites and $\xi_W$ moves by up to 0.57; the two sizes differ by up to 0.5. $\xi_W(V)$ needs larger $L$.
- **Consistency with Step 2:** $\omega_{\rm edge}$ equals Step 2's addition gap ($0.0226, 0.0241, 0.0288, 0.0608, 0.162$ at $L = 12$), as it must, because the edge line is the lowest $N+1$ state.
- **Window:** the next line is above $0.9\Delta_+$ everywhere, so $z_0$ is identical for $f = 0.6$ to $0.9$.

**A defensible statement** (in the guide's form): *for $L \le 12$, the edge spectral weight of the weak-bond chain stays in a single isolated line up to $V = 3$, while its weight $z_0$ falls from 0.713 to 0.464 at $L = 12$; for $V \le 2$ this does not change between $L = 8$ and $12$.*

## 7. Limits and what was not tested

- **Sizes.** Only $L = 8, 12$. At $\xi = 3.23$ that is $L/\xi \le 3.7$. The guide says only size dependence can tell whether a dominant line survives or fragments.
- **The window at $L = 16, 18$.** At $V = 3$, Step 2's second addition level falls with $L$: 1.502 at $L = 12$ (outside the window for every $f \le 0.9$, $0.9\Delta_+ = 1.410$); 1.255 at $L = 16$ (outside $0.8\Delta_+ = 1.185$, inside $0.9\Delta_+ = 1.333$); 1.181 at $L = 18$ (inside $0.8\Delta_+ = 1.245$). At $L = 12$ it carries $z = 0.163$ at $j = 0$. If it still carries weight at $j = 0$ for $L = 16, 18$, $z_0$ there depends on $f$, and we must say so.
- **Odd $M$ at $L = 18$.** $L = 8, 12, 16$ have even $M = L/2$, so their rings have $\Delta_+ = 2|\delta| = 0.600$ at $V = 0$. $L = 18$ has $M = 9$ and $\Delta_+ = 0.685$, so its window edge is higher for a reason that has nothing to do with the open chain.
- **Odd sites.** $\xi_W$ is fitted on even sites only. At $V = 0$ the odd sites carry the right edge's tail (0.356 at the last site at $L = 12$ in $W_j$), so an odd-site fit would measure the far edge.
- **One site for the spectra.** Spectra and $z_0$ are at $j = 0$. Edge and centre spectra are not yet compared.
- **No broadening.** Lines are shown as they are.

## 8. Next steps

1. Lanczos spectral measure from $c_j^\dagger\Omega_N$ at $L = 16, 18$, first checked against the full spectrum at $L = 12$ with 100, 200, 400 vectors.
2. At $L = 16, 18$, $V = 3$: measure the weight of the second addition level at $j = 0$ and report $z_0$ for each $f$.
3. $z_0(V)$ and $\xi_W(V)$ at $L = 16, 18$, against the same-size $V = 0$ values 3.2277 ($j = 0$ to 6) and 3.2285 ($j = 0$ to 8).
4. Edge and centre spectra side by side.
5. Broadened spectra (widths 0.02, 0.01, 0.005) for presentation only.

## 9. Files and how to reproduce

| File | Contents |
|---|---|
| [`project_A.ipynb`](../project_A.ipynb) | Step 3 sections 3.1 to 3.5; Steps 1 and 2 unchanged |
| [`report.tex`](../report.tex) | Section 4 |
| [`figures/step3_spectra.pdf`](../figures/step3_spectra.pdf) | addition spectra, $L = 12$, $j = 0$, $\delta = \pm0.30$, $V = 0, 1, 2, 3$ |
| [`figures/step3_z0_xi.pdf`](../figures/step3_z0_xi.pdf) | $z_0(V)$ and $\xi_W(V)$ side by side, $L = 8, 12$ |
| [`figures/step3_profile.pdf`](../figures/step3_profile.pdf) | $z_j$ of the window line and of the next line along the chain, $L = 12$ |
| [`data/step3_z0.csv`](../data/step3_z0.csv) | one row per $(L,\delta,V)$: $z_0$ for each $f$, weight split at $\Delta_+$, $\omega_{\rm edge}$, next line and its weight, second $N+1$ level, $\Delta_+$, total weight |
| [`data/step3_xi.csv`](../data/step3_xi.csv) | $\xi_W$ per $(L, V)$ with fit sites, the same-size $V = 0$ value and the limit |
| [`data/step3_profile.csv`](../data/step3_profile.csv) | $W_j$, $z_j$ of the window line and of the next line, every site |
| [`data/step3_checks.csv`](../data/step3_checks.csv) | every Step 3 check with its worst deviation |

Run the whole notebook, or the Step 3 cells after Steps 1 and 2. Step 3 takes under a minute and prints `Step 3: 9/9 checks passed`.

## 10. Record of LLM assistance

The guide asks for a short record of substantive assistance and how it was verified.

- **Tool.** Claude (Anthropic), in claude.ai sessions on 2026-11-04 and 2026-11-05.
- **What it did.** Explained the observable, the sum rule and the window. Wrote `addition_spectrum`, which I typed in and ran myself after a line-by-line explanation. Wrote the check cells A1 to A9, `merge_lines`, `brute_force_addition`, the profile and $\xi_W$ cells, the figure and table cells, the notebook text, report Section 4 and this summary, built on my first runs. On 2026-11-05 it made the changes asked for in the group leader's review: same-size benchmarks for $z_0$ and $\xi_W$, the whole-profile check, the weight split at $\Delta_+$, the window notes for $L = 16, 18$, and moving the $L = 4$ comparison into the notebook as A9.
- **How it was verified.** A1 to A4 and A7, A8 compare against results that do not use the new code: the sum rule, a hand-derived $L = 2$ formula, and Step 1's exact solution at the same $L$. A9 compares with a separate construction that shares no code with Step 2, but it was written by the same LLM as `addition_spectrum`, so it guards against errors in assembling Step 2's pieces, not against a shared misunderstanding of the definition. The numbers in this summary are the ones the notebook prints.
- **What I still need to do.** Be able to derive the sum rule and the $L = 2$ check by hand, explain each line of `addition_spectrum`, and decide how to interpret the drop of $z_0$ and the line at $\omega = 1.502$ once larger sizes are available.
