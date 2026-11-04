# Step 3 summary: localisation and coherence (in progress)

Read this after [`step2.md`](step2.md), or at least its Sections 3 and 6. Each section answers one question. Sources: the guide is `QFT_student_research_guides.pdf`, Project A, p. 8; the code is Step 3 of [`project_A.ipynb`](../project_A.ipynb); the write-up is Section 4 of [`report.tex`](../report.tex).

## Short version

Steps 1 and 2 built the edge state and the many-body machine. Step 3 asks a physical question with it: **if you add one electron at the end of the chain, how much of it becomes one single eigenstate, the edge excitation?** That fraction is $z_0$.

**Done and verified, for $L = 8, 12$, $j = 0$:**

- The addition spectrum is built from Step 2's matrices and operators, and passes 5 checks (Section 5), including a hand-derived $L = 2$ result at $V \neq 0$ and the Step 1 orbitals at $V = 0$.
- Weak end bonds ($\delta = -0.30$): exactly one line lies inside the edge window at every $V \le 3$. Strong end bonds ($\delta = +0.30$): none.
- That line's weight falls with repulsion: $z_0 = 0.713, 0.700, 0.675, 0.588, 0.464$ at $V = 0, 0.5, 1, 2, 3$ ($L = 12$). The lost weight moves to bulk lines above $\Delta_+$.
- $L = 8$ and $L = 12$ agree for $V \le 2$ (to $0.006$ for $V > 0$) and differ at $V = 3$ ($0.503$ against $0.464$).
- $z_0$ does not depend on the window factor between 0.6 and 0.9.

**Not done yet:** the spatial profile $W_j$ and $\xi_W$ (guide item 4); $L = 16, 18$ with the Lanczos spectral measure (item 5); presentation broadening (item 1). Step 2 found that the window may stop being isolated at $L = 18$, $V = 3$, so that has to be tested there.

| Symbol | Meaning |
|---|---|
| $\Omega_N$ | the selected half-filled ground state (Step 2, Section 6) |
| $w_{jm}$ | probability that adding a particle at site $j$ gives exactly eigenstate $m$ of the $N+1$ sector |
| $\omega_m$ | $E_m^{N+1} - E_0^N$, the cost of that outcome |
| $z_{jm}$ | $w_{jm}/(1-\langle n_j\rangle) = 2w_{jm}$ |
| $\Delta_+$ | ring addition threshold at the same $(L,\delta,V)$ (Step 2) |
| $W$ | edge window, $0 \le \omega < f\Delta_+$, $f = 0.8$ |
| $z_0$ | total $z_{0m}$ inside $W$ |

---

## 1. What did Step 3 have to produce?

| Guide requirement (p. 8) | Where | Status |
|---|---|---|
| 1. Full addition spectrum at $L \le 12$; lines first; broadening 0.02, 0.01, 0.005 for presentation only | notebook 3.1, 3.3 | lines done at $L = 8, 12$; broadening not done |
| 2. $z_{jm}$ for isolated lines; summed weight for clusters | 3.1, 3.2 (A3 sums over clusters) | done |
| 3. Edge window from the ring, factor 0.8, varied 0.6 to 0.9; say so if no isolated window | 3.2 (A5), 3.4 | done at $L = 8, 12$: window isolated at every point |
| 4. $W_j$, fit $e^{-2j/\xi_W}$ on even and odd sites; compare edge and centre spectra | | not done |
| 5. $L = 16, 18$ with Lanczos (100, 200, 400 vectors, reorthogonalised) | | not done |
| Before moving on: one figure with $\xi_W(V)$ and $z_0(V)$ side by side | | $z_0(V)$ only |

## 2. What is being measured, and why?

Take the half-filled chain in its ground state and add one electron exactly on site $j$. Then ask which eigenstate it is in. Quantum mechanically there is no single answer: the result $c_j^\dagger|\Omega_N\rangle$ is a superposition of the eigenstates $|m, N+1\rangle$, and

$$w_{jm} = \big|\langle m,N+1|\,c_j^\dagger\,|\Omega_N\rangle\big|^2$$

is the probability of finding it in state $m$. Plotting $w_{jm}$ at $\omega_m$ gives the addition spectrum (guide Eq. A.5).

- **One tall line** means the electron operator creates mostly one eigen-excitation: the edge excitation is **coherent**.
- **Many small lines** mean the added electron spreads over many states, even if it stays near the edge.

That is why the guide separates coherence ($z_0$) from localisation ($\xi_W$): they answer different questions.

**At $V = 0$** (Step 1, Section 3.6): the added electron can only enter empty orbitals, with weight $\phi_n(j)^2$. The empty antibonding edge orbital takes $|\psi^L_0|^2/2$, so $z_0 = |\psi^L_0|^2 \to 1 - r^2 = 0.710$. Step 3 asks what repulsion does to this.

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

`edge_window` returns $f\Delta_+$ using Step 2's `ring_gap`.

## 5. How do we know it is right?

Step 2 verified the matrices and $c_j^\dagger$. What is new is how they are put together, so the checks target that.

| Id | What | Worst deviation | Independent of the new code because |
|---|---|---|---|
| A1 | sum rule $=\tfrac12$, every site, $L = 8, 12$, both $\delta$, every $V$ | $1.4\times10^{-13}$ | exact identity; fails if $\Omega$, $c^\dagger$ or the eigenbasis is wrong |
| A2 | $L = 2$: one line at $\omega = J_0 + V/2$, weight $\tfrac12$ | $8.9\times10^{-16}$ | derived by hand: $E_0 = -V/4 - J_0$, $E(\lvert 11\rangle) = V/4$; tests $V \neq 0$ |
| A3 | $V = 0$: lines at the empty orbital energies with weights $\phi_n(j)^2$, nothing elsewhere, every site | $4.9\times10^{-14}$ | uses only the $L\times L$ hopping matrix of Step 1 |
| A4 | $V = 0$: $z_0 = \lvert\psi^L_0\rvert^2$ | $2.8\times10^{-14}$ | Step 1's separated left edge state |
| A5 | one line in the window for $\delta = -0.30$, none for $+0.30$, every $V$ and $f$ | exact | a result as well as a check: the window is isolated |

During development the $L = 4$ spectrum was also compared with a separate brute-force calculation that shares no code with the notebook: weights $0.4201, 0.0799$ at $V = 0$ and $0.3723, 0.1194, 0.0067, 0.0016$ at $V = 2$. It agreed to all printed digits. That comparison was written by the same LLM, so by the guide's rule it supports A1 to A5 but is not counted as independent.

## 6. What do the results show?

Look at [`figures/step3_spectra.pdf`](../figures/step3_spectra.pdf). The full table is [`data/step3_z0.csv`](../data/step3_z0.csv).

| $V$ | $z_0$, $L = 8$ | $z_0$, $L = 12$ | $\omega_{\rm edge}$ ($L = 12$) | next line ($L = 12$) | $0.8\Delta_+$ ($L = 12$) |
|---|---|---|---|---|---|
| 0 | 0.731 | 0.713 | 0.023 | 0.828 | 0.480 |
| 0.5 | 0.706 | 0.700 | 0.024 | 0.982 | 0.609 |
| 1 | 0.673 | 0.675 | 0.029 | 1.116 | 0.735 |
| 2 | 0.589 | 0.588 | 0.061 | 1.342 | 0.983 |
| 3 | 0.503 | 0.464 | 0.162 | 1.502 | 1.253 |

Reading it:

- **One isolated edge line at every $V$** for weak end bonds; **none** for strong end bonds. The line is an edge effect.
- **Its weight falls** from 0.71 to 0.46 at $L = 12$. The weight goes to bulk lines above $\Delta_+$; at $V = 3$ the largest carries $z = 0.16$.
- **Sizes:** $L = 8$ and $12$ agree for $V \le 2$. At $V = 3$ the longer chain is lower by 0.039, so the $V = 3$ value is not converged in size.
- **$V = 0$:** 0.731 then 0.713, approaching $1 - r^2 = 0.710$.
- **Consistency with Step 2:** $\omega_{\rm edge}$ equals Step 2's addition gap ($0.0226, 0.0241, 0.0288, 0.0608, 0.162$ at $L = 12$), as it must, because the edge line is the lowest $N+1$ state.
- **Window:** the next line is above $0.9\Delta_+$ everywhere, so $z_0$ is identical for $f = 0.6$ to $0.9$.

**A defensible statement** (in the guide's form): *for $L \le 12$, the edge spectral weight of the weak-bond chain stays in a single isolated line up to $V = 3$, while its weight $z_0$ falls from 0.71 to 0.46; for $V \le 2$ this does not change between $L = 8$ and $12$.*

## 7. Limits and what was not tested

- **Sizes.** Only $L = 8, 12$. At $\xi = 3.23$ that is $L/\xi \le 3.7$. The guide says only size dependence can tell whether a dominant line survives or fragments.
- **One site.** Spectra and $z_0$ are at $j = 0$. The checks use every site, but the edge-versus-centre comparison is not yet made.
- **No $W_j$ or $\xi_W$.** Localisation is not yet measured, only coherence.
- **No broadening.** Lines are shown as they are.
- **The window at $L = 18$.** Step 2 found the second addition level inside $0.8\Delta_+$ at $L = 18$, $V = 3$. That has to be re-tested with Lanczos.
- **Odd $M$.** $L = 8, 12$ both have even $M = L/2$, so the rings have $\Delta_+ = 2|\delta|$ at $V = 0$. $L = 18$ will not.

## 8. Next steps

1. $W_j = \sum_{m\in W}w_{jm}$ for every $j$; fit $e^{-2j/\xi_W}$ separately on even and odd sites; plot edge and centre spectra.
2. Lanczos spectral measure from $c_j^\dagger\Omega_N$ at $L = 16, 18$, first checked against the full spectrum at $L = 12$ with 100, 200, 400 vectors.
3. The figure the guide asks for: $\xi_W(V)$ beside $z_0(V)$.
4. Broadened spectra (widths 0.02, 0.01, 0.005) for presentation only.

## 9. Files and how to reproduce

| File | Contents |
|---|---|
| [`project_A.ipynb`](../project_A.ipynb) | Step 3 sections 3.1 to 3.4; Steps 1 and 2 unchanged |
| [`report.tex`](../report.tex) | Section 4 added |
| [`figures/step3_spectra.pdf`](../figures/step3_spectra.pdf) | addition spectra, $L = 12$, $j = 0$, $\delta = \pm0.30$, $V = 0, 1, 2, 3$ |
| [`data/step3_z0.csv`](../data/step3_z0.csv) | one row per $(L,\delta,V)$: $z_0$ for each $f$, $\omega_{\rm edge}$, next line, $\Delta_+$, total weight |
| [`data/step3_checks.csv`](../data/step3_checks.csv) | every Step 3 check with its worst deviation |

Run the whole notebook, or the Step 3 cells after Steps 1 and 2. Step 3 takes about 25 seconds and prints `Step 3: 5/5 checks passed`.

## 10. Record of LLM assistance

The guide asks for a short record of substantive assistance and how it was verified.

- **Tool.** Claude (Anthropic), in a claude.ai session on 2026-11-04.
- **What it did.** Explained the observable, the sum rule and the window. Wrote `addition_spectrum`, which I typed in and ran myself after a line-by-line explanation. Wrote the check cells A1 to A5, the figure and table cells, the notebook text, report Section 4 and this summary, built on my first runs. Computed the $L = 4$ brute-force reference.
- **How it was verified.** A1 to A4 compare against results that do not use the new code: the sum rule, a hand-derived $L = 2$ formula, Step 1's orbitals, and Step 1's edge state. The numbers in this summary are the ones the notebook prints.
- **What I still need to do.** Be able to derive the sum rule and the $L = 2$ check by hand, explain each line of `addition_spectrum`, and decide how to interpret the drop of $z_0$ once larger sizes are available.
