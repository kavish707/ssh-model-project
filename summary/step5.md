# Step 5 summary: entropy resolved by particle number

Read this after [`step4.md`](step4.md). Each section answers one question.

**Sources:**
- the guide, `QFT_student_research_guides.pdf`, Project A, p. 10;
- the code, Step 5 of [`project_A.ipynb`](../project_A.ipynb), sections 5.1-5.10;
- the write-up, the Step 5 section of [`report.tex`](../report.tex).

## Short version

Steps 3 and 4 asked where an electron added at the edge goes and how coherent it stays. Step 5 asks a quantum-information question: **how much entanglement does it create between the edge region and the rest, and is that entanglement usable at fixed particle number ($E_{\rm acc}$), or is it only uncertainty about how many particles are in the region ($H_N$)?**

**Done and verified:**
- **Grid.** Static states at $L=8,12,16,18$; time evolution at $L=16$ ($t\le8$) and $L=20$ ($t\le16$); $\delta=\pm0.30$, $V=0$ to $3$.
- **Cuts.** $A=\{0,1,2,3\}$ (edge) and the half chain.
- **19 checks** (Section 5): the guide's two test states and three more; an independent full-Fock-space route at $L=4,6$; the correlation-matrix result at $V=0$ for every state and time; particle-hole symmetry, solver independence and Krylov tolerance at $V\neq0$; a negative control for the state selection.

**Main results** ($L=16$, edge cut, $\delta=-0.30$):
- **The edge pair carries about $\ln2$ across every cut.** In the selected ground state the two ends share one fermion. A single shared fermion is pure number entropy. $E_{\rm acc}=0.29$ appears only because a second uncertain fermion, the strong dimer cut in two, also crosses the edge cut.
- **Repulsion converts number entropy into accessible entanglement in the ground state:** $H_N$ goes $1.006\to0.746$ and $E_{\rm acc}$ goes $0.294\to0.474$ from $V=0$ to 3.
- **Adding the electron at the edge *lowers* the entanglement of the edge region:** $\Delta S_A=-0.30$ ($V=0$) to $-0.41$ ($V=3$) at $t=0$. The coherent part (the filtered edge eigenstate) fills the shared edge pair. It has the entropies of a left-loaded ground state, with almost no accessible entanglement ($E_{\rm acc}=0.011$).
- **In time, inside the certified intervals,** the entropies recover only partly: late $\Delta S_A=-0.13$ to $-0.24$. With $V$ the persistent change moves from $H_N$ into $E_{\rm acc}$. $L=16$ and $20$ agree to 0.011.
- **Controls raise the entropy instead.** $\delta=+0.30$ adds pure number entropy (up to $+0.27$, $\Delta E_{\rm acc}\le0.02$). Centre injection raises $E_{\rm acc}$ at the edge cut by 0.2-0.4 once the electron arrives.
- **Guide item 5:** $E_{\rm acc}(t)$ tracks the edge charge $Q_4(t)$ (correlation $-0.7$), not $P_{\rm return}(t)$ ($|r|\le0.25$). Coherence and accessible entanglement are not lost together: the late $P_{\rm return}$ falls from 0.54 to 0.22 while the late $E_{\rm acc}$ rises from 0.19 to 0.35.

| Symbol | Meaning |
|---|---|
| $A$, $\ell$ | left block $\{0,\dots,\ell-1\}$; edge cut $\ell=4$, half cut $\ell=L/2$ |
| $C^{(q)}$ | block of the coefficient matrix with $q$ particles in $A$ |
| $s_{q\alpha}$ | Schmidt values (singular values) of $C^{(q)}$ |
| $p_q$ | probability of $q$ particles in $A$ |
| $S_q$ | entanglement entropy of the normalised fixed-$q$ part |
| $S_A$, $H_N$, $E_{\rm acc}$ | total, number and accessible entropy, $S_A=H_N+E_{\rm acc}$ |
| $\Psi(t)$ | Step 4's edge-addition state, $c_0^\dagger\Omega_N/\sqrt{1-\langle n_0\rangle}$, evolved |
| filtered | Step 4's projection of $\Psi(0)$ onto the edge window: the edge eigenstate alone |
| $\Delta X$ | $X(\Psi)-X(\Omega_N)$ at the same $(L,\delta,V)$ and cut |

---

## 1. What did Step 5 have to produce?

| Guide item (p. 10) | Where | Status |
|---|---|---|
| 1. Test states: $(\lvert10\rangle+\lvert01\rangle)/\sqrt2$ gives $S_A=H_N=\ln2$, $E_{\rm acc}=0$; the swapped two-site state gives $H_N=0$, $E_{\rm acc}=\ln2$ | 5.2, A1 | exact, plus three more states |
| 2. All three quantities in the ground state and the edge-addition state, edge cut $\{0,1,2,3\}$ and half cut | 5.4, 5.6 | $L=8$-$18$ static, $L=16,20$ in time |
| 3. Changes relative to the ground state, absolute entropies kept | 5.4, 5.8 (`step5_dynamics`, `step5_dynamics_absolute`) | done |
| 4. $V=0$: correlation-matrix $S_A$; $S_A-H_N-E_{\rm acc}=0$ | 5.3 (B1-B4), A2 | done for $\Omega_N$, $\Psi(0)$ and $\Psi(t)$; $p_q$ also checked |
| 5. $E_{\rm acc}(t)$ beside $P_{\rm return}(t)$ and $Q_{\rm edge}(t)$ for each $V$ | 5.6-5.8 | done, with late averages and correlations |
| Before moving on: edge region vs the rest is not the two edges vs each other | intro, 5.10 | stated; no two-edge entropy computed |
| Extras | 5.4-5.5 | cut-position profile, filtered state, centre and $\delta=+0.30$ controls, left-loaded negative control, $L=20$ to $t=16$ |

## 2. What is being measured, and why?

The entanglement entropy $S_A$ measures how much is unknown about region $A$ when the whole chain is in a known pure state. Take the Schmidt decomposition (the SVD of the coefficient matrix $C_{ab}$); then $S_A=-\sum s^2\ln s^2$.

For particles there is a catch:
- **A single particle in a superposition of "in $A$" and "in $B$" already gives $S_A=\ln2$.** Local operations that conserve particle number cannot use that entanglement.
- **Total number is conserved,** so $C$ is block diagonal in $q=N_A$, and the entropy splits exactly:

$$S_A=\underbrace{-\sum_qp_q\ln p_q}_{H_N}+\underbrace{\sum_qp_qS_q}_{E_{\rm acc}} .$$

**$H_N$** answers "how many particles are in $A$?". **$E_{\rm acc}$** is what remains after that number is known: correlations in *how* the particles are arranged.

**The two-line proof.** $\rho_A=\bigoplus_qp_q\rho_q$ has eigenvalues $p_q\lambda_{q\alpha}$, and $-\sum p_q\lambda\ln(p_q\lambda)=-\sum p_q\ln p_q-\sum_qp_q\sum_\alpha\lambda\ln\lambda$.

**The guide's two test states show the extremes:**
- one shared fermion gives $H_N=\ln2$, $E_{\rm acc}=0$;
- two particles whose site positions are swapped between regions give $H_N=0$, $E_{\rm acc}=\ln2$.

The extra test state used here is **two fermions, each shared across the cut**: $(S,H,E)=(2\ln2,\tfrac32\ln2,\tfrac12\ln2)=(1.39,1.04,0.35)$. It is the key to the ground-state result, because $E_{\rm acc}$ needs at least two uncertain modes crossing the cut.

## 3. Why no fermion signs?

Basis states are $|n\rangle=(c_0^\dagger)^{n_0}\cdots(c_{L-1}^\dagger)^{n_{L-1}}|0\rangle$. For a left block, all of $A$'s operators already stand before $B$'s, so $|n\rangle=|a\rangle_A|b\rangle_B$ exactly. In code $a=s\,\&\,(2^\ell-1)$ and $b=s\gg\ell$. A region that is not a left block would need a reordering sign; none is used here.

## 4. How does the code compute it?

**`cut_plan` and `sector_entropies`** (notebook 5.1):
- `cut_plan(L, N, ell)` splits each basis integer into $a$ and $b$, groups by $q=$ popcount$(a)$, and records each state's row and column in $C^{(q)}$.
- `sector_entropies(psi, L, N, ell)` fills each block, takes its singular values, and returns $p_q$, $S_q$, $S_A$ (from the whole spectrum, not as $H_N+E_{\rm acc}$, so the identity is a real test), $H_N$, $E_{\rm acc}$, and the weight dropped below $p_q<10^{-14}$.
- The blocks are tiny: at most $252\times252$, at the half cut of $L=20$.

**Independent routes:**
- **`fock_route`** builds $c_i=Z\otimes\cdots\otimes a\otimes1\cdots$ as $2^L\times2^L$ Kronecker products, then $H$, $\Omega_N$, $c_0^\dagger\Omega_N$, and $\rho_A$ by reshaping, with site 0 the leading factor. It uses no Step 2 code.
- **`gaussian_entropies` and `slater_C`** handle $V=0$ through the correlation matrix. In the eigenbasis of $C_A$, the state of $A$ is a product of independent modes, each filled with probability $\nu_k$. That gives $S_A=-\sum[\nu\ln\nu+(1-\nu)\ln(1-\nu)]$, and $p_q$ as the distribution of a sum of independent coin flips (a repeated convolution). Adding a particle adds the orbital $\chi=(1-C)e_j/\sqrt{1-C_{jj}}$ (Step 4), and $\chi(t)=e^{-iht}\chi(0)$.

**`loaded_states`** gives the two combinations of the near-degenerate $N$-sector doublet that diagonalise the left-minus-right imbalance. They are a negative control.

**`evolve_entropies`** repeats Step 4's Krylov propagation exactly (`krylov_step`, the same tolerance, the same step splitting). It evaluates the entropies, $p_q$ and $S_q$ at both cuts at every $\Delta t=0.05$, so the result can be compared with Step 4's saved trajectory.

## 5. How do we know it is right?

| Id | What | Worst | Independent because |
|---|---|---|---|
| A1 | guide test states and three more | $2.2\times10^{-16}$ | exact by hand |
| A2 | $\sum p_q=1$, $S_A-H_N-E_{\rm acc}=0$, every state, cut and time (three checks) | $4.9\times10^{-15}$ | bookkeeping |
| A3 | blocked SVD = one SVD of the whole $C_{ab}$ | $1.5\times10^{-13}$ | no weight lost between blocks |
| A4 | $L=4,6$, every cut: $S_A$, $p_q$, $S_q$ of $\Omega_N$ and $c_0^\dagger\Omega_N$ = Fock-space partial trace; $\rho_A$ has no elements between sectors | $1.1\times10^{-14}$ | separate construction |
| B1 | $V=0$: $S_A(\Omega_N)$ = correlation matrix, every cut, $L\le18$ | $1.4\times10^{-13}$ | $L\times L$ matrices only |
| B2 | $V=0$: $p_q$, $H_N$, $E_{\rm acc}$ of $\Omega_N$ = independent modes | $1.5\times10^{-12}$ | tests the split, not only $S_A$ |
| B3 | $V=0$: the same for $c_j^\dagger\Omega_N$, edge and centre | $1.6\times10^{-12}$ | Slater determinant $C+\chi\chi^\dagger$ |
| B4 | $V=0$: the same for $\Psi(t)$ at every time, $L=16$ | $2.1\times10^{-12}$ | $\chi(t)=e^{-iht}\chi(0)$ |
| C1 | $\Omega_N$: $p_q=p_{\ell-q}$, $S_q=S_{\ell-q}$, every $V$ | $1.1\times10^{-11}$ | exact particle-hole symmetry |
| C2 | $\sum q\,p_q=\sum_{j<\ell}\langle n_j\rangle$ | $9.3\times10^{-15}$ | Step 2 densities, Step 4 $\Delta n_j(0)$ |
| C3 | the evolution reproduces Step 4's $G(t)$, $\Delta n_j(t)$ | 0 | same state as Step 4 |
| C4 | Krylov tolerance halved and $10^{-13}$ | $1.1\times10^{-13}$ | convergence |
| C5 | $L=12$: Lanczos ground state = dense lowest vector | $2.6\times10^{-15}$ | solver independence |
| C6 | negative control: a left-loaded doublet state changes $S_A$ by 0.59; at $V=0$ it equals its Slater determinant to $10^{-14}$ | – | shows that the state selection matters |

B1-B3 run at $L=8$ first, and the larger sizes run only if they pass.

## 6. What do the results show?

**Figures:**
- `step5_static`: entropies against $V$, both cuts;
- `step5_cut_profile`: entropy against the cut position;
- `step5_dynamics`: changes in time beside $P_{\rm return}$ and $Q_4$;
- `step5_dynamics_absolute`;
- `step5_sectors`: $p_q$ and $S_q$.

**Ground state and $t=0$** ($(S_A,H_N,E_{\rm acc})$, edge cut, $L=16$):

| | $V=0$ | $V=1$ | $V=3$ |
|---|---|---|---|
| $\Omega_N$, $\delta=-0.30$ | 1.300, 1.006, 0.294 | 1.310, 0.970, 0.340 | 1.220, 0.746, 0.474 |
| $\Psi(0)$ | 1.002, 0.853, 0.150 | 0.971, 0.806, 0.165 | 0.811, 0.625, 0.186 |
| filtered | 0.714, 0.703, 0.011 | 0.698, 0.687, 0.011 | 0.630, 0.620, 0.010 |
| left-loaded $\Omega$ | 0.714, 0.703, 0.011 | 0.698, 0.687, 0.011 | 0.586, 0.573, 0.013 |
| $\Omega_N$, $\delta=+0.30$ | 0.201, 0.197, 0.004 | 0.172, 0.144, 0.028 | 0.224, 0.091, 0.133 |

**Late averages** (edge cut, edge injection, $\delta=-0.30$):

| $V$ | 0 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| $\Delta S_A$, $L=16$ | -0.235 | -0.215 | -0.184 | -0.125 | -0.164 |
| $\Delta H_N$, $L=16$ | -0.127 | -0.113 | -0.095 | -0.056 | -0.036 |
| $\Delta E_{\rm acc}$, $L=16$ | -0.108 | -0.102 | -0.090 | -0.069 | -0.128 |
| $L=20$ ($\Delta S$, $\Delta H$, $\Delta E$) | -0.241, -0.125, -0.116 | | -0.192, -0.091, -0.101 | | -0.159, -0.027, -0.131 |
| $P_{\rm return}$, $L=16$ | 0.540 | 0.522 | 0.487 | 0.366 | 0.220 |

**Reading it:**
- **The cut-position profile** has two parts. The zigzag comes from the dimers: a cut through a strong bond splits a dimer and shares a fermion across the cut, about $\ln2$. On top of that, $\delta=-0.30$ sits $0.60$-$0.68\approx\ln2$ above $\delta=+0.30$ at the same bond type. That offset is the edge fermion shared between the two ends.
- **The ground state's $E_{\rm acc}$ is in the $q=2$ sector** and needs both shared fermions, the cut dimer and the edge pair. The two-fermion test state predicts $(1.39,1.04,0.35)$, close to $(1.30,1.01,0.29)$.
- **Adding at the edge fills the edge pair.** The filtered state equals a left-loaded ground state. The coherent part of the electron removes the shared uncertainty, so the entanglement drops.
- **In time** the incoherent part crosses the cut and partly restores the entropy. The edge region stays less entangled than before. With $V$, the persistent change shifts into $E_{\rm acc}$.
- **Controls:** a trivial edge or centre injection *raises* the entropy instead.
- **Guide item 5:** over the certified interval, $E_{\rm acc}(t)$ follows $Q_4(t)$ (correlation $-0.68$ to $-0.75$ for $V\le1$) and not $P_{\rm return}(t)$. Coherence falls with $V$ while the accessible entanglement does not.

**A defensible statement:** *on $L\le20$ and up to $V=3$, inside the certified intervals, adding an electron at the weak-bond edge lowers the entanglement between the first four sites and the rest, by 0.30-0.41 at once and 0.13-0.24 after the incoherent part has left. The coherent part does this by filling the edge pair the two ends share. With repulsion, the change moves from number entropy into accessible entanglement, and $E_{\rm acc}(t)$ tracks the edge charge, not the return probability.*

## 7. Limits

- **Absolute entropies depend on which bond the cut crosses.** Compare only the same bond type. The half cut is at a weak bond for $L=18$ and at a strong one otherwise.
- **The half-cut changes are not size-converged.** The edge cut is: $L=16$ and $18$ agree to 0.013 statically, and $L=16$ and $20$ to 0.011 in time.
- **$E_{\rm acc}$** is operational under local number conservation, not a general distillable rate. A difference of entropies is not an entanglement measure.
- **The correlation coefficients** describe co-movement over a few time units.
- **All values refer to the selected ($C=+1$) ground state.** The negative control shows the left-loaded state differs by 0.59.
- **No entanglement between the two edges is computed:** their reduced state is mixed.

## 8. Files and how to reproduce

| File | Contents |
|---|---|
| [`project_A.ipynb`](../project_A.ipynb) | Step 5 sections 5.1-5.10 after Step 4 |
| [`report.tex`](../report.tex) | Step 5 section |
| `figures/step5_static.pdf`, `step5_cut_profile.pdf`, `step5_dynamics.pdf`, `step5_dynamics_absolute.pdf`, `step5_sectors.pdf` | figures |
| `data/step5_static.csv` | every $(L,\delta,V)$, state and cut: entropies, changes, $p_q$, the bond crossed, check columns |
| `data/step5_profile.csv`, `step5_negative_control.csv` | cut profile; selected vs loaded states |
| `data/step5_dynamics.csv` | every run and time: $P_{\rm return}$, $G$, $Q_2$, $Q_4$, $\Delta n_j$, entropies, changes, $p_q$ and $S_q$ at both cuts |
| `data/step5_late.csv` | late averages and correlation coefficients |
| `data/step5_checks.csv`, `step5_run_info.json` | checks and run settings |

**To reproduce:** run the whole notebook. Step 5 needs Steps 1-4 in memory, takes about 3 minutes after Step 4 (the whole notebook about 15), and prints `Step 5: 19/19 checks passed`.

## 9. Record of LLM assistance

*(to be filled in by Kavish)*
