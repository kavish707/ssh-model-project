# Step 4 summary: watch an edge particle evolve

Personal study copy, kept outside the repo. Read it after `ssh-model-project/summary/step3.md`, or at least its Sections 2 and 6. Step 4 is the time-domain version of Step 3.

**Sources:**
- the guide, `QFT_student_research_guides.pdf`, Project A, Step 4;
- the code, Section 4 of `ssh-model-project/project_A_steps1-4.ipynb` (committed to `main` by Somnath on 17 Nov, cells 85–121);
- Somnath's short report, `projectA_step4_short_report.pdf`.

The numbers below are the ones printed in the saved notebook outputs. Somnath ran the notebook, and I have not rerun it. The same physics in his earlier script version (`somnath-project` branch) passed all 16 of its tests when I ran them on our environment.

---

## Short version

Step 3 decomposed "one electron added at the edge" into eigenstates of the $N+1$ sector. Step 4 **lets that state evolve in time** and asks three things:
- does the added electron stay at the edge;
- does it stay one coherent excitation;
- how long can a finite chain answer that before reflections from the far end spoil the result?

**Done and verified:**
- Sizes $L=8,12,16$ (plus $18$ and $20$ only to fix the clean time window), $\delta=\pm0.30$, $V=0,0.5,1,2,3$, $0\le t\le8$; $t\le16$ at $L=18,20$.
- 38/38 Step 4 checks pass (159/159 over the whole notebook).

**Main results:**
- **Part of the electron stays, and Step 3 predicts how much.** Inside the clean time window, the return probability settles near $z_0(\mathcal W)^2$, the square of Step 3's edge-line weight:
  - $L=20$: $0.499$ vs $0.504$ ($V=0$), $0.460$ vs $0.464$ ($V=1$), $0.228$ vs $0.208$ ($V=3$).
  - The edge charge on four sites settles near $z_0q_4$, the charge carried by the edge eigenstate: $0.327$ vs $0.325$ and $0.280$ vs $0.288$ at $V=0,1$.
- **Repulsion reduces what stays,** through $z_0$: late $P_{\rm return}$ falls from $0.54$ to $0.22$ at $L=16$.
- **The controls lose the electron.** With the opposite dimerisation ($\delta=+0.30$), $P_{\rm return}\approx0.01$–$0.02$. With centre injection it is $0.07$–$0.10$. The effect belongs to the topological edge.
- **About a quarter of the added charge is at the far end already at $t=0$** ($0.25$ at $V=0$, $0.34$ at $V=3$). This is not transport. In the selected ground state the two ends share one fermion (Step 1's bonding orbital), so the coherent part of the electron is half at each end.
- **No beating between the two edges.** At $L=8$, evolved exactly to $t=300$, there is no signal at the splitting frequency. At $L\ge16$ the would-be period (480 to 3300) is far beyond any run anyway.
- **Clean time windows**, before a reflection returns:
  - $L=16$: the whole run ($t\le8$) at $V\le2$, $t<4.90$ at $V=3$;
  - $L=18$: $t\le16$ ($V=0$), $14.45$ ($V=1$), $6.50$ ($V=3$);
  - $L=8$ and $12$: none, because they already differ from the next size at $t=0$.
- **The filtered excitation** (only the edge-window part of $\Psi(0)$) is an eigenstate, so it never changes. It shows why the guide insists on the unfiltered electron: filtering removes the decay by construction.

**Open:**
- why the edge charge at $V=3$ stays above both predictions;
- the large effective "return speeds" at $V=3$;
- what Step 3's second line is. It carries even more edge charge than the edge line.
- **Repo status:** Step 4 is not in `report.tex` yet, there are two notebooks, and the figure and data files are not consistent with the notebook (Section 8).

| Symbol | Meaning |
|---|---|
| $\Omega_N$, $E_0^N$ | selected half-filled ground state and its energy (Step 2 rule: $C=+1$) |
| $\Psi(0)$ | $c_j^\dagger\Omega_N/\sqrt{1-\langle n_j\rangle}$, the electron just after it is added at site $j$ |
| $\omega_m$, $z_{jm}$ | Step 3's lines: cost $E_m^{N+1}-E_0^N$, and share $\lvert\langle m\rvert\Psi(0)\rangle\rvert^2$ |
| $G(t)$ | return amplitude $\langle\Psi(0)\rvert\Psi(t)\rangle$ |
| $P_{\rm return}$ | $\lvert G(t)\rvert^2$ |
| $\Delta n_j(t)$ | density at $j$ in $\Psi(t)$ minus that in $\Omega_N$ |
| $Q_{\rm edge}=Q_{\ell_e}$ | $\sum_{j<\ell_e}\Delta n_j$, with $\ell_e=2$ or $4$ |
| $z_0(\mathcal W)$ | Step 3's window weight at $j=0$ (here always one line) |
| $q_4$ | excess charge on sites 0–3 carried by the edge eigenstate itself |
| $t_{\rm fs}$ | end of the finite-size-free ("certified") interval |
| $s$ | static offset: how much two chain lengths already differ at $t=0$ |

---

## 1. What did Step 4 have to produce?

| Guide item | Where in the notebook | Status |
|---|---|---|
| Prepare $\Psi(0)$ (Eq. A.6), keep the full state; compute $\Delta n_j$, $Q_{\rm edge}$ for $\ell_e=2,4$, $P_{\rm return}$ (A.7–A.9) | 4.1, 4.2 | done |
| 1. Krylov evolution at $L=8,12,16$, $\Delta t=0.05$, $t\le8$; full diagonalisation at $L=8$; halve the tolerance | 4.2, 4.3 | done; full diagonalisation also at $L=12$; tolerance scan $10^{-3}$ to $5\times10^{-11}$ |
| 2. Repeat with injection near the centre and with the opposite dimerisation | 4.2 (every run) | done |
| 3. Check norm, energy and $\sum_j\Delta n_j=1$; plot the density profile and $P_{\rm return}$ | 4.2, 4.10, 4.11 | done, plus closed-form energy moments |
| 4. Inspect reflections and two-edge beating; stop before the return reaches the observation region; extend $t$ only after increasing $L$ | 4.4–4.6 | done: certified interval for every run, $t\le16$ only at $L=18,20$, beating checked to $t=300$ |
| 5. Project onto the edge window; call it a filtered excitation; report $z_0(\mathcal W)$; compare with the unfiltered electron | 4.2, 4.9 | done |

---

## 2. The physics: what happens after you add an electron at the edge

### 2.1 The prepared state is a superposition, and Step 3 already knows it

$\Psi(0)$ is not an eigenstate of $H_{N+1}$. Expand it in eigenstates:

$$|\Psi(0)\rangle=\sum_m a_m|m\rangle,\qquad |a_m|^2=z_{0m},\qquad \sum_m z_{0m}=1 .$$

The $z_{0m}$ are exactly Step 3's normalised line weights. Each component only picks up a phase:

$$|\Psi(t)\rangle=\sum_m a_m e^{-i\omega_m t}|m\rangle,\qquad G(t)=\sum_m z_{0m}\,e^{-i\omega_m t}.$$

So **Step 4 contains no new information beyond the full Step 3 spectrum and eigenstates**. It shows the same physics in time, and it is also a strong cross-check: $G(t)$ from time evolution must equal Step 3's spectral sum, as a complex number, phase included. That is why $E_0^N$ is kept in the propagator $e^{-i(H_{N+1}-E_0^N)t}$.

### 2.2 Why the return probability settles at about $z_0^2$

The edge line, with weight $z_0$, sits far below the bulk threshold $\Delta_+$. The other lines are bulk lines spread over a wide band.

$$P_{\rm return}(t)=\Big|z_0e^{-i\omega_0t}+\sum_{\rm bulk}z_{0m}e^{-i\omega_mt}\Big|^2 .$$

- **The bulk part dephases.** Many frequencies run out of phase, so their sum averages to nearly zero. In an infinite chain the bulk lines form a continuum and that part decays to zero.
- **The edge part never decays.** It is one eigenstate, and it is below the continuum, so there is nothing at its energy for it to decay into.

Averaged over time, $\overline{P}=\sum_m z_{0m}^2=z_0^2+\sum_{\rm bulk}z_{0m}^2$. The bulk sum shrinks as $L$ grows because the bulk weight is shared among more lines. So

$$\overline{P_{\rm return}}\;\longrightarrow\;z_0^2 .$$

The physical reading: **the added electron behaves like a stable edge quasiparticle with weight $z_0$** (the quasiparticle residue) **plus a part $1-z_0$ that leaves into the bulk.** Repulsion lowers $z_0$ (Step 3), so less of the electron stays.

### 2.3 Charge: what $\Delta n_j$ and $Q_{\rm edge}$ measure

- **$\sum_j\Delta n_j(t)=1$ at all times.** $\Psi$ has $N+1$ particles and $\Omega$ has $N$. This is a check on the code.
- **$\Delta n_j$ can be negative.** The background rearranges around the added electron; with repulsion, neighbours are pushed away. This is why the profiles alternate in sign at $V=3$.
- **$Q_{\rm edge}$ is an excess charge, not a probability.** It can exceed $P_{\rm return}$ or fall below it.

### 2.4 The key subtlety: the edge excitation is half at each end

At $t=0$ the excess charge is not all at site 0:
- $\Delta n_0(0)=1-\langle n_0\rangle=\tfrac12$;
- $\Delta n_{L-1}(0)\approx0.25$ at $V=0$, at the **far** end, 15–19 sites away.

**Why:**
- At $V=0$ the end pair (Step 1) is one bonding orbital $\phi_-$, filled, and one antibonding orbital $\phi_+$, empty. Each has half its weight at each end.
- So in $\Omega_N$ the two ends **share** one fermion, $\tfrac12$ each.
- The coherent part of the added electron goes into $\phi_+$. Once both end modes are full, each end has gained $\tfrac12$.

**So the edge eigenstate carries half its excess charge at each end.** The notebook confirms this: the edge line's excess on the left half is exactly $0.5$ at every $V$, with $q_4=0.458$ of it on sites 0–3 at $V=0$.

**The $V=0$ formula.** By Wick's theorem, $\Delta n_{L-1}(0)=2|C_{0,L-1}|^2$, where $C_{ij}=\langle c_i^\dagger c_j\rangle$. The bonding orbital gives $|C_{0,L-1}|\approx(1-r^2)/2$, so

$$\Delta n_{L-1}(0)\to\frac{(1-r^2)^2}{2}=0.252\quad(\delta=-0.30).$$

The same number is $z_0\times\tfrac{1-r^2}{2}$: the edge-line weight times that eigenstate's charge on the last site.

**Consequences:**
- $Q_4(0)$ is only $0.65$ at $V=0$ and $0.77$ at $V=3$. The rest is mostly at the far end ($0.32$ and $0.21$ on the last four sites).
- **The most charge that can stay on the left four sites is about $z_0q_4\approx z_0/2$, not $z_0$.**
- This is the Step 1/Step 2 state-selection question again. In a left- or right-loaded combination of the near-degenerate doublet, the far-end response is below $10^{-3}$ at $L\ge16$. The selected state ($C=+1$) is the exact finite-$L$ ground state and the one Step 3 used, so it is the consistent choice. The half-and-half picture is a property of that state.
- $\Delta n_{L-1}(0)$ does not shrink with $L$ (0.249 at $L=20$, $V=0$). It is a correlation of the ground state, not a signal that has travelled.

### 2.5 What a finite chain can and cannot tell you

**Reflections.** The part of the electron that leaves the edge travels to the far end, reflects, and comes back. After that, the chain shows a finite-size revival, not a long lifetime. The guide therefore asks to stop interpreting before the reflection returns.

**The infinite-time average is meaningless on a finite chain.** Every $N+1$ eigenstate is reflection symmetric, so the infinite-time (diagonal-ensemble) average puts **exactly $\tfrac12$ of the excess charge in each half** (check E2). Wait long enough and the edge "keeps" half the electron for a trivial reason. Only the certified interval carries physics.

**Two-edge beating.** On a finite chain the $N$-sector doublet is split by $E_1^N-E_0^N\sim e^{-L/\xi}$. If the added electron reached two $N+1$ lines separated by that amount, charge would slosh between the ends with period $2\pi/(E_1^N-E_0^N)$. It doesn't: in the $N+1$ sector both end modes are filled, so there is one edge line, not a split pair. The partner weight is exactly 0.

---

## 3. Exact results used as checks: derivations

**$t=0$ identities, for any $V$.** For $k\neq j$, $n_k$ commutes with $c_j$, so

$$\langle\Psi(0)|n_k|\Psi(0)\rangle=\frac{\langle n_k(1-n_j)\rangle}{1-\langle n_j\rangle}\;\Rightarrow\;\Delta n_k(0)=-\frac{\langle\delta n_j\,\delta n_k\rangle_c}{1-\langle n_j\rangle}.$$

For $k=j$, $n_jc_j^\dagger=c_j^\dagger$ gives $\Delta n_j(0)=1-\langle n_j\rangle=\tfrac12$. So the $t=0$ profile is a ground-state density–density correlation (check B2).

**Energy moments and the initial decay (check B4).** From $[H,c_0^\dagger]=-J_0c_1^\dagger+Vc_0^\dagger(n_1-\tfrac12)$ and $\langle n_j\rangle=\tfrac12$:

$$m_1=\langle H-E_0^N\rangle_{\Psi(0)}=2\big[J_0\langle c_0^\dagger c_1\rangle-V\langle\delta n_0\delta n_1\rangle_c\big],\qquad m_2=\langle(H-E_0^N)^2\rangle_{\Psi(0)}=J_0^2+2J_0V\langle c_0^\dagger c_1\rangle+\tfrac{V^2}{4}.$$

Expanding $G(t)=\langle e^{-i(H-E_0)t}\rangle$ gives $P_{\rm return}=1-(m_2-m_1^2)\,t^2+O(t^4)$. **The initial drop is set by the energy spread $\sigma=\sqrt{m_2-m_1^2}$ of the added electron.** At $L=16$:

| | $V=0$ | $V=3$ |
|---|---|---|
| $\delta=-0.30$: $\sigma$ | 0.596 | 1.383 |
| $\delta=+0.30$: $\sigma$ | 0.354 | 0.662 |

- **$V=3$ drops faster** because the electron has a larger energy spread.
- **$\delta=+0.30$ starts slower.** Sites 0–1 form a strong dimer whose bonding orbital is full, so the added electron sits mostly in the antibonding dimer orbital at energy about $J_s$, a narrow distribution. Later it leaks out completely through the weak bond.

**Pauli blocking at $V=0$ (one-body route, check B3).** Write $C$ for the filled-orbital correlation matrix with the bonding orbital, `filled_correlations`. Then $c_j^\dagger\Omega$ creates one particle in

$$\chi(0)=\frac{(1-C)e_j}{\sqrt{1-C_{jj}}},$$

which is $e_j$ with its occupied part removed. Then $\chi(t)=e^{-iht}\chi(0)$, $G(t)=\chi(0)^\dagger\chi(t)$ and $\Delta n_j(t)=|\chi_j(t)|^2$. The whole Step 4 calculation at $V=0$ therefore needs only the $L\times L$ matrix, which makes it an independent check.

---

## 4. How the code computes it

**Inputs** (notebook 4.1; Steps 1–3 are not rebuilt):
- **$\Omega_N$:** Step 2's saved `STATES`. Only $L=20$ is new; it uses the same `lowest` and `select_ground_state`.
- **Operators:** Step 2's `c_dag`, which carries the sign string, and `hamiltonian`.
- **Step 2 grid:** $E_0^N$, the neutral gap and $\Delta_+$ from `data/step2_grid.csv`.
- **Step 3 lines:** `site_spectra` ($L\le12$) and the Lanczos measure ($L=16,18$).

**`prepare(L, delta, V, site)`.**
- Builds $\Psi(0)$ and $H_{N+1}$.
- Stores $\langle n_k\rangle$ and $\langle\delta n_j\delta n_k\rangle_c$ for the $t=0$ check.
- The norm before normalising is $1-\langle n_j\rangle$, the Step 3 sum rule.
- The centre site is $j_c=2\lfloor L/4\rfloor$, the even site nearest the centre (site 8 at $L=16$).

**`evolve` (Krylov, short-iterative Lanczos).** For each output step $\Delta t=0.05$:
1. Build orthonormal vectors $q_0=\psi/|\psi|,q_1,\dots$ by applying $H-E_0^N$, re-orthogonalising against all earlier vectors twice (as in Step 3's `lanczos_spectrum`).
2. In that basis $H$ is a small tridiagonal matrix $T_m$, and $e^{-i(H-E_0)\Delta t}\psi\approx|\psi|\,Q_m e^{-iT_m\Delta t}e_1$.
3. Add vectors until the error estimate $|\psi|\,\beta_m|e_m^Te^{-iT_m\Delta t}e_1|<\text{tol}\cdot\Delta t$. If the dimension cap of 40 is reached, the step is split in two.

**Notes on the tolerance:**
- The tolerance is an error **per unit time**, so the bound after time $t$ is $\text{tol}\cdot t$, which is $8\times10^{-10}$ at $t=8$ with tol $=10^{-10}$.
- In practice 9–10 vectors per step are enough.
- The $N+1$ sector has 56, 792, 11 440, 43 758 and 167 960 states at $L=8,12,16,18,20$.

**References.**
- `evolve_exact`: full diagonalisation of $H_{N+1}$, at $L=8,12$.
- `spectral_amplitude`: $\sum_m z_{jm}e^{-i\omega_mt}$ from Step 3's lines.
- `free_trajectory`: the $V=0$ one-body route of Section 3.

**`prepare_filtered` (item 5).**
- Projects $\Psi(0)$ onto the $N+1$ eigenstates with $0\le\omega<f\Delta_+$, then renormalises.
- $f=0.8$, unless that window holds more than one line; then the widest $f$ in $0.6$–$0.9$ holding one line. In practice $f=0.8$ everywhere except $L=18,20$, $V=3$, where $f=0.7$.
- The preparation probability $\|P_{\mathcal W}\Psi(0)\|^2$ must equal Step 3's $z_0(\mathcal W)$ at the same $f$.

**`finite_size_time(small, big)` (item 4).** Compares a chain with the next larger one, $8\to12\to16\to18\to20$, using $P_{\rm return}$, $Q_2$ and $Q_4$.
- **The two chains already differ at $t=0$** by a static offset $s$, because their ground states differ slightly near the ends.
  - If $s>\varepsilon=0.01$, the smaller chain is not size-converged at any time, and no interval is claimed. This happens for $L=8$ and $12$.
- **A return** is the first time the *change since $t=0$* differs between the two chains by more than $\varepsilon+2s$. The $2s$ allows for the relaxation of the initial difference itself; one transient at $L=16$, $V=3$, $t\approx1.85$ is correctly ignored.
- **$L=20$ has no larger partner.** Its interval is $L=18$'s, scaled by the round-trip distance $2(L-1)-(\ell_e-1)$, which assumes the same speed.

**`late_window`:** the second half of the certified interval, starting no earlier than $t=2$, and nothing if it is shorter than 1. All "late" numbers are averages over this window.

**`bound_lines`** (notebook 4.8): every $N+1$ line below $\Delta_+$ that the electron reaches, with its $z$ and its own charge $q_2$, $q_4$. $\sum z\,q_4$ is the edge charge that cannot leave.

**`far_edge_response`** (4.7): $\Delta n_{L-1}(0)$ for the selected $\Omega_N$, and for the left- and right-loaded mixtures of the $N$-sector doublet.

---

## 5. How do we know it is right?

There are 38 checks. Some are repeated for each block of sizes ([L=8], [L=8–16], [L=18, 20]). Worst deviations are from the notebook output.

| Id | What | Worst | Why it is independent |
|---|---|---|---|
| A0 | $\Omega_N$ is Step 2's state: $C=+1$, residual $<10^{-9}$, $E_0^N$ = `step2_grid.csv` | $8\times10^{-14}$ | Step 2's saved data |
| A1 | Krylov = full diagonalisation ($G$, $P$, $\Delta n_j$), $L=8,12$, edge and centre | $1.2\times10^{-12}$ | a different algorithm |
| A2 | Krylov $G(t)$ = Step 3 spectral sum, complex, phase included, $L\le16$ | $1.2\times10^{-12}$ | Step 3's lines only |
| A3 | halving the tolerance changes nothing beyond the bound | $1.3\times10^{-12}$ | convergence |
| A4 | tolerance scan: error falls $4\times10^{-6}\to2\times10^{-12}\to10^{-14}$ as tol goes $10^{-3}\to10^{-9}\to10^{-10}$, always below tol$\cdot t_{\max}$ | ratio $\le9\times10^{-4}$ | shows the tolerance actually controls the error |
| A5 | rerun bit-for-bit identical | exact | no random numbers |
| B1 | norm, $\langle H_{N+1}\rangle$, $\sum_j\Delta n_j=1$ conserved, every run up to $L=20$ | $1.5\times10^{-12}$ | exact conservation laws |
| B2 | $t=0$ identities (Section 3) | $1.9\times10^{-14}$ | ground-state correlations only |
| B3 | $V=0$: $G$, $\Delta n_j$ = one-body route | $2.9\times10^{-12}$ | $L\times L$ matrices, no Fock space |
| B4 | energy moments $m_1$, $m_2$ = closed forms | $8.8\times10^{-12}$ | derived by hand from commutators |
| C1 | filtered preparation probability = Step 3 $z_0(\mathcal W)$ at the same $f$ | $5\times10^{-15}$ | Step 3's tables |
| C2 | filtered state: $P_{\rm return}=1$, static profile, $\tfrac12$ of its charge per half | $2.2\times10^{-12}$ | it must be an eigenstate |
| C3 | $\delta=+0.30$: window empty | exact | the control |
| D1–D5 | finite-size criterion: never fires on itself; flags $L=8$ vs 12 ($s=0.080$); ignores the $V=3$ transient and finds the return at 4.90; gives $L=16$ an explicit interval; late-window logic | – | tests of the criterion itself |
| E1 | $L=8$, $t\le300$: no beat at the splitting (partner weight 0; Fourier at the splitting $\le2\%$ of the maximum) | 0.020 | exact long evolution |
| E2 | diagonal ensemble: $\tfrac12$ of the charge per half | $8\times10^{-15}$ | symmetry argument (Section 2.5) |
| E3 | $V=0$ far-edge excess $\to(1-r^2)^2/2$, monotonically; symmetry-broken states $<0.01$ at $L=16$ | 0.0076 at $L=16$ | Step 1 formula |
| E4 | $L=16$, $V=0,1,3$: late $P_{\rm return}$ within 0.05 of $z_0^2$; $\delta=+0.30$ below 0.05 | 0.036 | the physics prediction of Section 2.2 |

**The Krylov errors are about 600 times below the accumulated bound.** A1, A2 and B3 test the time evolution against three routes that share no propagation code.

---

## 6. What do the results show?

**Figures** (notebook 4.11):
- `step4_validation`;
- `step4_spacetime_density`: where the charge goes; $L=20$ maps on top, $L=16$ profiles below;
- `step4_time_traces`: $P_{\rm return}$ and charges against $t$ for edge, centre, $\delta=+0.30$ and filtered;
- `step4_finite_size_beating`;
- `step4_summary_vs_V`: late values against $V$ next to the Step 3 predictions.

Section 4.13 redraws the five figures of the short report.

### 6.1 Clean (certified) intervals, edge injection, $\delta=-0.30$

| $L$ | $V=0$ | $0.5$ | $1$ | $2$ | $3$ |
|---|---|---|---|---|---|
| 8 | $s=0.053$, none | 0.067, none | 0.080, none | 0.086, none | 0.062, none |
| 12 | $s=0.014$, none | 0.016, none | 0.020, none | 0.030, none | 0.031, none |
| 16 | all $t\le8$ | all $t\le8$ | all $t\le8$ | all $t\le8$ | $t<4.90$ |
| 18 | all $t\le16$ | – | $t<14.45$ | – | $t<6.50$ |
| 20 (scaled) | $\gtrsim18$ | – | $\approx16.3$ | – | $\approx7.3$ |

Returns come sooner as $V$ grows. The implied speed (round trip / return time at $L=18$) is $\le1.94$, $2.15$ and $4.77$ sites per unit time at $V=0,1,3$. The free band's maximum group velocity is $2\min(J_w,J_s)=1.40$. See Limits for how far to trust these speeds.

### 6.2 Coherence: the return probability

| $V$ | 0 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| $\overline{P}_{\rm late}$, edge, $L=16$ | 0.540 | 0.522 | 0.487 | 0.366 | 0.220 |
| $z_0(\mathcal W)^2$, $L=16$ | 0.505 | 0.492 | 0.462 | 0.354 | 0.209 |
| $\overline{P}_{\rm late}$, edge, $L=20$ ($t\le16$) | 0.499 | | 0.460 | | 0.228 |
| centre injection, $L=16$ | 0.098 | 0.096 | 0.072 | 0.077 | – (interval too short) |
| $\delta=+0.30$, $L=16$ | 0.020 | 0.013 | 0.010 | 0.008 | 0.008 |

- **$L=16$ sits a little above $z_0^2$.** Its averaging window ends at $t=8$, and $\sum_{\rm bulk}z^2$ is not negligible at 16 sites.
- **$L=20$ with $t\le16$ agrees to 0.006 at $V=0,1$.** This is the cleanest confirmation of Section 2.2.

### 6.3 Localisation: charge kept on the first four sites

| $V$ | 0 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| $\overline{Q_4}_{\rm late}$, $L=16$ | 0.366 | 0.356 | 0.332 | 0.275 | 0.278 |
| $\overline{Q_4}_{\rm late}$, $L=20$ ($t\le16$) | 0.327 | | 0.280 | | 0.270 |
| $z_0q_4$ (edge eigenstate), $L=16$ | 0.326 | 0.312 | 0.288 | 0.220 | 0.141 |
| $\sum z\,q_4$, all lines below $\Delta_+$, $L=16$ | 0.326 | 0.312 | 0.288 | 0.220 | 0.182 |
| $\delta=+0.30$, $L=16$ | 0.164 | 0.116 | 0.085 | 0.053 | 0.043 |

- **$L=20$ matches $z_0q_4$ to 0.008 at $V=0,1$.** $L=16$ is higher because charge is still in transit at $t\le8$.
- **$V=3$ is not understood.** About 0.27–0.28 stays at both sizes, above both predictions. The clean interval is short ($t<4.9$ and $t<7.3$), so this describes the first few time units, not a final value.
- **Step 3's second line** at $V=3$ ($\omega=1.255$, $z=0.104$) is also below $\Delta_+$. It carries $q_4=0.386$, *more* edge charge than the edge line's $0.309$. This is a new clue for Step 3's open question: the second line looks edge-bound.
- **Remember Section 2.4:** the coherent part keeps about $z_0q_4\approx z_0/2$ on the left, because the other half of it is at the far end.

### 6.4 Far edge at $t=0$

| $L$ | $V=0$ | $V=1$ | $V=3$ | left-loaded state, $V=0$ |
|---|---|---|---|---|
| 8 | 0.213 | 0.228 | 0.284 | 0.012 |
| 12 | 0.232 | 0.263 | 0.315 | 0.0017 |
| 16 | 0.244 | 0.279 | 0.332 | 0.0002 |
| 20 | 0.249 | 0.283 | 0.338 | 0.00002 |

The $V=0$ column tends to $(1-r^2)^2/2=0.252$. Repulsion increases the correlation between the ends.

### 6.5 Beating

The beat periods would be 34–40 ($L=8$), 139–144 ($L=12$), 481–658 ($L=16$), 894–1459 ($L=18$) and 1661–3283 ($L=20$). Even at $L=8$, evolved to $t=300$ (more than 7 periods), there is no peak at the splitting.

### 6.6 Filtered excitation

$z_0(\mathcal W)$ at $L=16$ is $0.710$, $0.702$, $0.680$, $0.595$, $0.457$; at $L=20$ it is $0.710$, $0.681$, $0.456$ ($V=0,1,3$). It equals Step 3 to $10^{-14}$. The filtered state has $P_{\rm return}\equiv1$ and a static profile, because it is one eigenstate. This is stable by construction, not evidence of a long lifetime. The unfiltered electron is the physical preparation.

### 6.7 How to read the main figures

- **`step4_time_traces`:**
  - top row: blue edge curves drop, then oscillate around the grey $z_0^2$ line; orange (centre) and green ($\delta=+0.30$) go to about 0;
  - bottom row: $Q_4$ settles near $z_0q_4$, and the grey dashed far-edge curve starts near 0.32 at $V=0$, not at 0;
  - shading marks times after the clean interval; ignore anything there.
- **`step4_spacetime_density`, top:** the change of charge per two-site cell at $L=20$.
  - A front leaves the edge and reflects off the far end; the dashed line is drawn at the implied speed.
  - The first cell turns blue (it loses charge) and stays blue at a fixed level: that is the part that stays.
- **`step4_spacetime_density`, bottom:** profiles at $t=0,1,2,4,8$.
  - Edge injection keeps a profile close to the dashed filtered (static) one, including the copy at sites 13–15.
  - Centre injection spreads out.
  - $\delta=+0.30$ moves into the bulk.

---

## 7. Limits

- **Sizes and times.** Only $L=16$ (to $t=8$), $L=18$ and $L=20$ (to $t=16$) have clean intervals. At $V=3$ these are short, $t<4.9$ to $7.3$. Everything is finite-size and finite-time; it is not an infinite-chain lifetime.
- **Effective speeds.** The implied return speeds (2.15 at $V=1$, 4.77 at $V=3$) exceed the free maximum of 1.40.
  - The notebook calls them upper-side effective speeds, not group velocities. A small reflected amplitude $\delta\psi$ interferes with the large static edge amplitude $A$, giving a signal of order $2A\,\delta\psi$ that is linear in $\delta\psi$. A faint precursor can therefore cross $\varepsilon$ before the main front arrives.
  - The $L=20$, $V=3$ map does show changes spreading at roughly that speed, so the intervals are a safe cut. Whether excitations really move that fast at $V=3$ is not established.
  - Centre injection at $V=3$ departs from $L=18$ already by $t\approx0.75$. That is too early to be a reflection, so it is probably a ground-state difference between sizes, and no late value is quoted there.
- **$V=3$ edge charge** stays above both predictions inside its short interval (6.3).
- **State selection.** The far-edge charge and the "half at each end" picture belong to the selected finite-$L$ ground state ($C=+1$, bonding). A symmetry-broken state would give a different, one-sided picture. Step 4 is consistent with Steps 1–3, which use the same state.
- **$L=20$** has no larger partner; its interval is scaled from $L=18$.

---

## 8. Files, reproduction, and repo status

| File | Contents |
|---|---|
| `ssh-model-project/project_A_steps1-4.ipynb` | Steps 1–3 (cells identical to `project_A.ipynb`) plus Step 4, sections 4.1–4.14 |
| `data/step4_*.csv`, `step4_trajectories.npz`, `step4_run_info.json`, `step4_checks.csv` | written by notebook 4.12 and the summary cell |
| `figures/step4_*.pdf` | written by 4.3, 4.11 and 4.13 |
| `figures/fig15_*`–`fig19_*.pdf` | committed on 17 Nov |
| `projectA_step4_short_report.pdf` | Somnath's standalone short report (in `acads/QFT/`) |

**To reproduce:** run `project_A_steps1-4.ipynb` top to bottom. Step 4 needs Steps 1–3 in memory and takes about 15 minutes; $L=18,20$ take about 5 of them, $L=20$ alone about 3. It should print `Step 4: 38/38 checks passed` and `all steps: 159/159`.

**What Somnath fixed relative to my review (`step4_review_for_somnath.md`):**
- Step 4 now reuses the notebook's Step 1–3 functions instead of a second codebase.
- It explains the far-edge charge (4.7).
- It states $P\to z_0^2$ and $Q_4\to z_0q_4$ (4.9).
- It has the $\tfrac12:\tfrac12$ diagonal ensemble (E2).
- It gives a physics conclusion (4.14).
- It states $f$ per point ($0.7$ at $L=18,20$, $V=3$).
- The $V=0$ route uses the bonding orbital (`filled_correlations`).
- Centre site $2\lfloor L/4\rfloor$ is used consistently.
- Both plain and Hann Fourier amplitudes are recorded.

**Still open in the repo:**
1. **Two notebooks.** `project_A.ipynb` (Steps 1–3) and `project_A_steps1-4.ipynb` (Steps 1–4) now coexist. The agreed layout is one notebook, so one should replace the other; probably the 1–4 notebook should become `project_A.ipynb`.
2. **`report.tex` has no Step 4 section.**
3. **Figure names don't match.** The committed PDFs are `fig15_dynamics_validation` … `fig19_summary_vs_V`, which look like they come from the earlier script version. The notebook saves `step4_validation`, `step4_spacetime_density`, `step4_time_traces`, `step4_finite_size_beating`, `step4_summary_vs_V` and `step4_s1`–`s5`, and none of those are in the repo.
4. **No Step 4 data in the repo.** `data/step4_*` is not committed. `step4_trajectories.npz` should probably go in `.gitignore`, like `step2_states.npz`.
5. **No `summary/step4.md`.** This file could be the starting point.
6. The `somnath-project` branch with the modular scripts still exists. It doesn't need to be merged.

---

## 9. Questions to be able to answer

1. Why is $G(t)=\sum_mz_{0m}e^{-i\omega_mt}$, and why must Krylov agree with it *including the phase*?
2. Why does $\overline{P_{\rm return}}\to z_0^2$? What would make it larger on a small chain?
3. Why does the edge line never decay, while the bulk part does? Answer in terms of $\Delta_+$.
4. Derive $\Delta n_j(0)=\tfrac12$ and $\Delta n_k(0)=-\langle\delta n_j\delta n_k\rangle_c/(1-\langle n_j\rangle)$.
5. Why is there charge at the far end at $t=0$? Why does it tend to $(1-r^2)^2/2$, and why does it vanish for a left-loaded state?
6. Why does the edge eigenstate carry exactly $\tfrac12$ of its charge on each half?
7. Why is the infinite-time edge charge meaningless on a finite chain?
8. Derive $[H,c_0^\dagger]=-J_0c_1^\dagger+Vc_0^\dagger(n_1-\tfrac12)$, and from it $m_2$ at $V=0$. Why does $\delta=+0.30$ start decaying more slowly?
9. How is the clean interval defined, why the allowance $2s$, and why do $L=8,12$ get none?
10. Why is there no two-edge beat, and why would you not see one at $L\ge16$ anyway?
11. Why is the filtered excitation stationary, and why is it not evidence of a long lifetime?
