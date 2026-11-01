# Step 2 summary: a many-body calculation you can trust

Read this after [`step1.md`](step1.md), or at least its first page. Each section answers one question. Sources: the guide is `QFT_student_research_guides.pdf`, Project A, p. 7; the code is Step 2 of [`project_A.ipynb`](../project_A.ipynb); the write-up is Section 3 of [`report.tex`](../report.tex).

## Short version

Step 1 solved the chain without interactions. Step 2 builds and tests the machine that solves it with interactions ($V\neq0$). Later steps depend on that machine, so the deliverable is mostly evidence that it is right.

**Verified, for $L\le18$ on this grid:**

- The matrices are the guide's Hamiltonian (A.1), including the fermion sign on a ring's closing bond.
- The eigenvalues and eigenvectors are accurate. All 61 checks pass in a fresh kernel, and all 37 of Step 1's still pass.
- The half-filled ground state is unique on the open chain, so there is an explicit rule for selecting it.
- For $\delta=-0.30$ (weak end bonds), the neutral and addition gaps fall with $L$ at every $V\le3$ and stay far below the bulk threshold. For $\delta=+0.30$ there is no low-energy pair.
- No level crossing involves the two lowest levels, and nothing special happens at $V=2$.

**Not established, and stated as hypotheses:**

- Whether an edge state survives. Energies cannot answer that; Step 3's $z_0$ and $\xi_W$ can.
- Why the neutral and addition gaps decay with different lengths at $V>0$ (Section 8).
- Whether Step 3's edge window stays isolated at $V=3$ (Section 8).

**How much to trust this.** All of this was written by an LLM, and the guide says agreement between LLM outputs is not an independent check. So the trust should rest on Section 4's checks that do not depend on this code being right, on Section 9's sabotage test, and on your own re-run (Section 12). Section 11 lists what was not tested. To present this you still need to understand it, and Sections 2 to 6 are the minimum.

| Symbol | Meaning |
|---|---|
| $N=L/2$, $M=L/2$ | half filling; number of two-site cells |
| $E_m^{N'}$ | level $m$ in the sector with $N'$ particles |
| neutral gap | $E_1^N-E_0^N$; at $V=0$, $\delta<0$ it is Step 1's end-pair splitting $\Delta E$ |
| addition gap | $E_0^{N+1}-E_0^N$, the cost of adding one particle; $\Delta E/2$ at $V=0$, $\delta<0$ |
| $\Delta_+$ | the addition gap on a ring, i.e. the bulk threshold |
| $C$, $R$ | particle-hole and reflection symmetries (Section 3) |
| $\Omega_N$ | the selected half-filled ground state |

---

## 1. What did Step 2 have to produce?

| Guide requirement (p. 7) | Where | Status |
|---|---|---|
| Integer states, fixed $N$, sparse matrices | notebook 2.1 | done |
| Sign rule, including a ring's closing bond | 2.1, 2.4 | done |
| 1. Independent matrices in $N$, $N\pm1$; diagonal interaction; hop only into an empty site | 2.1, 2.2 | done |
| 2. Whole matrix at $L=6,8,10$; sparse above; residuals $<10^{-10}$ and $<10^{-9}$ | 2.7 | worst $7\times10^{-14}$; $6\times10^{-13}$ |
| 3. Hermiticity, number conservation, eigenvalue vs independent $\langle\psi\lvert H\rvert\psi\rangle$ | 2.2, 2.7 | done |
| 4. $V=0$: energy and $\langle c_i^\dagger c_j\rangle$ vs filled orbitals | 2.3 | done |
| 5. Jordan-Wigner second implementation at $L=6$ | 2.2 | done, also $L=8$, all $N$, rings |
| Grid: $L=8$ first, then 12, 16, plus 18; refine by 0.25 where something changes | 2.9 | done, $L=18$ at every point |
| Record levels, residuals, boundary, symmetry; explicit state selection | 2.8, 2.9, `data/` | done |

## 2. Why is Step 1's $L\times L$ matrix no longer enough?

At $V=0$ the particles are independent, so one $L\times L$ matrix describes everything. The term $V(n_j-\tfrac12)(n_{j+1}-\tfrac12)$ couples particles, so we must work with whole occupation patterns.

Every hop has one $c^\dagger$ and one $c$, so $[H,\hat N]=0$, and $H$ splits into blocks of size $\binom{L}{N}$:

| $L$ | 6 | 8 | 10 | 12 | 14 | 16 | 18 |
|---|---|---|---|---|---|---|---|
| $N=L/2$ | 20 | 70 | 252 | 924 | 3,432 | 12,870 | 48,620 |
| $N\pm1$ | 15 | 56 | 210 | 792 | 3,003 | 11,440 | 43,758 |

A dense matrix at $L=18$ would need 18.9 GB. The sparse one takes about 6 MB, and solving it takes under a second.

## 3. Where do the fermion signs come from?

A state is an integer, $s=\sum_jn_j2^j$, with one fixed operator order:

$$|n\rangle=(c_0^\dagger)^{n_0}\cdots(c_{L-1}^\dagger)^{n_{L-1}}|0\rangle,\qquad c_j^\dagger|n\rangle=(-1)^{\sum_{l<j}n_l}(1-n_j)|n+e_j\rangle .$$

The sign counts the operators $c_j^\dagger$ has to pass. For a hop $c_a^\dagger c_b$, annihilate at $b$ first, then create at $a$ on the intermediate state. That gives two results:

- **Open-chain bond.** For $c_j^\dagger c_{j+1}$, the two strings count the same particles, so the sign is $+1$. Every off-diagonal element is exactly $-J_j$.
- **Ring's closing bond.** $c_{L-1}^\dagger c_0$ moves a particle past the other $N-1$, so the sign is $(-1)^{N-1}$. For odd $N$ this is $+1$, so a test in an odd-$N$ sector cannot catch a missing sign.

**Why it matters.** An open-chain energy cannot detect a missing string, because the hop signs are all $+1$ anyway. Signs enter only through $\langle c_i^\dagger c_j\rangle$ at $|i-j|\ge2$, through $c_j$ and $c_j^\dagger$ (Step 3 needs them), and through the closing bond. Those three places get their own tests.

**Symmetries**, which give free checks and state labels:

| Symmetry | Action | In the basis | Consequence |
|---|---|---|---|
| $C$ | $c_j\to(-1)^jc_j^\dagger$ | bit complement, sign $+1$ | $N\leftrightarrow L-N$: addition gap = removal gap, and a unique half-filled state has $\langle n_j\rangle=\tfrac12$, at every $V$ |
| $R$ | $c_j\to c_{L-1-j}$ | bit reversal $\times(-1)^{N(N-1)/2}$ | commutes with the open-chain $H$, because $J_{L-2-j}=J_j$ |
| ring translation | $c_j\to c_{j+1}$ | | maps $\delta\to-\delta$: ring spectra identical for $\pm\delta$ at every $V$ |

## 4. How do we know the matrix is really Eq. (A.1)?

$H=H_{\rm hop}(\delta)+VH_{\rm int}$ is linear in $V$, and $H_{\rm int}$ does not depend on $\delta$. So there are two questions.

**Is $H_{\rm hop}$ right? Test at $V=0$, against Step 1.**

- Every many-body level must be a sum of $N$ one-particle energies. This holds for the whole spectrum at $L\le10$, open and ring.
- The gaps must be Step 1's 80-digit $\Delta E$ and $\Delta E/2$. They agree to $5\times10^{-14}$.
- $\langle c_i^\dagger c_j\rangle$ must be the filled-orbital matrix built on Step 1's bonding orbital, for all $i,j$ at $L=6$ to 18. This confirms Step 1's sign $s=(-1)^{L/2}$ for both parities of $M$.

**Is $H_{\rm int}$ right? Test at $V\neq0$.** These checks are the important ones, because the $V=0$ checks cannot see $H_{\rm int}$ at all.

| Route | Reference | Worst deviation |
|---|---|---|
| Fermion algebra | $\{c_i,c_j^\dagger\}=\delta_{ij}$ holds exactly, so $H$ written as a literal product of $c$, $c^\dagger$ is correct by construction; the matrix must equal it | exact |
| Jordan-Wigner, Eq. (A.4) | a separate build on all $2^L$ states | exact (blocks) |
| By hand | $L=2$ and $L=4$, whole spectrum (below) | $7\times10^{-15}$ |
| Weak coupling | $dE/dV$ at $V=0$ from Wick's theorem and Step 1 orbitals | $4\times10^{-8}$ |
| Strong coupling | second-order perturbation theory (below) | scales as predicted |
| Bethe ansatz, $\delta=0$ | a different algorithm; solves equations, builds no matrix | $7\times10^{-14}$ |
| Literature, $\delta=0$, $V=2$ | energy per site $\to2(\tfrac14-\ln2)$ (Hulthén) | $10^{-3}$, finite size |

**You can check the $L=2$ case yourself.** With one particle, the two states $|10\rangle,|01\rangle$ have $(n_0-\tfrac12)(n_1-\tfrac12)=-\tfrac14$ and are joined by $-J_0$, so $E=-V/4\mp J_0$. Empty or full gives $+V/4$.

**$L=4$, $N=2$.** Write $a=J_0=J_2$, $b=J_1$. Label states by their occupied sites. The interaction energies are $\{0,1\},\{2,3\}$: $+V/4$; $\{0,2\},\{1,3\}$: $-3V/4$; $\{0,3\},\{1,2\}$: $-V/4$. Taking combinations with definite bit reversal and bit complement gives

$$\begin{pmatrix}V/4&-b&0\\-b&-3V/4&-2a\\0&-2a&-V/4\end{pmatrix},\qquad\begin{pmatrix}V/4&-b\\-b&-3V/4\end{pmatrix},\qquad(-V/4).$$

On the ring the closing bond has sign $(-1)^{N-1}=-1$. There $C$ alone already reduces the problem, to

$$\{0,\ 0,\ -\tfrac V2\pm\sqrt{\tfrac{V^2}{4}+4a^2},\ -\tfrac V2\pm\sqrt{\tfrac{V^2}{4}+4b^2}\}.$$

Dropping the sign changes this spectrum by about 1, so this one line tests the sign at $V\neq0$.

**Strong coupling.** At large $V$ the ground states are the two patterns $1010\ldots$ and $0101\ldots$, with every bond anti-aligned, at energy $-(L-1)V/4$. A hop on an interior bond spoils two neighbouring bonds and costs $V$. A hop on an end bond spoils one and costs $V/2$. Second-order perturbation theory gives

$$E_0^N=-\frac{(L-1)V}{4}-\frac1V\Big[\sum_{j=1}^{L-3}J_j^2+2(J_0^2+J_{L-2}^2)\Big]+O(V^{-3}).$$

An added particle creates one aligned bond $11$ (cost $V/2$). It moves two sites at a time through the odd bonds, with amplitude $1-\delta$, on $M$ positions, so

$$E_0^{N+1}-E_0^N=\frac V2-2(1-\delta)\cos\frac{\pi}{M+1}+O(V^{-1}).$$

The errors at $V=200$ and $400$ differ by the predicted factors of 8 and 2, which is stronger evidence than a small error alone.

**Bethe ansatz.** At $\delta=0$, a sublattice rotation gives $H=2H_{\rm XXZ}(\Delta=V/2)$, which is integrable. The notebook solves the Bethe equations and gets the same energies at $V=0.5,1,1.5,2$, up to $L=18$, on open chains and odd-$N$ rings. I wrote the equations from memory of the literature and derived the open-end term from the one-magnon problem. I did not read Alcaraz et al. (1987) in full. They are validated by the $L=2,4$ closed forms and the agreement itself: a wrong equation would not match to 14 digits.

## 5. How do we know the eigenpairs are right?

| Test | Worst deviation | Target |
|---|---|---|
| Residuals, all eigenpairs, $L\le10$ | $7.0\times10^{-14}$ | $10^{-10}$ |
| Residuals, production, open / ring | $1.3\times10^{-13}$ / $6.1\times10^{-13}$ | $10^{-9}$ |
| Lanczos vs dense, $L=12$ | $7.1\times10^{-14}$ | $10^{-11}$ |
| Lanczos finds exact multiplicities, $V=0$, $L\le18$ | $6.2\times10^{-14}$ | $10^{-11}$ |
| Eigenvalue vs $\langle H\rangle$ from correlators (never calls `hamiltonian()`) | $7.8\times10^{-14}$ | $10^{-9}$ |

**A failure worth knowing about.** My plan set a $10^{-12}$ tolerance on $\langle c_i^\dagger c_j\rangle$ and tested it only to $L=12$. Extending it to $L=18$ made it fail, at $2\times10^{-12}$, only for $\delta<0$.

The cause is that a residual $\|r\|$ only bounds the eigenvector error by $\|r\|/\text{gap}$ (Davis-Kahan), and an observable can shift by up to $2\|r\|/\text{gap}$. With $\delta<0$ the nearest level is the end partner, with gap $7\times10^{-3}$ at $L=18$. So a near-perfect residual still leaves a small admixture of the partner.

I did not loosen the tolerance. I tested the explanation, which makes a prediction: the partner has $C=-1$, so projecting onto $C=+1$ must remove the error. At $L=18$, in the saved run:

| | $\max\lvert C-C_{\rm filled}\rvert$ |
|---|---|
| raw Lanczos vector | $1.2\times10^{-12}$ |
| bound $2\|r\|/\text{gap}$ | $4.8\times10^{-12}$ |
| after $C$ projection | $2.0\times10^{-15}$ |

So the $10^{-12}$ tolerance stands, for the projected vector, and the raw vector is checked against the bound. **For later steps:** the guide's $10^{-9}$ residual target, with a gap of $4\times10^{-3}$, would allow an eigenvector error of $2.5\times10^{-7}$. Our residuals are about $10^{-13}$.

## 6. Which half-filled state do we use?

This is the guide's "before moving on" requirement, and it generalises Step 1's rule of filling the bonding orbital.

On the open chain, every off-diagonal element is $-J_j<0$ (Section 3), and any configuration can reach any other by hops. The **Perron-Frobenius theorem** then says the lowest level is non-degenerate, with all amplitudes positive, at every $V$ and $\delta$. $C$ and bit reversal send positive vectors to positive vectors, so

$$C=+1,\qquad R=(-1)^{N(N-1)/2},\qquad\langle n_j\rangle=\tfrac12 .$$

**Rule.** Take the lowest eigenvector with $\sum_s\psi_s>0$ and verify the above. If the neutral gap is below $\eta=10^{-8}$, the eigenvector order means nothing, so diagonalise $C$ inside the cluster and keep the $C=+1$ vector. If $C$ does not fix the state, raise an error.

On this grid the second branch never fires, since the smallest neutral gap is $4.3\times10^{-3}$. It is there for Step 7's long chains, so it is tested on a constructed rotated pair, which it maps back to the right state.

**Why bother?** An equal mixture of the end pair, which a solver could return once the pair is nearly degenerate, moves $\langle n_0\rangle$ from $\tfrac12$ by $0.38$ at $L=12$, $V=1$.

On a ring with even $N$ the closing sign is $-1$, so Perron-Frobenius does not apply. Rings are used here for energies only.

## 7. What does the grid show?

Look at [`figures/step2_gaps.pdf`](../figures/step2_gaps.pdf). Neutral and addition gaps (neutral / addition):

| $\delta$ | $V$ | $L=8$ | $L=12$ | $L=16$ | $L=18$ | ring $\Delta_+$, $L=16$ |
|---|---|---|---|---|---|---|
| $-0.30$ | 0 | 0.158 / 0.0791 | 0.0451 / 0.0226 | 0.0131 / 0.00653 | 0.00703 / 0.00351 | 0.6000 |
| | 0.5 | 0.172 / 0.100 | 0.0442 / 0.0241 | 0.0112 / 0.00586 | 0.00566 / 0.00291 | 0.7585 |
| | 1 | 0.184 / 0.131 | 0.0446 / 0.0288 | 0.0103 / 0.00605 | 0.00493 / 0.00279 | 0.9106 |
| | 2 | 0.194 / 0.238 | 0.0464 / 0.0608 | 0.00980 / 0.0142 | 0.00441 / 0.00675 | 1.1969 |
| | 3 | 0.177 / 0.423 | 0.0437 / 0.162 | 0.00954 / 0.0644 | 0.00431 / 0.0412 | 1.4806 |
| $+0.30$ | 0 | 1.57 / 0.783 | 1.41 / 0.706 | 1.34 / 0.670 | 1.32 / 0.658 | 0.6000 |
| | 0.5 | 1.68 / 0.960 | 1.54 / 0.873 | 1.48 / 0.832 | 1.47 / 0.820 | 0.7585 |
| | 1 | 1.70 / 1.13 | 1.57 / 1.03 | 1.52 / 0.989 | 1.50 / 0.975 | 0.9106 |
| | 2 | 1.57 / 1.46 | 1.42 / 1.34 | 1.35 / 1.28 | 1.33 / 1.26 | 1.1969 |
| | 3 | 1.34 / 1.77 | 1.15 / 1.61 | 1.05 / 1.54 | 1.02 / 1.52 | 1.4806 |

Reading it:

- At $V=0$ every curve starts on the Step 1 value (open circles).
- $\delta=-0.30$: both gaps fall with $L$ at every $V$. At $L=16$ the addition gap stays below 5 % of $\Delta_+$.
- $\delta=+0.30$: the neutral gap stays above 1, and at $L=16$ the addition gap lies 0.06 to 0.08 above $\Delta_+$. So reversing $\delta$ still removes the low-energy pair.
- Ground state and level 1 have fixed labels. The ground state has $C=+1$, $R=(-1)^{N(N-1)/2}$; level 1 has $C=-1$ and the opposite $R$, at every point. No crossing involves them.
- Levels 2 to 5 do cross, but stay at least 0.7 above the ground state.

**Refinement.** The guide says to refine "if a change occurs". I made that precise: refine an interval if (a) the labels of level 0 or 1 change, (b) a gap stops being monotonic in $V$ or changes its trend with $L$, or (c) neighbouring ground states overlap below 0.99. That reading is mine.

| Trigger | Coarse grid | 0.25 grid |
|---|---|---|
| (a) labels | never | never |
| (b) not monotonic in $V$ | many intervals | only smooth extrema (below) |
| (b) trend with $L$ | never | never |
| (c) overlap $<0.99$ | $1\le V\le3$ | never (minimum 0.9964) |

Every interval was refined except $0\le V\le0.5$ at $\delta=+0.30$. The remaining extrema: the $\delta=+0.30$ neutral gap peaks near $V=1$ at every $L$; the $\delta=-0.30$ addition gap has a minimum near $V=0.75$ at $L=16,18$; the $\delta=-0.30$ neutral gap has a broad maximum or plateau around $V=2$ at $L=8,12,16$, where at $L=16$ it changes by only about $10^{-5}$ per step of 0.25. These are smooth, well above the numerical noise, and show no sign of a transition. **Nothing singles out $V=2$**, which matches the guide's warning that $V=2$ belongs to the uniform chain.

## 8. What does it suggest for Step 3? (Hypotheses, not results.)

**The two gaps stop sharing one length.** At $V=0$ they are locked: addition = neutral/2. At $V>0$ they separate. Decay lengths at $\delta=-0.30$, from sizes with the same parity of $M$, so that an even/odd effect cannot fake this:

| Gap | sizes | $V=0$ | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|---|
| neutral | 12, 16 | 3.23 | 2.92 | 2.73 | 2.57 | 2.63 |
| | 14, 18 | 3.23 | 2.92 | 2.72 | 2.52 | 2.55 |
| addition | 12, 16 | 3.23 | 2.83 | 2.57 | 2.75 | 4.34 |
| | 14, 18 | 3.23 | 2.84 | 2.57 | 2.71 | 4.43 |
| $U_{\rm eff}$ | 12, 16 | | 1.87 | 2.03 | 2.86 | 4.68 |
| | 14, 18 | | 1.82 | 1.98 | 2.83 | 4.71 |

**A possible reading.** The most general number-conserving Hamiltonian of two end modes with symmetries $C$ and $R$ is $t(a_L^\dagger a_R+\text{h.c.})+U(n_L-\tfrac12)(n_R-\tfrac12)$. With one particle the energies are $-U/4\mp|t|$; empty or full gives $+U/4$. So

$$\text{neutral gap}=2|t|,\qquad\text{addition gap}=|t|+\frac U2,\qquad U_{\rm eff}\equiv2\times\text{addition}-\text{neutral}.$$

This suggests that tunnelling between the ends gets shorter-ranged while an effective repulsion between them gets longer-ranged and comes to dominate the addition gap at large $L$. It would also explain why the addition-gap length is still rising (4.17, 4.34, 4.47 for sizes (8,12), (12,16), (16,18)) while $U_{\rm eff}$'s moves less (4.59, 4.68, 4.72). The four low levels are well isolated: at $L=18$, $V=3$, $U_{\rm eff}=0.078$ and the next levels lie about 1.15 higher. It remains a hypothesis: nothing here shows that this two-mode model is the right description.

**The edge window may not stay isolated.** Step 3 uses $0\le\omega<0.8\Delta_+$. The second addition level at $V=3$, $\delta=-0.30$:

| $L$ | second level $\omega$ | $0.8\Delta_+$ (same-$L$ ring) | inside the window? | inside at factor 0.9? |
|---|---|---|---|---|
| 12 | 1.502 | 1.253 | no | no |
| 16 | 1.255 | 1.185 | no | yes |
| 18 | 1.181 | 1.245 | **yes** | yes |

$L=18$ has an odd-$M$ ring, whose $\Delta_+$ overshoots, but the level is also just below the $L=16$ ring's 1.185. At $V\le2$ the second level stays outside even the 0.9 window. Step 3 has to test this directly and, as the guide says, report it if no isolated window exists.

## 9. Can the checks actually fail?

A test suite that cannot fail proves nothing. Three negative controls are in the notebook:

| Removed | Result |
|---|---|
| ring's closing sign | fails for every even $N$ (by at least 0.061), for no odd $N$ |
| sign string in $c_j$ | adjacent correlators unchanged; $\lvert i-j\rvert\ge2$ off by up to 0.68 |
| state selection | $\langle n_0\rangle$ off by 0.38 |

I also broke the code on purpose, one error at a time, running the checks up to the $L=8$ gate each time: 16 errors, of which 15 were caught. [`checks/step2_mutation_test.py`](../checks/step2_mutation_test.py) reproduces this. The unmodified code fails 0 checks.

| Error introduced | Checks failing | Examples |
|---|---|---|
| interaction twice too strong | 14 | G1, G3, G4, D1 |
| $n_jn_{j+1}$ instead of $(n_j-\tfrac12)(n_{j+1}-\tfrac12)$ | 19 | D1, G1, G4 |
| ring interaction misses the closing bond | 6 | D1, G1c, G4b |
| ring hop without its fermion sign | 8 | C1, E1, G1c |
| $\delta\to-\delta$ in the matrix | 16 | C1, C2, G1 |
| second sign on the original state | 26, and the rule raises | A1, and most others |
| string includes site $j$ itself | 12 | A5, C2, F0 |
| $c_j$ without the string | 6 | A5, C2 |
| $R$ without its sign | 3 | F0, F4, F5c |
| $C$ with sign $-1$ | 3 | C2, F5b, S1 |
| Lanczos targets smallest $\lvert E\rvert$ | solver raises | no convergence; with loose settings at $L=12$ it returns $-0.02$ against $E_0=-8.59$, which B3 catches |
| Lanczos stopped early | 10 | B3, C1b, C3 |
| Jordan-Wigner ring bond without its string | 2 | D1, D2 |
| selection keeps the $C=-1$ vector | 1 | S1 |
| selection accepts the first vector | 2 | S1, S3 |
| Jordan-Wigner bit order reversed | **0** | harmless, see below |

The one uncaught change is not a bug. Reversing the bit order means relabelling $j\to L-1-j$, and $H$ is reflection symmetric, so the matrix stays the same to $2\times10^{-16}$. With asymmetric bonds the same change alters the matrix by 0.2. So my plan's warning about this "trap" does not apply to this model.

## 10. What changed relative to the plan (`step2.md`)?

- **Added $V\neq0$ anchors.** Fermion algebra, literal-product comparison, $L=2,4$ closed forms, first-order Wick, strong coupling, Bethe ansatz, Hulthén. The plan's $V\neq0$ evidence was mainly Jordan-Wigner, written by the same author with the same bonds, plus symmetries.
- **Added** the sabotage test, a test for Lanczos multiplicities, and direct tests that `apply_C`, `apply_R` implement their operator maps.
- **Extended** the correlator check to $L=18$ (Section 5), and added $L=10,14$ so sizes of the same $M$ parity can be compared. Rings run at every point at $L=8$ to 18.
- **Corrected.** $L=12$ now uses Lanczos, as the guide asks (the plan used dense up to dimension 1000). Trigger (c) does fire on the coarse grid (the plan quoted only the 0.25 grid). The bit-order trap does not apply. The plan's `STEP2_FIRST = len(CHECKS)` no longer works, because Kavish made `CHECKS` a dict, so Step 2 checks are named `S2 ...` and counted by prefix.
- **Files.** `step2_grid.csv` has one row per point and `step2_levels.csv` one per level, instead of one combined table. There is no `runtime_s` column. `step2_states.npz` keeps vectors on the coarse $V$ grid only (4.5 MB), since refined points recompute in about a second.
- Energies match the plan's reference values, which is a sanity check only and not independent evidence.

## 11. Limits and what was not tested

- **Report build.** `report.tex` builds with pdfLaTeX (TeX Live 2026) with no errors and no overfull boxes, and I looked at the rendered Step 2 pages. A plain build on this machine stops in the preamble because `texlive-ec` is not installed (Section 12). It has not been tried on Overleaf.
- **Report length.** It is already 10 pages with Steps 3 to 8 still to come, and the guide asks for about 8 to 12. Sections 2 and 3 will need trimming, or their check tables moving to an appendix.
- **Exact references at $V\neq0$ cover limits.** Closed forms stop at $L=4$, perturbation theory holds at small and large $V$, and the Bethe ansatz needs $\delta=0$. At $\delta=\pm0.30$, $V\sim1$, $L=18$ no exact number exists. There the argument is: $H$ is linear in $V$, with each part tested separately, plus the symmetries, residuals and an independent $\langle H\rangle$.
- **Jordan-Wigner shares `bonds()`** with the main code. Step 1's checks and the $L=2,4$ closed forms cover that function.
- **The guide's benchmark script was not run.** `checks/benchmark_checks.py` from the guide's package is not in this repo. My `checks/` folder contains only the sabotage test.
- **Size limit.** $L\le18$, so $L/\xi\le5.6$ at $\delta=-0.3$. Two-size lengths are estimates, and they still drift.
- **Refinement triggers are my reading** of the guide's "if a change occurs".
- **Mutation test** stops at the $L=8$ gate, so it does not exercise the grid bookkeeping or the figure. Its list of errors is my choice, so a pass shows the checks catch those errors, not every possible one.
- **Round-off values change slightly between runs**, depending on thread count. Quoted deviations are from the saved run in `data/step2_checks.csv`.
- **Environment.** The notebook asks for a `qiskit-env` kernel, which is not on this machine, so I ran it with the project `.venv` (Python 3.12.14). Step 1's committed outputs are from Kavish's run with slightly older packages. A fresh run here reproduces all of their printed text exactly.
- `\author{[your name]}` in `report.tex` is still a placeholder.

## 12. Files and how to reproduce

| File | Contents |
|---|---|
| [`project_A.ipynb`](../project_A.ipynb) | Step 2 appended; Step 0 and 1 cells unchanged |
| [`report.tex`](../report.tex) | Section 3 added, two references added |
| [`figures/step2_gaps.pdf`](../figures/step2_gaps.pdf) | gaps against $V$ |
| [`data/step2_grid.csv`](../data/step2_grid.csv) | one row per $(L,\delta,V,\text{boundary})$: gaps, labels, state choice, overlap, $\Delta_+$ |
| [`data/step2_levels.csv`](../data/step2_levels.csv) | one row per level: energy, residual, labels, solver settings |
| [`data/step2_checks.csv`](../data/step2_checks.csv) | every check with its worst deviation |
| [`data/step2_run_info.json`](../data/step2_run_info.json) | conventions, tolerances, seed, versions, command |
| [`checks/step2_mutation_test.py`](../checks/step2_mutation_test.py) | the sabotage test |
| [`requirements.txt`](../requirements.txt) | package versions used |

`step2_states.npz` holds $\Omega_N$ on the coarse $V$ grid, with keys like `open_L16_d-0.30_V2.00`. It is not meant to be committed: it is 4.5 MB, nothing reads it, and the notebook writes a new copy to `data/` on every run. The copy from the saved run is outside the repo, in `../step2_extras/`.

Everything is uncommitted, as you asked. To reproduce, from the repo root, with the project `.venv`. I ran exactly this on a copy: all energies and labels came out identical to the 15 digits saved, and the state vectors matched to $2\times10^{-15}$.

```bash
jupyter nbconvert --to notebook --execute --inplace project_A.ipynb   # 3 minutes on an idle machine; prints 37/37 and 61/61
python checks/step2_mutation_test.py                                  # sabotage test; up to about 20 minutes
```

Running the notebook in place rewrites Step 1's outputs and figure PDFs too, with the same content but new metadata, so run it on a copy if you want a clean `git diff`.

To build the report here, keeping build files out of the repo:

```bash
latexmk -pdf -jobname=report -interaction=nonstopmode -usepretex='\RequirePackage{lmodern}' -outdir=../step2_extras/report_build report.tex
```

The `-usepretex` part is needed only because the preamble loads `fontenc` before `lmodern`, which briefly asks for the font `ecrm1095` from `texlive-ec`. That package isn't installed here, and the PDF doesn't use the font.

**If you want to check this yourself without trusting me, in order of value:**

1. Derive the $L=2$ result in Section 4 by hand, then the $L=4$ blocks. That takes a few minutes and tests $V\neq0$ directly.
2. Run the notebook and read the `[pass]` lines.
3. Break something in `hamiltonian()` and watch the checks fail.
4. Read the sign derivation in Section 3, and notebook cells 2.1 and 2.4.

## 13. Record of LLM assistance

The guide asks for a short record of substantive assistance and how it was verified.

- **Tool.** Claude Opus 5.5 (Anthropic), through Claude Code, on 2026-11-01. An earlier session of the same model wrote `step2.md`, the plan.
- **What it did.** All of Step 2: re-derived the conventions and limits, wrote the notebook cells, the sabotage test, the report section and this summary. No part of Step 2 was written independently by a person.
- **How it was verified.** The checks in Sections 4, 5 and 9. Agreement with the plan's reference values does not count, because both came from the same model.
- **What a person still needs to do.** Read and understand the derivations and code, and decide whether the two hypotheses in Section 8 belong in the final conclusions.

## Appendix: all 61 checks

`data/step2_checks.csv` lists each one by name. "Exact" means an integer or yes/no comparison, or a deviation of exactly 0.

| Id | What | Worst deviation | Tolerance |
|---|---|---|---|
| A1 | $H$ symmetric, $N$ and $N\pm1$, open and ring | exact | exact |
| A2 | wrong-number state rejected at construction | exact | exact |
| A3 | Jordan-Wigner $H$ conserves $N$ on all $2^L$ states | exact | exact |
| A4 | sector sizes $\binom{L}{N}$, $L\le18$ | exact | exact |
| A5 | $\{c_i,c_j^\dagger\}=\delta_{ij}$, $\{c_i,c_j\}=0$, $L=6,8$, all $N$ | exact | exact |
| A6 | `hamiltonian()` = literal product of $c,c^\dagger$ | exact | $10^{-13}$ |
| D1 | Jordan-Wigner blocks = occupation matrices, all $N$, open and ring | exact | $10^{-13}$ |
| D2 | Jordan-Wigner full spectrum = union of sectors | $4.9\times10^{-14}$ | $10^{-12}$ |
| C1 | $V=0$ whole spectrum = sums of $N$ one-particle energies, $L\le10$ | $3.6\times10^{-14}$ | $10^{-12}$ |
| C1b | $V=0$ Lanczos levels with exact multiplicities, $L=12$ to 18 | $6.2\times10^{-14}$ | $10^{-11}$ |
| C2 | $V=0$ $\langle c_i^\dagger c_j\rangle$ = filled orbitals, $L\le18$, $C$-projected | $2.0\times10^{-15}$ | $10^{-12}$ |
| C2a | same, raw vector / bound $(2\|r\|/\text{gap}+10^{-13})$ | ratio 0.23 | $<1$ |
| C2n | negative control: no string breaks $\lvert i-j\rvert\ge2$ only | 0.68 | $>10^{-2}$ |
| C3 | $V=0$ gaps = 80-digit $\Delta E$, $\Delta E/2$ | $4.8\times10^{-14}$ | $10^{-10}$ |
| C4 | $V=0$ addition weights = lowest empty orbital | $6.4\times10^{-14}$ | $10^{-12}$ |
| C5 | sum rule $\sum_mw_{jm}=1-\langle n_j\rangle$, every $V$ | $2.2\times10^{-15}$ | $10^{-12}$ |
| E1 | ring $V=0$, $E_0$ = free fermions, every $N$ | $4.4\times10^{-15}$ | $10^{-12}$ |
| E2 | negative control: no closing sign fails for even $N$ only | 0.061 | $>10^{-3}$ |
| E3, E4 | ring $V=0$, $\Delta_+$ for even and odd $M$ | $2.9\times10^{-14}$, $3.7\times10^{-14}$ | $10^{-11}$ |
| F0 | operator maps of $C$ and $R$ | exact | $10^{-13}$ |
| F1 | spectrum$(N+1)$ = spectrum$(N-1)$ | $3.6\times10^{-14}$ | $10^{-12}$ |
| F2 | $[C,H]=[R,H]=0$ | $3.6\times10^{-15}$ | $10^{-12}$ |
| F3, F3b | ring spectra and $\Delta_+$ for $\pm\delta$ | $6.3\times10^{-14}$, $5.9\times10^{-14}$ | $10^{-12}$, $10^{-10}$ |
| F4 | $R$ of $\Omega_N$ = product of orbital parities | $4.9\times10^{-15}$ | $10^{-12}$ |
| F5a–d | positivity, $C$, $R$ (also $N\pm1$), $\langle n_j\rangle=\tfrac12$; at $L=8$ and on the whole grid (8 checks) | $-8.9\times10^{-16}$; exact; exact; $4.2\times10^{-12}$ | $-10^{-12}$; exact; exact; $10^{-9}$ |
| F6 | every kept level has definite labels | exact | exact |
| G1a–c | $L=2$; $L=4$ open; $L=4$ ring | $2.2\times10^{-16}$, $6.2\times10^{-15}$, $7.1\times10^{-15}$ | $10^{-14}$, $10^{-13}$, $10^{-13}$ |
| G2 | weak coupling vs Wick, three energies | $4.1\times10^{-8}$ | $10^{-6}$ |
| G3a, b | strong coupling at $V=400$, with error scaling | $3.0\times10^{-6}$, $6.4\times10^{-3}$ | $10^{-5}$, $10^{-2}$ |
| G4a, b | Bethe ansatz, open; ring | $6.2\times10^{-14}$, $6.8\times10^{-14}$ | $10^{-10}$ |
| G5 | Hulthén energy per site | $1.0\times10^{-3}$ | $2\times10^{-3}$ |
| B1 | residuals, all eigenpairs, $L\le10$ | $7.0\times10^{-14}$ | $10^{-10}$ |
| B2 | production residuals: $L=8$; whole grid; rings (3 checks) | $1.4\times10^{-14}$, $1.3\times10^{-13}$, $6.1\times10^{-13}$ | $10^{-9}$ |
| B3, B3b | Lanczos vs dense; vs another seed | $7.1\times10^{-14}$, $8.7\times10^{-14}$ | $10^{-11}$ |
| B4a, b | eigenvalue vs independent $\langle H\rangle$: $L\le10$; $L=8$ and whole grid | $2.3\times10^{-14}$; $7.8\times10^{-14}$ | $10^{-12}$; $10^{-9}$ |
| S1 | rotated pair mapped back to the selected state | $1.2\times10^{-14}$ | $10^{-12}$ |
| S2 | negative control: unselected mixture shifts $\langle n_0\rangle$ | 0.38 | $>0.1$ |
| S3 | rule raises when $C$ does not fix the state | exact | exact |
| R1 | labels of levels 0, 1 never change | exact | exact |
| R2 | neighbouring overlap on the refined grid | minimum 0.9964 | $>0.99$ |
| R3 | cluster branch unused; smallest neutral gap | $4.3\times10^{-3}$ | $>10^{-8}$ |
| H1 | two-size lengths reproduce $\xi$ at $V=0$ | $5.5\times10^{-3}$ | $10^{-2}$ |
| H2 | $U_{\rm eff}=0$ at $V=0$ | $7.5\times10^{-14}$ | $10^{-11}$ |
