# Step 6 summary: the bulk and the competing ordered state

Read this after [`step5.md`](step5.md). Each section answers one question.

**Sources:**
- the guide, `QFT_student_research_guides.pdf`, Project A, p. 11;
- the code, Step 6 of [`project_A.ipynb`](../project_A.ipynb), sections 6.1-6.11;
- the write-up, the Step 6 section of [`report.tex`](../report.tex).

## Short version

Steps 3-5 found a localised, partly coherent excitation at the weak-bond end. But a density bump at a boundary does not prove a topological phase: a charge-density wave (CDW) or a plain weak end bond can also put charge there. Step 6 goes back to the bulk of the same Hamiltonian, on a ring, and asks: **is it an insulator with the quantised polarisation that goes with the edge, and where does the competing CDW take over?**

**Done and verified:**
- **Grid.** Rings $L=10,14,18$ (odd $M$), $\delta=\pm0.30$, $V=0$ to $10$ in steps of $0.25$. Free fermions up to $L=202$. Open chains $L=12,16,20$ up to $V=10$.
- **Perturbations.** At $V=1$ and $3$ (plus $V=8$ for contrast): five seeds of reflection-paired bond disorder of amplitude $0.05$, and a potential $0.1\,n_0$. Ring and open chain.
- **20 checks** (Section 5), all passed.

**Main results:**
- **$V=0$:** $P=\tfrac12$ for $\delta=-0.30$ (the weak-end chain) and $P=0$ for $\delta=+0.30$. $|z|=0.62,\,0.72,\,0.78$ at $L=10,14,18$, and $1-|z|\simeq4.4/L$. The Zak phase gives the same $P$: the electrons' Wannier centre sits half a cell away from the background's.
- **$V\le3$** (the range of Steps 3-5): gaps open and nearly size-independent, $|z|$ growing with $L$, no charge order ($L\langle M^2\rangle$ saturates), and $P$ stays exactly $\tfrac12$ and $0$.
- **The CDW takes over between $V\approx4$ and $5$.** $L\Delta_0$ curves cross at $V=4.8$-$4.9$. Beyond that, $\langle M^2\rangle\to0.22$-$0.23$ independent of $L$, the partner splitting closes exponentially, $|z|$ collapses (below 0.1 from $V\approx7.3$), and $P$ carries no information.
- **Symmetry-preserving bond disorder** leaves $P$ exactly quantised, with $\Delta P=\tfrac12$, even though translation is broken. **The end potential $0.1\,n_0$** breaks $C$ and $R$: $P$ shifts by about $0.09/L$ at $V=1$, and at $V=8$ it pins one CDW pattern.
- **The edge line** (open $L=12$) survives the bond disorder: it stays near $\omega=0$ with $z_0\approx0.63$-$0.71$ at $V=1$. The end potential moves it up (0.029 to 0.088), while the bulk $P$ is untouched.
- **Open chains at strong $V$:** the weak-end chain's near-degenerate pair is an edge pair at $V\le3$ and turns into the CDW doublet from $V\approx5$. The other dimerisation gets its own CDW doublet at $V\gtrsim6$. At $V=3$ the edge charge already induces staggering in the middle of the chain (0.13-0.23), but it decreases with $L$.

| Symbol | Meaning |
|---|---|
| $x_j$ | $n_j-\tfrac12$ |
| $X$ | $\sum_j j\,x_j$, the dipole relative to the background, origin at site 0 |
| $z$, $P$ | $\langle e^{2\pi iX/L}\rangle$; $P=\arg z/2\pi$ mod 1, folded into $[-\tfrac14,\tfrac34)$ |
| $M_{\rm CDW}$ | $\frac1L\sum_j(-1)^jx_j$ |
| $\Delta_0$, $\Delta_\pm$ | neutral gap $E_1^N-E_0^N$; addition/removal gap $E_0^{N\pm1}-E_0^N$ |
| $C$, $R$, $T$ | particle-hole $c_j\to(-1)^jc_j^\dagger$; reflection $j\to L-1-j$; translation $j\to j+1$ |
| partner | the first excited half-filled level |
| $z_0(\mathcal W)$, $\omega_{\rm edge}$ | Step 3's window weight at site 0 and the lowest addition line |
| $Q_4$ | $\sum_{j<4}(\langle n_j\rangle-\tfrac12)$ of a $C$-breaking combination |

---

## 1. What did Step 6 have to produce?

| Guide item (p. 11) | Where | Status |
|---|---|---|
| Ring with the fermionic sign on the closing bond; $L=10,14,18$ | 6.1, A1, A3 | done; A3 shows why these sizes |
| 1. $V=0$: $\Delta P=\tfrac12$ mod 1; Berry phase if desired; plot $\lvert z\rvert$ with $P$ | 6.3, `step6_polarisation` | done, with the Slater determinant, the Berry product and the Zak phase |
| 2. Finite $V$: neutral gap, addition/removal gaps, $P$; ground-state overlap | 6.4, `step6_bulk_vs_V` | $V=0$-$10$ in steps of 0.25 |
| 3. $\langle M^2_{\rm CDW}\rangle$, its size dependence and the near-degenerate partner | 6.5, `step6_charge_order` | done, with the partner's symmetry and matrix element |
| 4. Paired bond perturbations of 0.05, five seeds, two points; $0.1\,n_0$; which symmetry each keeps | 6.6, 6.7, `step6_perturbations` | $V=1,3$ (plus 8); ring and edge line; symmetries computed (D2, D4) |
| Before moving on: compare explicitly with the charge-ordered regime | 6.8, `step6_open_doublets` | open chains $L=12$-$20$ to $V=10$ |

## 2. What is being measured, and why?

**Why the bulk.** The edge state of Step 1 exists because the two bond patterns are different phases of the bulk. An end that cuts the chain where the invariant changes must carry a zero mode. With interactions there are no bands, so the bulk invariant has to be measured on the many-body ground state. On top of that, at strong $V$ there is a different, symmetry-broken phase, the CDW. The task is to tell which of the two the bulk is in.

**Why a ring with odd $M$.**
- In operator language the ring is periodic, $c_L\equiv c_0$. In the occupation basis the closing hop passes the other $N-1$ particles, so its matrix element carries $(-1)^{N-1}$.
- At $\delta=V=0$ the energies $-2\cos k$ vanish at $k=\pm\pi/2$, and these momenta are on the grid $2\pi m/L$ only if $4\mid L$. In that case half filling leaves one particle in two zero modes, and the ground state is degenerate. For $L=10,14,18$ the negative levels fill exactly (check A3).
- For odd $N$ every off-diagonal matrix element is $-J_j<0$. So Perron-Frobenius applies: the ground state is unique and positive at every $V$, and no selection rule is needed (A4).

**Why $z$ and not $\langle X\rangle$.** On a ring, a hop across the closing bond changes $X$ by $L$, so $X$ is not an observable. $e^{2\pi iX/L}$ is periodic in every coordinate and therefore well defined. It is diagonal in the occupation basis: each configuration is a point on the unit circle at angle $2\pi X/L$, and $z$ is their weighted mean.
- $\arg z$ locates the charge centre.
- $|z|$ says how sharply that centre is defined: close to 1 in an insulator, close to 0 when the configurations are spread or sit at opposite points.

**What makes $P$ quantised.** $C$ and $R$ both send $X\to-X$ at half filling. For a ground state that is an eigenstate of either, $z=z^*$, so $P$ is $0$ or $\tfrac12$.

**What is automatic and what is not.** Moving the origin by one site multiplies $z$ by $e^{2\pi ix_0}=-1$, which shifts $P$ by $\tfrac12$. On the clean ring, $\delta\to-\delta$ is a translation by one site, so $\Delta P=\tfrac12$ there is exact and automatic (B4). The physics is in three other things:
- the values are pinned to $0$ and $\tfrac12$;
- $|z|$ is appreciable;
- both survive when translation is broken (item 4).

**Charge order.** $\langle M\rangle=0$ in any $C$-symmetric state, including an ordered finite ring. There the ground state is $(|A\rangle+|B\rangle)/\sqrt2$ of the two patterns, so the test is $\langle M^2\rangle$:
- **no order:** $\langle M^2\rangle\sim1/L$;
- **order:** $\langle M^2\rangle\to$ const, and the partner $(|A\rangle-|B\rangle)/\sqrt2$ comes down exponentially in $L$.

The same superposition has $z=0$, because the patterns contribute $z=\mp i$. So the CDW makes $P$ meaningless rather than changing its value.

## 3. How does the code compute it?

**Ring and polarisation** (notebook 6.1):
- `apply_T` is translation, with the sign $(-1)^{N-1}$ for the particle that wraps around.
- `ring_diagonals` gives $X$ and $M$ on every basis state.
- `ring_point` diagonalises $N$ and $N\pm1$, selects the ground state with Step 2's rule, labels levels with $C$ and $R$, and returns gaps, $z$, $P$, $\langle M\rangle$, $\langle M^2\rangle$ and the partner's data.

**Independent $V=0$ routes:**
- **`slater_z`** computes $\det(\Phi^\dagger D\Phi)\,e^{-i\pi(L-1)/2}$ from the filled orbitals of the $L\times L$ hopping matrix, with no Step 2 code.
- **`berry_product`** uses the fact that $\langle\phi_k|D|\phi_{k'}\rangle$ is nonzero only for $k'=k-2\pi/M$. The determinant is then $(-1)^{M-1}$ times $\prod_k\langle u_k|\mathrm{diag}(1,e^{i\pi/M})|u_{k-2\pi/M}\rangle$, the discrete Berry phase factor of the filled band with site $2m+1$ half a cell from $2m$. Its phase converges to the Zak phase $\gamma$, and $P=\gamma/2\pi-\tfrac14$ for odd $M$. The $-\tfrac14$ is the centre of the $-\tfrac12$-per-site background.

**Perturbations:**
- `hamiltonian_bonds` is Step 2's `hamiltonian` with explicit bonds and an on-site term.
- `paired_bonds` draws $\Delta J_j=\Delta J_{L-2-j}$; reflection maps bond $j$ to bond $L-2-j$, and the ring's closing bond to itself.
- `commutator_norms` measures which of $C$ and $R$ each perturbation keeps, instead of assuming it.

**Open chain:**
- `edge_lines` and `edge_weight` are Step 3's addition spectrum and window for explicit bonds and potentials.
- `doublet_profile` takes the two lowest half-filled states and the density of their $C$-breaking combination.

## 4. How do we know it is right?

| Id | What | Worst | Independent because |
|---|---|---|---|
| A1 | $TH(\delta)=H(-\delta)T$, every sector of $L=8,10$ | $3.6\times10^{-15}$ | even $N$ tests the sign of $T$ |
| A2 | $[C,H]=[R,H]=0$ on the ring | $1.3\times10^{-15}$ | operator maps of Step 2 |
| A3 | $\delta=V=0$: degenerate for $L=8,12,16$, unique for $10,14,18$ | $1.9\times10^{-14}$ | the reason for the sizes |
| A4 | every ring ground state unique, positive, $C=+1$ | min amplitude $-1.8\times10^{-16}$ | Perron-Frobenius |
| A5, A6 | residuals of every level used, ring and open | $3.4\times10^{-13}$ | |
| B1 | $V=0$: many-body $z$ = Slater determinant | $8.6\times10^{-15}$ | no Step 2 code |
| B2 | Slater determinant = Berry product, odd and even $M$ | $1.0\times10^{-14}$ | $k$-space route |
| B3, B3b | $z$ real, clean ring and under paired bonds | $1.3\times10^{-10}$ | symmetry; Lanczos vector accuracy |
| B4 | $z(-\delta)=-z(\delta)$ | $6.5\times10^{-11}$ | translation |
| B5 | $z=0$ at $\delta=0$ | $1.0\times10^{-10}$ | translation maps $H$ to itself |
| B6 | $\langle M\rangle=0$ | $6.4\times10^{-11}$ | $C$ |
| C1 | $\Delta_+=\Delta_-$ | $2.0\times10^{-13}$ | $C$ maps $N+1$ to $N-1$ |
| C2 | $L=14$: Lanczos = dense ($E_0$, $\Delta_0$, $z$, $\langle M^2\rangle$) | $1.1\times10^{-11}$ | solver independence |
| C3 | $L=18$: Lanczos seed independence | $1.5\times10^{-11}$ | |
| D1 | `hamiltonian_bonds` = `hamiltonian` | 0 | |
| D2, D4 | paired bonds keep $C$, $R$; $0.1n_0$ keeps neither (ring, open) | $1.6\times10^{-15}$ | computed commutators |
| D3 | $L=12$ edge weight and $\Delta_+$ = Step 3's saved values | $3.4\times10^{-12}$ | `data/step3_z0.csv` |

The symmetry checks B3-B6 are exact for exact eigenvectors. A Lanczos vector is off by about residual/gap, which reaches $10^{-10}$ where the partner gap is $2\times10^{-3}$, so their tolerance is $10^{-9}$.

## 5. What do the results show?

**Figures:**
- `step6_polarisation`: $P$ and $|z|$ at $V=0$, and $|z|\to1$;
- `step6_bulk_vs_V`: gaps, $L\Delta_0$, $|z|$, $P$ and the overlap against $V$;
- `step6_charge_order`: $L\langle M^2\rangle$, the partner splitting, ground state against partner;
- `step6_perturbations`: ring $P$ and the edge line under the perturbations;
- `step6_open_doublets`: open chains, edge doublet against CDW doublet.

**Ring, $\delta=-0.30$** (the spectra are identical for $+0.30$; $P=\tfrac12$ and $0$ at every point):

| $V$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| $\Delta_0$, $L=10/14/18$ | 1.68/1.47/1.37 | 1.61/1.50/1.46 | 1.32/1.25/1.23 | 0.96/0.89/0.87 | 0.64/0.53/0.49 | 0.40/0.27/0.21 | 0.25/0.13/0.076 | 0.10/0.032/0.011 | 0.046/0.0095/0.0021 |
| $\Delta_\pm$, $L=18$ | 0.69 | 1.00 | 1.29 | 1.56 | 1.87 | 2.35 | 3.04 | 4.76 | 6.63 |
| $\lvert z\rvert$, $L=10/14/18$ | 0.62/0.72/0.78 | 0.61/0.71/0.77 | 0.56/0.67/0.74 | 0.47/0.59/0.67 | 0.36/0.45/0.52 | 0.25/0.29/0.33 | 0.17/0.18/0.18 | 0.080/0.077/0.076 | 0.045/0.043/0.043 |
| $L\langle M^2\rangle$, $L=10/14/18$ | 0.50/0.50/0.50 | 0.63/0.63/0.63 | 0.80/0.81/0.81 | 1.05/1.11/1.14 | 1.36/1.59/1.75 | 1.66/2.16/2.62 | 1.91/2.61/3.31 | 2.18/3.05/3.91 | 2.31/3.23/4.14 |
| slope of $\ln\Delta_0$ against $L$ | -0.026 | -0.012 | -0.008 | -0.013 | -0.035 | -0.083 | -0.15 | -0.28 | -0.38 |

**Reading it:**
- **$V\le3$ is an insulator without charge order.** The gaps barely depend on $L$. $L\langle M^2\rangle$ saturates: the ratio $\langle M^2\rangle_{18}/\langle M^2\rangle_{10}$ is $0.556=10/18$ at $V\le1$ and $0.60$ at $V=3$. $|z|$ grows with $L$.
- **Crossover at $V\approx4.8$-$4.9$.** $L\Delta_0$ curves cross there, and the ground-state overlap between neighbouring $V$ is smallest at $V\approx4$-4.75, with a dip that deepens with $L$. All overlaps are positive, so this is not an eigenvector swap.
- **$V\gtrsim5$ is the CDW.** $\langle M^2\rangle$ is independent of $L$ (ratio 0.97-1.00 for $V\ge6$) and close to $\tfrac14$. The partner ($C=-1$, $R=-1$) has the same $\langle M^2\rangle$ and $4|\langle0|M|1\rangle|^2=0.85$-$0.91$, and its splitting closes exponentially. Its $C$-breaking combinations have $|z|=0.98$ and $P\approx\pm\tfrac14$: sharp but not quantised. The symmetric ground state averages them, $|z|\to0.04$-$0.08$, and $P$ means nothing.
- **Bond disorder:** $P$ stays exactly $\tfrac12$ and $0$ for all seeds, sizes and $V$, with $\Delta P=\tfrac12$ although translation is broken. $|z|$ and $\Delta_0$ change by a few percent (up to 9% and 14% at $V=8$).
- **End potential $0.1\,n_0$:** at $V=1$, $P$ shifts by $\pm0.0046$, $0.0031$, $0.0024$ ($L=10,14,18$); at $V=3$ by $\pm0.016$, $0.012$, $0.009$. Both go roughly as $1/L$, because a local potential moves a fixed charge. At $V=8$ it pins the pattern with site 0 empty for both dimerisations: $P=0.23$-$0.27$, $|z|=0.65$-$0.97$.
- **Edge line, open $L=12$, $\delta=-0.30$:**
  - clean: $\omega_{\rm edge}=0.029$, $z_0=0.67$ at $V=1$; $0.16$, $0.46$ at $V=3$;
  - paired bonds: $\omega_{\rm edge}=0.021$-$0.036$, $z_0=0.63$-$0.71$ at $V=1$; $0.15$-$0.18$, $0.43$-$0.49$ at $V=3$;
  - $0.1\,n_0$: $\omega_{\rm edge}=0.088$ and $0.23$, with $\langle n_0\rangle$ down to $0.16$ and $0.13$. The potential (0.1) exceeds the edge-pair splitting (0.045) and picks the end-empty state, while the bulk $P$ is unchanged.
- **Open chains:**
  - $\delta=-0.30$ has a near-degenerate pair at every $V$. Its $C$-breaking combination is an edge state at $V=0$ ($Q_4=0.46$, middle stagger 0.002). At $V=3$ it has $Q_4=0.27$ and an induced stagger of 0.23, 0.18, 0.13 at $L=12,16,20$. From $V\approx5$ it is a CDW (stagger 0.37-0.48, $Q_4\le0.09$).
  - $\delta=+0.30$ gets a CDW doublet at $V\gtrsim6$ and, at $V=8$, edge-window weight $z_0=0.24$ at $L=12$.

**A defensible statement:** *on rings $L\le18$ and open chains $L\le20$, for $\delta=\pm0.30$ and $V\le3$, the bulk is a gapped, charge-disordered insulator. Its polarisation is quantised to $\tfrac12$ for the weak-end bond pattern and to $0$ for the other, and the quantisation survives symmetry-preserving bond disorder. These bulk diagnostics support the edge interpretation of Steps 3-5. The CDW takes over between $V\approx4$ and 5. There $|z|$ collapses, the neutral gap closes exponentially, and the weak-end chain's edge pair turns into the CDW doublet. A boundary potential moves the edge line without touching the bulk $P$.*

## 6. Limits

- **Sizes.** The rings go only up to $L=18$ and the open chains up to 20. The crossover $V\approx4.8$-$4.9$ is a finite-size marker. The $L\Delta_0$ crossing assumes $\Delta_0\propto1/L$ at the transition. Neither the transition point nor its universality class is determined.
- **Only $\delta=\pm0.30$.** The crossover moves with $|\delta|$; at $\delta=0$ it is at $V=2$.
- **$P$ refers to a fixed origin and cell** (site 0, cell $(0,1)$). Only the difference between the two patterns, and its stability, is convention-free.
- **One disorder strength** ($0.05$), five seeds, two points (plus one in the CDW).
- **The edge line under perturbations is computed at $L=12$ only**, with Step 3's clean-ring window.
- **The $V=3$ staggering** around the edge charge decreases with $L$ on $L\le20$, but larger chains (Step 7) are needed to see how far it extends.

## 7. Files and how to reproduce

| File | Contents |
|---|---|
| [`project_A.ipynb`](../project_A.ipynb) | Step 6 sections 6.1-6.11 after Step 5 |
| [`report.tex`](../report.tex) | Step 6 section |
| `figures/step6_polarisation.pdf`, `step6_bulk_vs_V.pdf`, `step6_charge_order.pdf`, `step6_perturbations.pdf`, `step6_open_doublets.pdf` | figures |
| `data/step6_ring.csv` | every $(L,\delta,V)$ on the clean ring: levels, residuals, gaps, $z$, $P$, $\langle M\rangle$, $\langle M^2\rangle$, partner data, overlap |
| `data/step6_free.csv` | free-fermion $\lvert z\rvert$, $P$, $L=10$-202 |
| `data/step6_partner.csv`, `step6_cdw_combinations.csv` | partner table; $C$-breaking combinations |
| `data/step6_perturbations_ring.csv`, `step6_perturbations_edge.csv` | item 4, with commutator norms |
| `data/step6_open_doublets.csv` | open-chain doublets and densities |
| `data/step6_checks.csv`, `step6_run_info.json` | checks and run settings |

**To reproduce:** run the whole notebook. Step 6 uses only the Step 1-3 functions. It takes about 6 minutes, most of it for the ring scan (6.4), and prints `Step 6: 20/20 checks passed`.

## 8. Record of LLM assistance

Step 6 (code, data, figures, this summary and the report section) was produced by Claude (Anthropic) at Daniyaal's request, in one session on 25 November 2026, from the guide and the existing Steps 1-5. The Step 6 cells were first executed after the notebook's Step 0-3 definition cells only. On 27 November 2026 Daniyaal ran the whole notebook top to bottom (Python 3.12.16, numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, matplotlib 3.11.2): no errors, 198/198 checks (Steps 1-6: 46, 61, 14, 38, 19, 20). That run's outputs are the ones committed. Its Step 2-3 data are identical to the earlier commits; its Step 6 data agree with the first run to $10^{-10}$, except for quantities that depend on an arbitrary choice the solver makes (the partner vector inside the fourfold-degenerate level at $V=0$, and the overall sign of the partner, which only swaps the two labels of the $C$-breaking combinations), none of which changes a number quoted here.
