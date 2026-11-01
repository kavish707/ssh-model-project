"""Can the Step 2 checks fail? Break the code one way at a time and list which checks notice.

Runs the notebook's Step 2 cells up to the L = 8 gate on a modified copy of the source, in a temporary directory,
so nothing in the repository changes. Run from the repository root: python checks/step2_mutation_test.py
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

NOTEBOOK = Path(__file__).resolve().parent.parent / "project_A.ipynb"
TIMEOUT_S = 900

MUTATIONS = {
    "interaction_x2": ("(H + sp.diags(V * interaction_diagonal", "(H + sp.diags(2 * V * interaction_diagonal",
                       "interaction twice too strong"),
    "interaction_no_half": ("x = occupations(basis, L) - 0.5", "x = occupations(basis, L) - 0.0",
                            "n_j n_(j+1) instead of (n_j - 1/2)(n_(j+1) - 1/2)"),
    "ring_no_int_bond": ("return total + x[:, -1] * x[:, 0] if periodic else total", "return total",
                         "ring: interaction misses the closing bond"),
    "ring_spin_bond": ("def hamiltonian(L, N, delta, V, periodic=False, closing_sign=True):",
                       "def hamiltonian(L, N, delta, V, periodic=False, closing_sign=False):",
                       "ring: closing bond without its fermion sign"),
    "bond_flip": ("for j, Jj in enumerate(bonds(L, delta, periodic)):\n        k = (j + 1) % L\n        for a, b in",
                  "for j, Jj in enumerate(bonds(L, -delta, periodic)):\n        k = (j + 1) % L\n        for a, b in",
                  "delta -> -delta in hamiltonian()"),
    "sign_on_original": ("jw_sign(states, b) * jw_sign(mid, a)", "jw_sign(states, b) * jw_sign(states, a)",
                         "second sign evaluated on the original state, not the intermediate one"),
    "sign_off_by_one": ("states & ((1 << j) - 1)", "states & ((1 << (j + 1)) - 1)", "sign string includes site j itself"),
    "c_no_string": ("sign = jw_sign(src[col], j) if string else np.ones(len(col), dtype=np.int64)",
                    "sign = np.ones(len(col), dtype=np.int64)", "c_j without the sign string"),
    "R_no_sign": ("= (-1) ** (N * (N - 1) // 2) * psi", "= psi", "R without (-1)^(N(N-1)/2)"),
    "C_wrong_sign": ("sector_basis(L, N) ^ ((1 << L) - 1))] = psi", "sector_basis(L, N) ^ ((1 << L) - 1))] = -psi",
                     "C with overall sign -1"),
    "lanczos_SM": ('which="SA"', 'which="SM", maxiter=50', "Lanczos targets smallest |E| instead of lowest E"),
    "lanczos_tol": ("tol=0, v0=v0", "tol=1e-2, v0=v0", "Lanczos stopped far too early"),
    "jw_no_parity": ('string = {} if k == j + 1 else {l: "Z" for l in range(1, L - 1)}', "string = {}",
                     "Jordan-Wigner ring bond without the parity string"),
    "select_wrong_C": ("plus = np.flatnonzero(vals > 0.5)", "plus = np.flatnonzero(vals < -0.5)",
                       "state selection keeps the C = -1 vector"),
    "select_no_rule": ("size = int(np.sum(E - E[0] < ETA_DEGENERATE))", "size = 1",
                       "state selection silently accepts the first vector"),
    "jw_bit_order": ("for j in reversed(range(L)):\n        out = sp.kron", "for j in range(L):\n        out = sp.kron",
                     "Jordan-Wigner bit order reversed (expected harmless: H is reflection symmetric, so same matrix)"),
}


def run_one(mutation):
    """Execute the Step 0 setup, Step 1 functions and Step 2 cells through the L = 8 gate; print failed check ids."""
    os.environ["MPLBACKEND"] = "Agg"
    os.chdir(tempfile.mkdtemp())
    cells = [c for c in json.loads(NOTEBOOK.read_text())["cells"] if c["cell_type"] == "code"]
    first = next(i for i, c in enumerate(cells) if "Step 2 parameters" in "".join(c["source"]))
    last = next(i for i, c in enumerate(cells) if "def refinement_triggers" in "".join(c["source"]))
    namespace, errors, applied = {}, [], mutation == "baseline"
    for cell in cells[:3] + cells[first:last]:
        source = "".join(cell["source"]).replace("from IPython.display import display", "display = lambda *a, **k: None")
        if mutation != "baseline" and MUTATIONS[mutation][0] in source:
            source, applied = source.replace(*MUTATIONS[mutation][:2]), True
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(source, "cell", "exec"), namespace)
        except Exception as error:
            errors.append(type(error).__name__)
    if not applied:
        raise SystemExit(f"mutation {mutation} matches nothing in the notebook source")
    failed = sorted({name.split()[1].rstrip(":") for name, c in namespace["CHECKS"].items()
                     if name.startswith("S2 ") and not c["pass"]})
    print(json.dumps({"failed": failed, "errors": sorted(set(errors))}))


def main():
    """Run the unmodified code and every mutation, each in its own process, and print a table."""
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    missed = []
    for mutation in ["baseline", *MUTATIONS]:
        try:
            out = subprocess.run([sys.executable, __file__, mutation], capture_output=True, text=True, env=env,
                                 timeout=TIMEOUT_S, check=True).stdout
            result = json.loads(out.splitlines()[-1])
            caught = result["failed"] + [f"raises {e}" for e in result["errors"] if e != "RuntimeError"]
        except subprocess.TimeoutExpired:
            caught = [f"no result in {TIMEOUT_S} s"]
        if mutation != "baseline" and not caught:
            missed.append(mutation)
        what = "unmodified code" if mutation == "baseline" else MUTATIONS[mutation][2]
        print(f"{mutation:20s} {len(caught):3d}  {' '.join(caught):60s}  {what}", flush=True)
    print(f"\nnot caught: {missed or 'none'}")


if __name__ == "__main__":
    main() if len(sys.argv) == 1 else run_one(sys.argv[1])
