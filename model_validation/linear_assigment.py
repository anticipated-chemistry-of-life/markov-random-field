"""
Mass Spectrometry Feature-to-Molecule Assignment — True Sparse
==============================================================
Solves a rectangular assignment problem entirely in sparse format.
No dense matrix is ever constructed — suitable for 1000s × 1000s problems.

Solver
------
scipy.sparse.csgraph.min_weight_full_bipartite_matching
  Implements LAPJVsp (Jonker-Volgenant sparse), operates directly on
  scipy CSR sparse arrays. Rectangular matrices are handled natively.

Input CSV  (COO edge list)
--------------------------
    feature,molecule,score
    Feature_0001,Mol_A,0.823
    Feature_0001,Mol_C,0.341
    Feature_0002,Mol_B,0.910

  - One row per existing candidate pair.
  - score : match score / probability  (higher = better, must be > 0).
  - Absent pairs are simply absent — no placeholder needed.

Usage
-----
    python ms_assignment_sparse.py                    # synthetic demo
    python ms_assignment_sparse.py edges.csv          # your data
    python ms_assignment_sparse.py edges.csv 0.15     # custom min-score threshold

Output
------
    assignment_results.csv  with columns:
        Feature, Assigned_Molecule, Score, Low_Confidence, Unassigned
"""

import sys
import time

import numpy as np
import pandas as pd
from scipy.sparse import coo_array
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

# ---------------------------------------------------------------------------
# 1.  Synthetic COO data generator
# ---------------------------------------------------------------------------


def make_synthetic_edges(
    n_features: int = 2000,
    n_molecules: int = 150000,
    avg_candidates: int = 8,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Generate a realistic sparse MS edge list.
    Each feature gets ~avg_candidates randomly sampled molecule candidates.
    Density ≈ avg_candidates / n_molecules  (e.g. 8/1500 ≈ 0.53 %).
    """
    rng = np.random.default_rng(seed)

    rows, cols, scores = [], [], []
    for f in range(n_features):
        k = min(int(rng.poisson(avg_candidates)), n_molecules)
        if k == 0:
            continue
        chosen = rng.choice(n_molecules, size=k, replace=False)
        rows.extend([f] * k)
        cols.extend(chosen.tolist())
        scores.extend(rng.beta(2, 3, size=k).tolist())

    return pd.DataFrame(
        {
            "feature": [f"Feature_{i + 1:05d}" for i in rows],
            "molecule": [f"Mol_{j + 1:05d}" for j in cols],
            "score": np.round(scores, 6),
        }
    )


# ---------------------------------------------------------------------------
# 2.  Build sparse biadjacency matrix (no dense array, ever)
# ---------------------------------------------------------------------------


def edges_to_sparse(edges: pd.DataFrame):
    """
    Convert COO edge list to a scipy sparse biadjacency array.

    Returns
    -------
    biadj      : coo_array  (features × molecules), scores as weights
    feat_index : list[str]  row index → feature name
    mol_index  : list[str]  col index → molecule name
    """
    feat_index = sorted(edges["feature"].unique())
    mol_index = sorted(edges["molecule"].unique())

    f_map = {name: i for i, name in enumerate(feat_index)}
    m_map = {name: i for i, name in enumerate(mol_index)}

    row = edges["feature"].map(f_map).to_numpy(dtype=np.int32)
    col = edges["molecule"].map(m_map).to_numpy(dtype=np.int32)
    data = edges["score"].to_numpy(dtype=np.float64)

    shape = (len(feat_index), len(mol_index))
    biadj = coo_array((data, (row, col)), shape=shape)
    return biadj, feat_index, mol_index


# ---------------------------------------------------------------------------
# 3.  Solver
# ---------------------------------------------------------------------------


def solve_assignment(edges: pd.DataFrame, min_score: float = 0.0) -> pd.DataFrame:
    n_f = edges["feature"].nunique()
    n_m = edges["molecule"].nunique()

    print(f"  Features  : {n_f:,}")
    print(f"  Molecules : {n_m:,}")
    print(
        f"  Edges     : {len(edges):,}  "
        f"(density {100 * len(edges) / (n_f * n_m):.3f} %)"
    )

    biadj, feat_index, mol_index = edges_to_sparse(edges)

    # min_weight_full_bipartite_matching minimises — pass maximize=True
    # to maximise scores directly. No negation, no dense padding needed.
    print("  Building CSR …", flush=True)
    biadj_csr = biadj.tocsr()

    print("  Solving (LAPJVsp, fully sparse) …", flush=True)
    t0 = time.perf_counter()
    row_ind, col_ind = min_weight_full_bipartite_matching(biadj_csr, maximize=True)
    elapsed = time.perf_counter() - t0
    print(f"  Done in {elapsed:.2f} s")

    # Build a lookup from assigned (row, col) pairs → score
    assigned = dict(zip(row_ind, col_ind))

    records = []
    for r, feat in enumerate(feat_index):
        if r in assigned:
            c = assigned[r]
            score = float(biadj_csr[r, c])
            records.append(
                {
                    "Feature": feat,
                    "Assigned_Molecule": mol_index[c],
                    "Score": round(score, 6),
                    "Low_Confidence": score < min_score,
                    "Unassigned": False,
                }
            )
        else:
            records.append(
                {
                    "Feature": feat,
                    "Assigned_Molecule": None,
                    "Score": 0.0,
                    "Low_Confidence": False,
                    "Unassigned": True,
                }
            )

    return pd.DataFrame(records).sort_values("Feature").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 4.  Reporting
# ---------------------------------------------------------------------------


def print_results(result: pd.DataFrame) -> None:
    n_total = len(result)
    n_assigned = (~result["Unassigned"]).sum()
    n_low = result["Low_Confidence"].sum()
    total_sc = result["Score"].sum()
    mean_sc = result.loc[~result["Unassigned"], "Score"].mean()

    print(f"\n── Assignment summary ──")
    print(f"  Features total               : {n_total:,}")
    print(
        f"  Assigned                     : {n_assigned:,}  ({100 * n_assigned / n_total:.1f} %)"
    )
    print(f"  Unassigned (no candidate)    : {n_total - n_assigned:,}")
    print(f"  Low-confidence (< threshold) : {n_low:,}")
    print(f"  Total assignment score       : {total_sc:.4f}")
    print(f"  Mean score (assigned only)   : {mean_sc:.4f}")
    print(f"\nFirst 20 rows:")
    print(result.head(20).to_string(index=False))


# ---------------------------------------------------------------------------
# 5.  Entry point
# ---------------------------------------------------------------------------


def load_edges(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.lower().str.strip()
    rename = {}
    for col in df.columns:
        if col in ("feat", "feature_id", "feature"):
            rename[col] = "feature"
        elif col in ("mol", "molecule_id", "molecule", "compound"):
            rename[col] = "molecule"
        elif col in ("score", "probability", "prob", "p", "value"):
            rename[col] = "score"
    df = df.rename(columns=rename)[["feature", "molecule", "score"]]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    if (df["score"] <= 0).any():
        print(
            f"  Warning: {(df['score'] <= 0).sum():,} edges with score ≤ 0 removed "
            "(zero weights are not supported by the sparse solver)."
        )
        df = df[df["score"] > 0]
    return df


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else None
    min_score = float(sys.argv[2]) if len(sys.argv) > 2 else 0.10

    if csv_path:
        print(f"Loading edges from: {csv_path}")
        edges = load_edges(csv_path)
    else:
        print("No CSV provided — generating synthetic sparse MS data …")
        edges = make_synthetic_edges(
            n_features=2000, n_molecules=150000, avg_candidates=10
        )
        edges.to_csv("synthetic_ms_edges.csv", index=False)
        print("Edge list saved → synthetic_ms_edges.csv")

    print(f"Min-score threshold : {min_score}\n")

    result = solve_assignment(edges, min_score=min_score)
    print_results(result)

    out = "assignment_results.csv"
    result.to_csv(out, index=False)
    print(f"\nFull results saved → {out}")


if __name__ == "__main__":
    main()
