"""Compute the per-branch state-transition density of the simulated latent field.

For each tree, a "branch" is a (parent, child) edge, and the field is replicated
across cliques (= leaves of the *other* tree). A branch carries information about
its length only when parent and child differ in state. This script reports the
fraction of (branch, clique) pairs where parent state != child state:

    - internal node states come from acol_simulated_Z_<tree>.txt
    - leaf node states come from the Y field (acol_simulated_Y.txt)

A density near 0 means the latent process is almost frozen, so individual branch
lengths have essentially no likelihood signal and cannot be inferred.
"""

import pathlib

import click
import pandas as pd


def _load_edges(path: pathlib.Path) -> tuple[pd.DataFrame, set[str]]:
    """Return (edges[child, parent], set_of_leaf_ids)."""
    edges = pd.read_csv(path, sep="\t")[["child", "parent"]].astype(str)
    leaves = set(edges["child"]) - set(edges["parent"])
    return edges, leaves


def _build_states(
    z_path: pathlib.Path,
    y_path: pathlib.Path,
    node_col: str,
    clique_col: str,
) -> pd.DataFrame:
    """Long table [clique, node, state] for one tree.

    Internal-node states come from the Z file; leaf states from the Y field.
    `node_col` is this tree's dimension, `clique_col` is the other tree's leaf.
    """
    z = pd.read_csv(z_path, sep="\t", usecols=[node_col, clique_col, "Z_state"])
    z = z.rename(columns={node_col: "node", clique_col: "clique", "Z_state": "state"})

    y = pd.read_csv(y_path, sep="\t", usecols=[node_col, clique_col, "Y_state"])
    y = y.rename(columns={node_col: "node", clique_col: "clique", "Y_state": "state"})

    states = pd.concat([z, y], ignore_index=True)
    states["node"] = states["node"].astype(str)
    states["clique"] = states["clique"].astype(str)
    states["state"] = states["state"].astype(int)
    return states


def _transition_stats(
    edges: pd.DataFrame,
    leaves: set[str],
    states: pd.DataFrame,
) -> dict[str, float]:
    """Join parent/child states across all cliques and summarise transitions."""
    # attach child state (one row per branch x clique)
    child = edges.merge(
        states.rename(columns={"node": "child", "state": "state_child"}),
        on="child",
    )
    # attach parent state for the same clique
    full = child.merge(
        states.rename(columns={"node": "parent", "state": "state_parent"}),
        on=["parent", "clique"],
    )

    full["differs"] = full["state_child"] != full["state_parent"]
    full["is_leaf_branch"] = full["child"].isin(leaves)

    n_cliques = states["clique"].nunique()
    n_branches = len(edges)
    expected_pairs = n_branches * n_cliques

    leaf = full[full["is_leaf_branch"]]
    internal = full[~full["is_leaf_branch"]]

    return {
        "n_cliques": n_cliques,
        "n_branches": n_branches,
        "n_leaf_branches": int(edges["child"].isin(leaves).sum()),
        "n_internal_branches": int((~edges["child"].isin(leaves)).sum()),
        "matched_pairs": len(full),
        "coverage": len(full) / expected_pairs if expected_pairs else float("nan"),
        "density_overall": full["differs"].mean() if len(full) else float("nan"),
        "density_leaf": leaf["differs"].mean() if len(leaf) else float("nan"),
        "density_internal": internal["differs"].mean()
        if len(internal)
        else float("nan"),
    }


def _report(tree_name: str, stats: dict[str, float]) -> None:
    click.echo(f"\n=== {tree_name} tree ===")
    click.echo(
        f"  cliques (leaves of other tree): {stats['n_cliques']}\n"
        f"  branches: {stats['n_branches']} "
        f"({stats['n_internal_branches']} internal, {stats['n_leaf_branches']} leaf)\n"
        f"  coverage (matched / expected pairs): {stats['coverage']:.3f} "
        f"({stats['matched_pairs']} pairs)\n"
        f"  transition density  overall : {stats['density_overall']:.4f}\n"
        f"                      internal: {stats['density_internal']:.4f}\n"
        f"                      leaf    : {stats['density_leaf']:.4f}"
    )


@click.command()
@click.argument("scenario_dir", type=click.Path(exists=True, file_okay=False))
def main(scenario_dir: str) -> None:
    """Report per-branch transition density for both trees of a scenario."""
    base = pathlib.Path(scenario_dir)
    y_path = base / "acol_simulated_Y.txt"
    if not y_path.exists():
        raise click.ClickException(f"Missing Y field: {y_path}")

    configs = [
        ("species", "species.txt", "acol_simulated_Z_species.txt", "species", "molecules"),
        ("molecules", "molecules.txt", "acol_simulated_Z_molecules.txt", "molecules", "species"),
    ]

    for tree_name, edge_file, z_file, node_col, clique_col in configs:
        edge_path = base / edge_file
        z_path = base / z_file
        if not edge_path.exists() or not z_path.exists():
            click.echo(f"\nSkipping {tree_name}: missing {edge_file} or {z_file}.")
            continue

        edges, leaves = _load_edges(edge_path)
        states = _build_states(z_path, y_path, node_col, clique_col)
        stats = _transition_stats(edges, leaves, states)
        _report(tree_name, stats)


if __name__ == "__main__":
    main()
