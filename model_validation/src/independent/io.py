"""The C++ file formats, and the linear-index encoding they depend on.

Everything the C++ stores sparsely is addressed by a single linear index into a
row-major multi-dimensional space (`getSubscriptsAsArray`,
`coretools/core/coretools/algorithms.h:693`). Row-major means the *last*
dimension varies fastest, which for this model is always the molecules tree:

    field (Y)             dims [n_species_leaves,    n_molecule_leaves]
                          linear = species_leaf    * n_molecule_leaves    + molecule_leaf

    species internal (Z)  dims [n_species_internals, n_molecule_leaves]
                          linear = species_internal * n_molecule_leaves    + molecule_leaf

    molecules internal (Z) dims [n_species_leaves,   n_molecule_internals]
                          linear = species_leaf    * n_molecule_internals + molecule_internal

A tree's own dimension carries *internal* nodes in its Z space and *leaf* nodes
in the field; the other dimension always carries leaves (`TTree::_initialize_Z`,
src/tree/TTree.cpp:123). The node orderings involved are defined in `indexing`.

File shapes, as the C++ readers expect them:

- field: 5 columns, `position Y_state species molecules fraction_of_one`. The
  reader (`TMarkovField::_read_Y_from_file`, src/TMarkovField.cpp:200) checks the
  column count but reads only `position` and `Y_state`. Every cell is written,
  not just the ones. So the field file is addressed purely positionally, and a
  change to leaf order changes what it means.
- internal states: 4 columns, `species molecules position Z_state`. The reader
  (`read_Z_from_file`, src/tree/io/read_Z.cpp:33) checks for 4 and resolves the
  cell from the two *name* columns, ignoring `position` -- which is what makes a
  file written before a node reordering still mean what it said (ADR-0004).
- observations: 2 columns naming one leaf per tree, one row per *positive* cell.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

from .indexing import TreeIndex


def write_tree(
    path: pathlib.Path, edges: list[tuple[str, str]], lengths: np.ndarray
) -> None:
    """Write a `child parent length` tree file."""
    pd.DataFrame(
        {
            "child": [c for c, _ in edges],
            "parent": [p for _, p in edges],
            "length": lengths,
        }
    ).to_csv(path, sep="\t", index=False)


def write_paper_counts(
    path: pathlib.Path, tree_name: str, leaf_names: list[str], counts: np.ndarray
) -> None:
    """Write the per-leaf paper counts that research effort is derived from."""
    pd.DataFrame({tree_name: leaf_names, "number_of_papers": counts}).to_csv(
        path, sep="\t", index=False
    )


def write_parameters(path: pathlib.Path, values: dict[str, object]) -> None:
    """Write a `name value` file, the format the C++ reads initial values from."""
    pd.DataFrame(
        {"name": list(values.keys()), "value": [str(v) for v in values.values()]}
    ).to_csv(path, sep="\t", index=False)


def write_field(
    path: pathlib.Path,
    field: np.ndarray,
    species: TreeIndex,
    molecules: TreeIndex,
) -> None:
    """Write the field as the C++ 5-column format, every cell included.

    `field` is indexed `[species_leaf, molecule_leaf]`.
    """
    n_molecule_leaves = molecules.n_leaves
    if field.shape != (species.n_leaves, n_molecule_leaves):
        raise ValueError(
            f"Field is {field.shape}, expected {(species.n_leaves, n_molecule_leaves)}."
        )

    flat = field.reshape(-1).astype(np.int8)
    species_ix, molecule_ix = np.divmod(np.arange(flat.size), n_molecule_leaves)

    pd.DataFrame(
        {
            "position": np.arange(flat.size),
            "Y_state": flat,
            "species": np.array(species.leaf_names())[species_ix],
            "molecules": np.array(molecules.leaf_names())[molecule_ix],
            # The C++ writes a posterior frequency here; a single draw has none,
            # so the state itself keeps the file self-consistent when eyeballed.
            "fraction_of_one": flat.astype(float),
        }
    ).to_csv(path, sep="\t", index=False)


def write_internal_states(
    path: pathlib.Path,
    states: np.ndarray,
    species: TreeIndex,
    molecules: TreeIndex,
    dimension: int,
) -> None:
    """Write one tree's internal states as the C++ 4-column format.

    `dimension` is 0 for the species tree and 1 for the molecules tree; that tree
    contributes internal nodes to `states`, the other contributes leaves.
    """
    if dimension == 0:
        row_names, col_names = species.internal_names(), molecules.leaf_names()
    else:
        row_names, col_names = species.leaf_names(), molecules.internal_names()

    if states.shape != (len(row_names), len(col_names)):
        raise ValueError(
            f"States are {states.shape}, expected {(len(row_names), len(col_names))}."
        )

    flat = states.reshape(-1).astype(np.int8)
    row_ix, col_ix = np.divmod(np.arange(flat.size), len(col_names))

    pd.DataFrame(
        {
            "species": np.array(row_names)[row_ix],
            "molecules": np.array(col_names)[col_ix],
            "position": np.arange(flat.size),
            "Z_state": flat,
        }
    ).to_csv(path, sep="\t", index=False)


def write_observations(
    path: pathlib.Path,
    observed: np.ndarray,
    species: TreeIndex,
    molecules: TreeIndex,
) -> None:
    """Write the positive cells of an observation as a 2-column sparse file."""
    species_ix, molecule_ix = np.nonzero(observed)
    pd.DataFrame(
        {
            "species": np.array(species.leaf_names())[species_ix],
            "molecules": np.array(molecules.leaf_names())[molecule_ix],
        }
    ).to_csv(path, sep="\t", index=False)


def read_field(path: pathlib.Path, n_molecule_leaves: int) -> np.ndarray:
    """Read a 5-column field file back into a `[species_leaf, molecule_leaf]` array."""
    frame = pd.read_csv(path, sep="\t", usecols=["position", "Y_state"])
    frame = frame.sort_values("position")
    return frame["Y_state"].to_numpy(dtype=bool).reshape(-1, n_molecule_leaves)
