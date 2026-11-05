"""A tiny LOTUS+MS scenario, and the output of one inference rung on it.

The files have the names and columns the simulator and the C++ binary write, so
`validate_independent.score` reads them as it reads a real run. The scenario is
small enough to check by hand:

            m1  m2  m3  m4        species s1 and s2 have a run, s3 has none.
    s1       1   0   0   0        LOTUS records: (s1, m1), (s1, m3), (s3, m4).
    s2       0   1   0   0        (s1, m3) and (s3, m4) are false records.
    s3       1   0   0   0

The field false-presence score counts the cells that are absent, have no LOTUS record
and lie in a species with a run. The rung's field posterior gives four of them the
values in `COUNTED` and omits the fifth, (s2, m3), which is a posterior of 0. The five
mean 0.2. Every other listed cell has a high posterior, so a score that counted one of
them would move.

The rung's trace draws every parameter at its truth except those in `MISSED`, so the
coverage scores differ from row to row. The assignment posterior puts `molecule_mass`
on the true assignment of the one molecule feature and `unknown_mass` on that of each
of the two unknown features.
"""

from __future__ import annotations

import json
import pathlib
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

SPECIES = ("s1", "s2", "s3")
MOLECULES = ("m1", "m2", "m3", "m4")
PRESENT = {("s1", "m1"), ("s2", "m2"), ("s3", "m1")}
LOTUS = (("s1", "m1"), ("s1", "m3"), ("s3", "m4"))
RUN_SPECIES = ("s1", "s2")

COUNTED = {
    ("s1", "m2"): 0.10,
    ("s1", "m4"): 0.20,
    ("s2", "m1"): 0.30,
    ("s2", "m4"): 0.40,
}
"""The listed cells that the false-presence score counts. (s2, m3) is the fifth, unlisted."""
UNCOUNTED = {
    ("s1", "m1"): 1.0,  # present
    ("s2", "m2"): 1.0,  # present
    # Present, no run. Under (s3, m2), so the AUC for species without runs is 2/3.
    ("s3", "m1"): 0.7,
    ("s1", "m3"): 0.9,  # absent, with a false LOTUS record
    ("s3", "m2"): 0.8,  # absent, no LOTUS record, species without a run
}

# (species, feature, candidates, true molecule or None for an unknown feature)
FEATURES = (
    ("s1", "feature_0", ("m1", "m2"), "m1"),
    ("s1", "feature_1", ("m3", "m4"), None),
    ("s2", "feature_0", ("m2", "m4"), None),
)
EXTRACT = "extract_1"
UNLISTED = "<unlisted>"

FILTER_INDICES = {"m1": 100, "m2": 200, "m3": 50, "m4": 150}
CONTAMINATION = 0.01
# One alpha and one log nu per clique. A species-tree clique is named by a molecule
# leaf and a molecule-tree clique by a species leaf. The values differ, so that the
# validator can correlate the posterior means with them.
SPECIES_TREE_TRUTH = {
    **{f"species_alpha_{m}": 0.2 * (i + 1) for i, m in enumerate(MOLECULES)},
    **{f"species_log_nu_{m}": -0.2 * (i + 1) for i, m in enumerate(MOLECULES)},
}
MOLECULE_TREE_TRUTH = {
    **{f"molecules_alpha_{s}": 0.25 * (i + 1) for i, s in enumerate(SPECIES)},
    **{f"molecules_log_nu_{s}": -0.25 * (i + 1) for i, s in enumerate(SPECIES)},
}
MISSED = frozenset(
    {
        "species_alpha_m4",
        "species_log_nu_m2",
        "species_log_nu_m3",
        "species_log_nu_m4",
        "molecules_alpha_s2",
        "molecules_alpha_s3",
        "filter_proba_m4",
        "contamination_proba",
    }
)
"""The parameters whose draws miss their truth, by `MISS`. Every other interval holds it.

The coverage is then 3/4 for species alpha, 1/4 for species log nu, 1/3 for molecules
alpha and 1 for molecules log nu. The filter probabilities of the molecules present in
no run (m3, m4) have 1/2, and those of molecules present in one run (m1, m2) have 1.
Kappa misses.
"""
MISS = 5.0
N_ITERATIONS = 60
FEATURE_KEY = ["species", "extract", "feature"]


def _write_tsv(path: pathlib.Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def _name_value(path: pathlib.Path, values: Mapping[str, object]) -> None:
    _write_tsv(
        path, pd.DataFrame({"name": list(values), "value": list(values.values())})
    )


def write_scenario(
    base: pathlib.Path,
    *,
    seed: int = 1,
    layout: str = "dense",
    lotus: Sequence[tuple[str, str]] = LOTUS,
) -> pathlib.Path:
    """The scenario's truth and data, as `simulate_independent.py` writes them.

    `lotus` lists the (species, molecule) cells that have a LOTUS record.
    """
    base.mkdir(parents=True, exist_ok=True)
    meta = {
        "seed": seed,
        "n_burnin_rounds": 0,
        "burnin": 0,
        "iterations": N_ITERATIONS,
        "mass_spec": {"layout": layout, "n_runs": len(RUN_SPECIES)},
    }
    (base / "meta.json").write_text(json.dumps(meta) + "\n")

    _name_value(base / "truth_species.txt", SPECIES_TREE_TRUTH)
    _name_value(base / "truth_molecules.txt", MOLECULE_TREE_TRUTH)
    _name_value(
        base / "simulated_parameters.txt",
        {
            **{f"filter_proba_{m}": index for m, index in FILTER_INDICES.items()},
            "contamination_proba": CONTAMINATION,
        },
    )

    cells = [(s, m) for s in SPECIES for m in MOLECULES]
    _write_tsv(
        base / "simulated_Y.txt",
        pd.DataFrame(
            {
                "position": range(len(cells)),
                "Y_state": [int(cell in PRESENT) for cell in cells],
                "species": [s for s, _ in cells],
                "molecules": [m for _, m in cells],
                "fraction_of_one": [float(cell in PRESENT) for cell in cells],
            }
        ),
    )
    _write_tsv(
        base / "simulated_lotus.tsv",
        pd.DataFrame(list(lotus), columns=["species", "molecules"]),
    )

    data_rows = [
        (species, EXTRACT, feature, molecule, 0.0, 0.0)
        for species, feature, candidates, _ in FEATURES
        for molecule in candidates
    ]
    _write_tsv(
        base / "simulated_ms_data.tsv",
        pd.DataFrame(
            data_rows,
            columns=[
                *FEATURE_KEY,
                "molecule",
                "log_correct_density",
                "log_incorrect_density",
            ],
        ),
    )
    _write_tsv(
        base / "simulated_ms_assignment.tsv",
        _assignment_frame(molecule_mass=1.0, unknown_mass=1.0),
    )
    return base


def _assignment_frame(*, molecule_mass: float, unknown_mass: float) -> pd.DataFrame:
    """One row per candidate and one for `<unlisted>`, for every feature.

    A feature of a molecule puts `molecule_mass` on that molecule and the rest on the
    unknown molecule. An unknown feature puts `unknown_mass` on the unknown molecule
    and the rest on its first candidate. Both at 1 are the true assignment. No feature
    lists an outside structure, so `<unlisted>` is the whole unknown molecule.
    """
    rows = []
    for species, feature, candidates, truth in FEATURES:
        explanations = [*candidates, UNLISTED]
        true = UNLISTED if truth is None else truth
        other = UNLISTED if truth is not None else candidates[0]
        mass = unknown_mass if truth is None else molecule_mass
        for molecule in explanations:
            fraction = (
                mass if molecule == true else (1.0 - mass if molecule == other else 0.0)
            )
            rows.append((species, EXTRACT, feature, molecule, fraction))
    return pd.DataFrame(rows, columns=[*FEATURE_KEY, "molecule", "fraction"])


def write_rung(
    base: pathlib.Path,
    rung: str,
    *,
    field_scale: float = 1.0,
    molecule_mass: float = 1.0,
    unknown_mass: float = 1.0,
) -> pathlib.Path:
    """The output of one inference rung on the scenario in `base`.

    `field_scale` multiplies the posterior on the cells that the false-presence score
    counts. `molecule_mass` and `unknown_mass` are the posterior mass on the true
    assignment of a molecule feature and of an unknown feature. The trace holds every
    parameter at its truth, except those in `MISSED`.
    """
    rung_dir = base / rung

    posterior = {
        **{cell: value * field_scale for cell, value in COUNTED.items()},
        **UNCOUNTED,
    }
    _write_tsv(
        rung_dir / "acol_Y_posterior.tsv",
        pd.DataFrame(
            {
                "position": range(len(posterior)),
                "Y_state": [int(cell in PRESENT) for cell in posterior],
                "species": [s for s, _ in posterior],
                "molecules": [m for _, m in posterior],
                "fraction_of_one": list(posterior.values()),
            }
        ),
    )
    _write_tsv(
        rung_dir / "acol_assignment_posterior.txt",
        _assignment_frame(molecule_mass=molecule_mass, unknown_mass=unknown_mass),
    )

    truth = {
        **SPECIES_TREE_TRUTH,
        **MOLECULE_TREE_TRUTH,
        **{f"filter_proba_{m}": float(index) for m, index in FILTER_INDICES.items()},
        "contamination_proba": CONTAMINATION,
    }
    drawn = {
        name: value + (MISS if name in MISSED else 0.0) for name, value in truth.items()
    }
    trace = pd.DataFrame(
        {name: np.full(N_ITERATIONS, value) for name, value in drawn.items()}
    )
    _write_tsv(rung_dir / "acol_trace.txt", trace)
    # The validator cross-checks its burn-in trim against the binary's own means.
    _write_tsv(
        rung_dir / "acol_meanVar.txt",
        pd.DataFrame(
            {
                "name": list(drawn),
                "posterior_mean": list(drawn.values()),
                "posterior_variance": 0.0,
            }
        ),
    )
    return rung_dir
