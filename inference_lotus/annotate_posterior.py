"""Annotate the Y posterior with LOTUS metadata.

Reads the Y posterior (``test_out/acol_Y_posterior.tsv``), keeps the rows where
``fraction_of_one`` exceeds ``FRACTION_THRESHOLD``, and enriches each surviving
organism/molecule pair with:

- the organism lineage: one column per rank above the modelled organism rank
- the molecule lineage: one column per NPClassifier rank above the modelled
  molecule rank
- the molecule's 2D SMILES, only when the model ran on individual molecules
  (inchikey_2D); a coarser rank has no single structure
- a flag telling whether the pair is already known, i.e. present in
  ``lotus.tsv`` (the pairs the tree/model were built from)

The posterior's ``species``/``molecules`` columns hold nodes of whichever ranks
``data_prep.py`` built the trees at, so ``--taxon-level`` and ``--molecule-level``
take the same values and MUST match the ``data_prep.py`` run that produced the
inputs. Metadata is pulled from ``data/lotus.csv.gz``, deduplicated per node of
the chosen ranks, and joined onto the filtered posterior. The result is written
to ``test_out/acol_Y_posterior_annotated.tsv``.
"""

import click
import polars as pl
from data_prep import (
    MOLECULE_COLUMNS,
    NA_VALUES,
    TAXON_COLUMNS,
    MoleculeLevel,
    TaxonLevel,
    columns_up_to,
)

Y_POSTERIOR_PATH = "test_out/acol_Y_posterior.tsv"
LOTUS_RAW_PATH = "data/lotus.csv.gz"
LOTUS_PAIRS_PATH = "lotus.tsv"
OUT_PATH = "test_out/acol_Y_posterior_annotated.tsv"

FRACTION_THRESHOLD = 0.8


def _underscored(col: str) -> pl.Expr:
    # Node names throughout the tree/tsv files have spaces replaced by underscores.
    return pl.col(col).str.replace_all(" ", "_").alias(col)


def lineage(
    lotus_raw: pl.DataFrame,
    columns: list[str],
    leaf_name: str,
    extra: list[str],
) -> pl.DataFrame:
    """One row per node of the last rank in `columns`, carrying its ancestors and `extra`.

    The leaf column is renamed to `leaf_name`, the posterior's column for that tree.
    """
    leaf = columns[-1]
    return (
        lotus_raw.select([*columns, *extra])
        .with_columns([_underscored(c) for c in columns])
        .drop_nulls(subset=[leaf])
        .unique(subset=[leaf], keep="first")
        .rename({leaf: leaf_name})
    )


@click.command()
@click.option(
    "--taxon-level",
    type=click.Choice([level.name.lower() for level in TaxonLevel], case_sensitive=False),
    default=TaxonLevel.SPECIES.name.lower(),
    show_default=True,
    help="Organism taxonomy rank data_prep.py built the inputs at.",
)
@click.option(
    "--molecule-level",
    type=click.Choice(
        [level.name.lower() for level in MoleculeLevel], case_sensitive=False
    ),
    default=MoleculeLevel.MOLECULE.name.lower(),
    show_default=True,
    help="Molecule taxonomy rank data_prep.py built the inputs at.",
)
def main(taxon_level: str, molecule_level: str) -> None:
    taxon = TaxonLevel[taxon_level.upper()]
    molecule = MoleculeLevel[molecule_level.upper()]

    taxon_columns = columns_up_to(TAXON_COLUMNS, taxon.value)
    molecule_columns = columns_up_to(MOLECULE_COLUMNS, molecule.value)
    # A structure only exists for an individual molecule, not for a class of them.
    smiles_columns = ["smiles"] if molecule is MoleculeLevel.MOLECULE else []

    posterior = pl.read_csv(Y_POSTERIOR_PATH, separator="\t").filter(
        pl.col("fraction_of_one") > FRACTION_THRESHOLD
    )

    lotus_raw = (
        pl.read_csv(
            LOTUS_RAW_PATH,
            columns=[
                "structure_inchikey",
                "structure_smiles_2D",
                *[c for c in MOLECULE_COLUMNS if c != MoleculeLevel.MOLECULE.value],
                *TAXON_COLUMNS,
            ],
            null_values=NA_VALUES,
            infer_schema_length=10000,
        )
        .with_columns(
            pl.col("structure_inchikey").str.slice(0, 14).alias("inchikey_2D")
        )
        .rename({"structure_smiles_2D": "smiles"})
    )

    species_meta = lineage(lotus_raw, taxon_columns, "species", extra=[])
    molecule_meta = lineage(
        lotus_raw, molecule_columns, "molecules", extra=smiles_columns
    )

    lotus_pairs = (
        pl.read_csv(LOTUS_PAIRS_PATH, separator="\t")
        .unique()
        .with_columns(pl.lit(True).alias("already_in_lotus"))
    )

    annotated = (
        posterior.join(species_meta, on="species", how="left")
        .join(molecule_meta, on="molecules", how="left")
        .join(lotus_pairs, on=["species", "molecules"], how="left")
        .with_columns(pl.col("already_in_lotus").fill_null(False))
        .select(
            [
                "position",
                "Y_state",
                "fraction_of_one",
                *taxon_columns[:-1],
                "species",
                "molecules",
                *smiles_columns,
                *molecule_columns[:-1],
                "already_in_lotus",
            ]
        )
        .sort("fraction_of_one", descending=True)
    )

    annotated.write_csv(OUT_PATH, separator="\t")
    print(f"wrote {annotated.height} rows to {OUT_PATH}")


if __name__ == "__main__":
    main()
