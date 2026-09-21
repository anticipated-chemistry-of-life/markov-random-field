"""Annotate the Y posterior with LOTUS metadata.

Reads the Y posterior (``test_out/acol_Y_posterior.txt``), keeps the rows where
``fraction_of_one`` exceeds ``FRACTION_THRESHOLD``, and enriches each surviving
species/molecule pair with:

- the full organism taxonomy (one column per rank, domain through species)
- the molecule's 2D SMILES (matching the inchikey_2D already in the posterior)
- the molecule taxonomy (NPClassifier pathway/superclass/class)
- a flag telling whether the species-molecule pair is already known, i.e.
  present in ``lotus.tsv`` (the pairs the tree/model were built from)

Metadata is pulled from ``data/lotus.csv.gz``, deduplicated per species and per
inchikey_2D, and joined onto the filtered posterior. The result is written to
``test_out/acol_Y_posterior_annotated.tsv``.
"""

import polars as pl

Y_POSTERIOR_PATH = "test_out/acol_Y_posterior.txt"
LOTUS_RAW_PATH = "data/lotus.csv.gz"
LOTUS_PAIRS_PATH = "lotus.tsv"
OUT_PATH = "test_out/acol_Y_posterior_annotated.tsv"

FRACTION_THRESHOLD = 0.8

ORGANISM_TAXONOMY_COLUMNS = [
    "organism_taxonomy_01domain",
    "organism_taxonomy_02kingdom",
    "organism_taxonomy_03phylum",
    "organism_taxonomy_04class",
    "organism_taxonomy_05order",
    "organism_taxonomy_06family",
    "organism_taxonomy_07tribe",
    "organism_taxonomy_08genus",
    "organism_taxonomy_09species",
]

MOLECULE_TAXONOMY_COLUMNS = [
    "structure_taxonomy_npclassifier_01pathway",
    "structure_taxonomy_npclassifier_02superclass",
    "structure_taxonomy_npclassifier_03class",
]


def _underscored(col: str) -> pl.Expr:
    # Node names throughout the tree/tsv files have spaces replaced by underscores.
    return pl.col(col).str.replace_all(" ", "_")


posterior = pl.read_csv(Y_POSTERIOR_PATH, separator="\t").filter(
    pl.col("fraction_of_one") > FRACTION_THRESHOLD
)

lotus_raw = pl.read_csv(
    LOTUS_RAW_PATH,
    columns=[
        "structure_inchikey",
        "structure_smiles_2D",
        *MOLECULE_TAXONOMY_COLUMNS,
        *ORGANISM_TAXONOMY_COLUMNS,
    ],
    null_values=["NA"],
    infer_schema_length=10000,
).with_columns(pl.col("structure_inchikey").str.slice(0, 14).alias("inchikey_2D"))

species_meta = (
    lotus_raw.select(ORGANISM_TAXONOMY_COLUMNS)
    .with_columns([_underscored(c).alias(c) for c in ORGANISM_TAXONOMY_COLUMNS])
    .drop_nulls(subset=["organism_taxonomy_09species"])
    .unique(subset=["organism_taxonomy_09species"], keep="first")
    .rename({"organism_taxonomy_09species": "species"})
)

molecule_meta = (
    lotus_raw.select(["inchikey_2D", "structure_smiles_2D", *MOLECULE_TAXONOMY_COLUMNS])
    .with_columns([_underscored(c).alias(c) for c in MOLECULE_TAXONOMY_COLUMNS])
    .drop_nulls(subset=["inchikey_2D"])
    .unique(subset=["inchikey_2D"], keep="first")
    .rename({"inchikey_2D": "molecules", "structure_smiles_2D": "smiles"})
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
            *ORGANISM_TAXONOMY_COLUMNS[:-1],
            "species",
            "molecules",
            "smiles",
            *MOLECULE_TAXONOMY_COLUMNS,
            "already_in_lotus",
        ]
    )
    .sort("fraction_of_one", descending=True)
)

annotated.write_csv(OUT_PATH, separator="\t")
print(f"wrote {annotated.height} rows to {OUT_PATH}")
