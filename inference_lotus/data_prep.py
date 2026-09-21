"""Build the organism/molecule taxonomy trees, paper counts, and species-molecule pairs
that feed the ACOL model.

Both taxonomies are hierarchies (domain -> ... -> species for organisms, pathway ->
... -> individual molecule for structures). By default the model runs at the finest
resolution (species / molecule), but `--taxon-level` and `--molecule-level` pick any
coarser rank instead: the tree is truncated at that rank (it becomes the leaf level),
and paper counts / species-molecule pairs are aggregated up to it. Aggregation falls
out of counting at the coarser column directly -- e.g. every row naming a species also
names its genus, so counting rows per genus already sums the paper counts of every
species sharing that genus.
"""

from enum import Enum

import click
import polars as pl
from downloaders import BaseDownloader
from tqdm.auto import tqdm

# pandas' default NA-token list: polars only treats an empty field as null unless told otherwise
NA_VALUES = [
    "#N/A",
    "#N/A N/A",
    "#NA",
    "-1.#IND",
    "-1.#QNAN",
    "-NaN",
    "-nan",
    "1.#IND",
    "1.#QNAN",
    "<NA>",
    "N/A",
    "NA",
    "NULL",
    "NaN",
    "None",
    "n/a",
    "nan",
    "null",
]


class TaxonLevel(Enum):
    """Organism taxonomy ranks, ordered root (domain) to leaf (species)."""

    DOMAIN = "organism_taxonomy_01domain"
    KINGDOM = "organism_taxonomy_02kingdom"
    PHYLUM = "organism_taxonomy_03phylum"
    CLASS = "organism_taxonomy_04class"
    ORDER = "organism_taxonomy_05order"
    FAMILY = "organism_taxonomy_06family"
    TRIBE = "organism_taxonomy_07tribe"
    GENUS = "organism_taxonomy_08genus"
    SPECIES = "organism_taxonomy_09species"


class MoleculeLevel(Enum):
    """Molecule taxonomy ranks, ordered root (pathway) to leaf (individual molecule)."""

    PATHWAY = "structure_taxonomy_npclassifier_01pathway"
    SUPERCLASS = "structure_taxonomy_npclassifier_02superclass"
    CLASS = "structure_taxonomy_npclassifier_03class"
    MOLECULE = "inchikey_2D"


# Full root-to-leaf column orderings, used to truncate a hierarchy at a chosen rank.
TAXON_COLUMNS = [level.value for level in TaxonLevel]
MOLECULE_COLUMNS = [level.value for level in MoleculeLevel]


def columns_up_to(all_columns: list[str], level_column: str) -> list[str]:
    """Return the root-to-leaf prefix of `all_columns` ending at `level_column`."""
    return all_columns[: all_columns.index(level_column) + 1]


def collect_edges(
    frame: pl.DataFrame, columns: list[str], desc: str
) -> set[tuple[str, str]]:
    """Pair each rank in a row with the next non-null rank to its left."""
    edges = set()
    for row in tqdm(frame.select(columns).iter_rows(), total=frame.height, desc=desc):
        # Drop Nones while keeping their positions for correct order
        non_na_values = [value for value in row if value is not None]

        # Go from right to left, pairing each with the next available non-None on the left
        for i in range(len(non_na_values) - 1, 0, -1):
            child = non_na_values[i]
            parent = non_na_values[i - 1]
            # Some ranks reuse the same label at two adjacent levels (e.g. NPClassifier's
            # "Miscellaneous alkaloids" superclass and class). A node can't be its own
            # parent, so skip the self-loop; a different, real parent from another row
            # (or no parent at all, i.e. a root) still applies.
            if child == parent:
                continue
            edges.add((child, parent))
    return edges


def edges_to_tree(edges: set[tuple[str, str]]) -> pl.DataFrame:
    edge_df = pl.DataFrame(list(edges), schema=["child", "parent"], orient="row")
    # replace all whitespace with underscores
    edge_df = edge_df.with_columns(
        pl.col("child").str.replace_all(" ", "_"),
        pl.col("parent").str.replace_all(" ", "_"),
        pl.lit(0.5).alias("branch_length"),
    )
    return edge_df.unique(subset="child", keep="first")


def paper_counts(
    frame: pl.DataFrame, level_column: str, valid_nodes: list[str]
) -> pl.DataFrame:
    """Count papers per node of `level_column`, restricted to nodes present in the tree.

    Counting rows directly at the chosen rank aggregates finer-grained counts for free:
    every row naming e.g. a species also names its genus, so counting per genus already
    sums the paper counts of every species sharing that genus.
    """
    counts = frame[level_column].drop_nulls().value_counts(sort=True)
    counts = counts.with_columns(pl.col(level_column).str.replace_all(" ", "_"))
    return counts.filter(pl.col(level_column).is_in(valid_nodes))


@click.command()
@click.option(
    "--taxon-level",
    type=click.Choice([level.name.lower() for level in TaxonLevel], case_sensitive=False),
    default=TaxonLevel.SPECIES.name.lower(),
    show_default=True,
    help="Organism taxonomy rank the tree, paper counts, and pairs are built at.",
)
@click.option(
    "--molecule-level",
    type=click.Choice(
        [level.name.lower() for level in MoleculeLevel], case_sensitive=False
    ),
    default=MoleculeLevel.MOLECULE.name.lower(),
    show_default=True,
    help="Molecule taxonomy rank the tree, paper counts, and pairs are built at.",
)
def main(taxon_level: str, molecule_level: str) -> None:
    taxon = TaxonLevel[taxon_level.upper()]
    molecule = MoleculeLevel[molecule_level.upper()]

    _ = BaseDownloader(auto_extract=False).download(
        "https://zenodo.org/records/7534071/files/230106_frozen_metadata.csv.gz",
        "data/lotus.csv.gz",
    )

    df = pl.read_csv(
        "data/lotus.csv.gz", null_values=NA_VALUES, infer_schema_length=None
    )
    # we get the 14 first characters of the InChIKey
    df = df.with_columns(
        pl.col("structure_inchikey").str.slice(0, 14).alias("inchikey_2D")
    )

    taxon_columns = columns_up_to(TAXON_COLUMNS, taxon.value)
    molecule_columns = columns_up_to(MOLECULE_COLUMNS, molecule.value)

    # Collect edges up to the chosen organism rank; that rank becomes the tree's leaf level.
    edges = collect_edges(
        df, taxon_columns, desc=f"Processing organism taxonomy up to {taxon.name.lower()}"
    )
    edge_df_species = edges_to_tree(edges)
    edge_df_species.write_csv("species.tsv", separator="\t")

    # Collect edges up to the chosen molecule rank; that rank becomes the tree's leaf level.
    edges = collect_edges(
        df,
        molecule_columns,
        desc=f"Processing molecule taxonomy up to {molecule.name.lower()}",
    )
    edge_df_molecules = edges_to_tree(edges)
    edge_df_molecules.write_csv("molecules.tsv", separator="\t")

    mol_papers = paper_counts(df, molecule.value, edge_df_molecules["child"].to_list())
    mol_papers.write_csv("molecules_paper_counts.tsv", separator="\t")

    species_papers = paper_counts(df, taxon.value, edge_df_species["child"].to_list())
    species_papers.write_csv("species_paper_counts.tsv", separator="\t")

    lotus_subset = (
        df.select([taxon.value, molecule.value])
        .unique()
        .with_columns(
            pl.col(taxon.value).str.replace_all(" ", "_"),
            pl.col(molecule.value).str.replace_all(" ", "_"),
        )
        .rename({taxon.value: "species", molecule.value: "molecules"})
    )
    # remove all rows where the molecule or the species is not in the edges
    lotus_subset = lotus_subset.filter(
        pl.col("molecules").is_in(edge_df_molecules["child"].to_list())
        & pl.col("species").is_in(edge_df_species["child"].to_list())
    )
    lotus_subset.write_csv("lotus.tsv", separator="\t")


if __name__ == "__main__":
    main()
