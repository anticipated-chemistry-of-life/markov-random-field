"""Turn the feature-probabilities export into the mass-spec file acol reads.

The export is ``feature_molecule_scores.csv``, the output of feature-probabilities'
``queries/feature_molecule_scores.sql``. Its columns are ``species, extract_id,
feature_id, inchikey_2d, log_correct_density, log_incorrect_density``. The output is a
tab-separated file with the header
``species extract feature molecule log_correct_density log_incorrect_density``, which
``--ms_data`` takes.

Spaces in species and molecule names become underscores, as in the trees
``data_prep.py`` writes. So names match the tree leaves exactly, and acol needs no
naming rule of its own.

The script keeps every row. A feature's unknown weight scales with its number of rows,
candidate or not, so acol must see the rows it later skips. The script reads every value
as text, so a log density reaches acol exactly as feature-probabilities wrote it.

The export leaves both densities empty for a proposed structure the fitted model never
calibrated. The script stops on such a row instead of writing it, because acol cannot
weigh a candidate it has no likelihood ratio for.
"""

import click
import polars as pl

SCORES_TO_MS_DATA = {
    "species": "species",
    "extract_id": "extract",
    "feature_id": "feature",
    "inchikey_2d": "molecule",
    "log_correct_density": "log_correct_density",
    "log_incorrect_density": "log_incorrect_density",
}

DENSITY_COLUMNS = [name for name in SCORES_TO_MS_DATA.values() if name.startswith("log_")]


def missing_density_rows(ms_data: pl.DataFrame) -> pl.DataFrame:
    """The rows with an empty or null density in either column."""
    is_missing = pl.any_horizontal(
        pl.col(column).is_null() | (pl.col(column).str.strip_chars() == "")
        for column in DENSITY_COLUMNS
    )
    return ms_data.filter(is_missing)


@click.command()
@click.option(
    "--scores",
    type=click.Path(exists=True, dir_okay=False),
    default="feature_molecule_scores.csv",
    show_default=True,
    help="feature-probabilities export, comma-separated.",
)
@click.option(
    "--out",
    type=click.Path(dir_okay=False, writable=True),
    default="ms_data.tsv",
    show_default=True,
    help="Mass-spec file to write, tab-separated.",
)
def main(scores: str, out: str) -> None:
    ms_data = (
        pl.read_csv(scores, infer_schema=False)
        .select(list(SCORES_TO_MS_DATA))
        .rename(SCORES_TO_MS_DATA)
        .with_columns(
            pl.col("species").str.replace_all(" ", "_"),
            pl.col("molecule").str.replace_all(" ", "_"),
        )
    )
    missing = missing_density_rows(ms_data)
    if missing.height:
        raise click.ClickException(
            f"{missing.height} of {ms_data.height} rows in {scores} have no log density, "
            "so the fitted model never calibrated them. First rows:\n"
            f"{missing.head(5)}"
        )
    ms_data.write_csv(out, separator="\t")
    print(f"wrote {ms_data.height} rows to {out}")


if __name__ == "__main__":
    main()
