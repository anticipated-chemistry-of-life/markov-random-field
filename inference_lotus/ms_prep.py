"""Turn the SIRIUS-side scores table into the mass-spec file acol reads.

The scores table has the columns ``species, extract_id, feature_id, inchikey_2d,
calibration_score``. The output is a tab-separated file with the header
``species extract feature molecule calibration_score``, which ``--ms_data`` takes.

Spaces in species and molecule names become underscores, as in the trees
``data_prep.py`` writes. So names match the tree leaves exactly, and acol needs no
naming rule of its own.

The script keeps every row. A feature's reference is its best row, candidate or
not, so acol must see the rows it later skips. The script reads every value as
text, so a calibration score reaches acol exactly as SIRIUS wrote it.
"""

import click
import polars as pl

SCORES_TO_MS_DATA = {
    "species": "species",
    "extract_id": "extract",
    "feature_id": "feature",
    "inchikey_2d": "molecule",
    "calibration_score": "calibration_score",
}


@click.command()
@click.option(
    "--scores",
    type=click.Path(exists=True, dir_okay=False),
    default="feature_molecule_scores.csv",
    show_default=True,
    help="SIRIUS-side scores table, comma-separated.",
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
    ms_data.write_csv(out, separator="\t")
    print(f"wrote {ms_data.height} rows to {out}")


if __name__ == "__main__":
    main()
