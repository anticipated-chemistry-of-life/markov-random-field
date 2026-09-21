import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    import polars as pl

    return (pl,)


@app.cell
def _(pl):
    df = pl.read_csv("./test_out/acol_Y_posterior.txt", separator="\t")
    trace = pl.read_csv("./test_out/acol_trace.txt", separator="\t")
    return df, trace


@app.cell
def _(trace):
    trace["species_alpha_LFSIKHRAVDEULJ"]
    return


@app.cell
def _(pl):
    lotus = pl.read_csv("lotus.tsv", separator="\t")
    return (lotus,)


@app.cell
def _(lotus, pl):
    lotus.filter(pl.col("molecules") == "LFSIKHRAVDEULJ")
    return


@app.cell
def _(df, lotus):
    posterior = df.join(lotus, on=["species", "molecules"], how="anti")
    return (posterior,)


@app.cell
def _(pl, posterior):
    posterior.filter(pl.col("fraction_of_one") > 0.8).sort(
        by="fraction_of_one", descending=True
    )
    return


if __name__ == "__main__":
    app.run()
