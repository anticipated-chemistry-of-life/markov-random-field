"""Path resolution and true-value layering for `visualize.py`.

These cover the two functions the filenames flow through — `_resolve_inputs` and
`_load_true_values` — rather than the rendering, because the decisions worth
locking are which file gets opened and which value wins when two of them carry
the same parameter. That the plots themselves still come out is checked by
running the command against a real scenario directory.
"""

from __future__ import annotations

import pathlib

import click
import pytest

import visualize as V

COMMANDLINE = click.core.ParameterSource.COMMANDLINE


def _ctx(*explicit: str) -> click.Context:
    """A context in which exactly the named options count as typed by the user."""
    ctx = click.Context(V.main)
    for name in explicit:
        ctx.set_parameter_source(name, COMMANDLINE)
    return ctx


def _write_tsv(
    path: pathlib.Path, header: tuple[str, ...], rows: list[tuple[object, ...]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines += ["\t".join(str(cell) for cell in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def _resolve(scenario: pathlib.Path, ctx: click.Context, **overrides: object) -> object:
    kwargs: dict[str, object] = {
        "run_prefix": V.DEFAULT_RUN_PREFIX,
        "truth_dir": None,
        "true_values": V.DEFAULT_TRUE_VALUES,
        "true_branch_lengths": V.DEFAULT_TRUE_BRANCH_LENGTHS,
        "true_y": V.DEFAULT_TRUE_Y,
        "true_scalars": {},
    }
    kwargs.update(overrides)
    return V._resolve_inputs(ctx, scenario, **kwargs)  # type: ignore[arg-type]


@pytest.fixture
def scenario(tmp_path: pathlib.Path) -> pathlib.Path:
    """A run directory laid out exactly the way the defaults expect.

    Two trees, so the {tree} template has something to expand against, and the
    two default true-value files disagree on `gamma_species` so their layering
    order is observable.
    """
    run = tmp_path / "rung"
    _write_tsv(
        run / "acol_meanVar.txt",
        ("name", "posterior_mean", "posterior_variance"),
        [
            ("gamma_species", 1.0, 0.01),
            ("species_branch_lengths_a", 2.0, 0.5),
            ("molecules_branch_lengths_b", 1.0, 0.25),
            ("weird_parameter", 0.0, 0.0),
        ],
    )
    _write_tsv(run / "acol_trace.txt", ("gamma_species",), [(1.0,), (1.2,)])
    _write_tsv(
        run / "acol_Y_posterior.txt", ("position", "fraction_of_one"), [(0, 0.9)]
    )
    for tree, branch in (("species", "species_branch_lengths_a"), ("molecules", "molecules_branch_lengths_b")):
        _write_tsv(
            run / f"acol_{tree}_statePosteriors.txt",
            ("name", "state_0", "state_1"),
            [(branch, 0.25, 0.75)],
        )
        _write_tsv(
            tmp_path / f"acol_{tree}_simulated.txt", ("name", "value"), [(branch, 1)]
        )

    _write_tsv(
        tmp_path / "acol_input_simulated.txt",
        ("name", "value"),
        [("gamma_species", 9.9), ("species_branch_lengths_a", 1)],
    )
    _write_tsv(
        tmp_path / "acol_simulated.txt", ("name", "value"), [("gamma_species", 1.1)]
    )
    _write_tsv(
        tmp_path / "acol_simulated_Y.txt",
        ("position", "Y_state", "fraction_of_one"),
        [(0, 1, 1.0)],
    )
    return run


# --- 1. the defaults still name what the hardcoded strings named ---


def test_defaults_resolve_to_the_historical_names(scenario):
    inputs = _resolve(scenario, _ctx())
    truth = scenario.parent

    assert inputs.mean_var == scenario / "acol_meanVar.txt"
    assert inputs.trace == scenario / "acol_trace.txt"
    assert inputs.y_posterior == scenario / "acol_Y_posterior.txt"
    assert inputs.true_y == truth / "acol_simulated_Y.txt"
    assert inputs.true_values == (
        truth / "acol_input_simulated.txt",
        truth / "acol_simulated.txt",
    )
    assert inputs.branch_truth == (
        (
            "molecules",
            scenario / "acol_molecules_statePosteriors.txt",
            truth / "acol_molecules_simulated.txt",
        ),
        (
            "species",
            scenario / "acol_species_statePosteriors.txt",
            truth / "acol_species_simulated.txt",
        ),
    )


# --- 2. --run-prefix and --truth-dir redirect their own halves ---


def test_run_prefix_and_truth_dir_redirect_independently(tmp_path):
    run = tmp_path / "rung"
    elsewhere = tmp_path / "elsewhere"
    _write_tsv(run / "zzz_meanVar.txt", ("name", "m", "v"), [("gamma_species", 1.0, 0.1)])
    _write_tsv(
        run / "zzz_species_statePosteriors.txt",
        ("name", "state_0"),
        [("species_branch_lengths_a", 1.0)],
    )
    _write_tsv(elsewhere / "truth.txt", ("name", "value"), [("gamma_species", 1.1)])
    _write_tsv(elsewhere / "y.txt", ("position", "Y_state"), [(0, 1)])
    _write_tsv(
        elsewhere / "acol_species_simulated.txt",
        ("name", "value"),
        [("species_branch_lengths_a", 0)],
    )

    inputs = _resolve(
        run,
        _ctx("run_prefix", "truth_dir", "true_values", "true_y"),
        run_prefix="zzz",
        truth_dir=str(elsewhere),
        true_values=("truth.txt",),
        true_y="y.txt",
    )

    assert inputs.mean_var == run / "zzz_meanVar.txt"
    assert inputs.trace == run / "zzz_trace.txt"
    assert inputs.y_posterior == run / "zzz_Y_posterior.txt"
    assert inputs.true_values == (elsewhere / "truth.txt",)
    assert inputs.true_y == elsewhere / "y.txt"
    # the tree name is read off the file, so it survives a changed prefix
    assert [tree for tree, _, _ in inputs.branch_truth] == ["species"]


def test_absolute_truth_paths_ignore_truth_dir(scenario, tmp_path):
    absolute = tmp_path / "away" / "truth.txt"
    _write_tsv(absolute, ("name", "value"), [("gamma_species", 1.1)])

    inputs = _resolve(
        scenario, _ctx("true_values"), true_values=(str(absolute),)
    )

    assert inputs.true_values == (absolute,)


# --- 3. --true-values replaces the defaults and layers left to right ---


def test_true_values_layer_in_order_with_the_last_file_winning(scenario):
    truth = scenario.parent
    ordered = (truth / "acol_input_simulated.txt", truth / "acol_simulated.txt")

    frame = V._load_true_values(ordered, {}, {"gamma_species"})
    assert frame.set_index("name")["true_value"]["gamma_species"] == 1.1

    reversed_frame = V._load_true_values(ordered[::-1], {}, {"gamma_species"})
    assert reversed_frame.set_index("name")["true_value"]["gamma_species"] == 9.9


def test_a_name_only_the_earlier_file_carries_survives(scenario):
    truth = scenario.parent
    frame = V._load_true_values(
        (truth / "acol_input_simulated.txt", truth / "acol_simulated.txt"),
        {},
        {"species_branch_lengths_a"},
    )

    values = frame.set_index("name")["true_value"]
    assert values["species_branch_lengths_a"] == 1


# --- 4-6. --true-scalar precedence and validation ---


def test_a_scalar_overrides_the_file_value_for_the_same_name(scenario):
    truth = scenario.parent
    frame = V._load_true_values(
        (truth / "acol_simulated.txt",), {"gamma_species": 2.5}, {"gamma_species"}
    )

    assert frame.set_index("name")["true_value"]["gamma_species"] == 2.5


def test_a_scalar_for_an_uninferred_name_is_an_error(scenario):
    truth = scenario.parent
    with pytest.raises(click.ClickException, match="not a parameter"):
        V._load_true_values(
            (truth / "acol_simulated.txt",), {"gamma_typo": 1.1}, {"gamma_species"}
        )


def test_a_scalar_for_an_unplottable_name_is_an_error(scenario):
    truth = scenario.parent
    with pytest.raises(click.ClickException, match="never plotted"):
        V._load_true_values(
            (truth / "acol_simulated.txt",),
            {"weird_parameter": 1.1},
            {"gamma_species", "weird_parameter"},
        )


# --- 7. explicit paths must exist; defaulted ones keep their old behaviour ---


def test_an_explicit_missing_true_y_is_an_error(scenario):
    with pytest.raises(click.ClickException, match="--true-y"):
        _resolve(scenario, _ctx("true_y"), true_y="nope.txt")


def test_a_defaulted_missing_true_y_is_left_for_the_plot_to_skip(scenario):
    (scenario.parent / "acol_simulated_Y.txt").unlink()

    inputs = _resolve(scenario, _ctx())

    assert inputs.true_y == scenario.parent / "acol_simulated_Y.txt"
    assert not inputs.true_y.exists()


def test_an_explicit_missing_true_values_file_is_an_error(scenario):
    with pytest.raises(click.ClickException, match="--true-values"):
        _resolve(scenario, _ctx("true_values"), true_values=("nope.txt",))


def test_defaulted_true_values_that_all_vanish_fail_at_load(scenario):
    truth = scenario.parent
    with pytest.raises(click.ClickException, match="No true-value file found"):
        V._load_true_values((truth / "gone.txt",), {}, set())


def test_a_missing_inference_output_is_an_error(scenario):
    (scenario / "acol_meanVar.txt").unlink()

    with pytest.raises(click.ClickException, match="No inference output"):
        _resolve(scenario, _ctx())


# --- 8-9. the {tree} template ---


def test_a_template_without_the_placeholder_is_shared_by_every_tree(scenario):
    truth = scenario.parent
    _write_tsv(
        truth / "pooled.txt",
        ("name", "value"),
        [("species_branch_lengths_a", 1), ("molecules_branch_lengths_b", 0)],
    )

    inputs = _resolve(
        scenario, _ctx("true_branch_lengths"), true_branch_lengths="pooled.txt"
    )

    assert [tree for tree, _, _ in inputs.branch_truth] == ["molecules", "species"]
    assert {path for _, _, path in inputs.branch_truth} == {truth / "pooled.txt"}


def test_an_explicit_template_missing_one_tree_is_an_error(scenario):
    (scenario.parent / "acol_molecules_simulated.txt").unlink()

    with pytest.raises(click.ClickException, match="--true-branch-lengths"):
        _resolve(
            scenario,
            _ctx("true_branch_lengths"),
            true_branch_lengths="acol_{tree}_simulated.txt",
        )


def test_a_defaulted_template_missing_one_tree_drops_that_tree(scenario):
    (scenario.parent / "acol_molecules_simulated.txt").unlink()

    inputs = _resolve(scenario, _ctx())

    assert [tree for tree, _, _ in inputs.branch_truth] == ["species"]
