import os
import pathlib

import click
import pandas as pd

from src.cli import common_tree_options
from src.params import ParameterGenerator
from src.tree import Tree, TreeType


def _simulate_script(
    species_tree_name: str, molecules_tree_name: str, out_dir: str
) -> str:
    return f"""\
#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."

cd "$ROOT" && xmake
ACOL=$(find "$ROOT/build" -name "acol" -not -path "*/obj/*" -type f | head -1)

cd "$SCRIPT_DIR"
"$ACOL" simulate \\
    --tree_species {species_tree_name}.txt \\
    --tree_molecules {molecules_tree_name}.txt \\
    --species_paper_counts {species_tree_name}_papers.txt \\
    --molecules_paper_counts {molecules_tree_name}_papers.txt \\
    --num_iterations 2000 \\
    --species_var_log_nu acol_input_simulated.txt \\
    --species_mean_log_nu acol_input_simulated.txt \\
    --species_alpha acol_input_simulated.txt \\
    --species_log_nu acol_input_simulated.txt \\
    --species_branch_lengths acol_input_simulated.txt \\
    --molecules_branch_lengths acol_input_simulated.txt \\
    --molecules_var_log_nu acol_input_simulated.txt \\
    --molecules_mean_log_nu acol_input_simulated.txt \\
    --molecules_log_nu acol_input_simulated.txt \\
    --molecules_alpha acol_input_simulated.txt \\
    --numThreads 1 \\
    --write_Y \\
    --write_Z \\
    --write_joint_log_prob_density \\
    --fixedSeed 12345
"""


def _infer_script(
    species_tree_name: str, molecules_tree_name: str, out_dir: str
) -> str:
    return f"""\
#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
ROOT="$SCRIPT_DIR/../.."

cd "$ROOT" && xmake
ACOL=$(find "$ROOT/build" -name "acol" -not -path "*/obj/*" -type f | head -1)

cd "$SCRIPT_DIR"
mkdir -p test_out
"$ACOL" infer \\
    --out ./test_out/acol \\
    --tree_species {species_tree_name}.txt \\
    --tree_molecules {molecules_tree_name}.txt \\
    --lotus acol_simulated_lotus.tsv \\
    --species_paper_counts {species_tree_name}_papers.txt \\
    --molecules_paper_counts {molecules_tree_name}_papers.txt \\
    --iterations 3000 \\
    --numThreads 1 \\
    --writeTrace \\
    --writeBurnin \\
    --molecules_branch_lengths acol_molecules_simulated.txt \\
    --molecules_branch_lengths.update false \\
    --molecules_log_nu acol_input_simulated.txt \\
    --molecules_mean_log_nu acol_input_simulated.txt \\
    --molecules_var_log_nu acol_input_simulated.txt \\
    --molecules_alpha acol_input_simulated.txt \\
    --species_branch_lengths acol_species_simulated.txt \\
    --species_var_log_nu acol_input_simulated.txt \\
    --species_mean_log_nu acol_input_simulated.txt \\
    --species_log_nu acol_input_simulated.txt \\
    --species_alpha acol_input_simulated.txt \\
    --molecules_var_log_nu.update false \\
    --molecules_mean_log_nu.update false \\
    --molecules_log_nu.update false \\
    --molecules_alpha.update false
"""


@click.command()
@common_tree_options("species", "species")
@common_tree_options("molecules", "molecules")
@click.option(
    "--out",
    type=click.Path(),
    default=None,
    help="Output directory (default: s_<type>_<n>_m_<type>_<n>).",
)
def main(
    species_tree_type: str,
    number_of_species: int,
    species_tree_name: str,
    species_branch_length: int,
    species_mean_log_nu: float,
    species_var_log_nu: float,
    species_log_nu: float,
    species_alpha: float,
    molecules_tree_type: str,
    number_of_molecules: int,
    molecules_tree_name: str,
    molecules_branch_length: int,
    molecules_mean_log_nu: float,
    molecules_var_log_nu: float,
    molecules_log_nu: float,
    molecules_alpha: float,
    out: str | None,
) -> None:
    if out is None:
        out = f"s_{species_tree_type}_{number_of_species}_m_{molecules_tree_type}_{number_of_molecules}"

    out_path = pathlib.Path(out)
    out_path.mkdir(parents=True, exist_ok=True)

    species_tree = Tree(
        number_of_nodes=number_of_species,
        tree_type=TreeType(species_tree_type),
        tree_name=species_tree_name,
    )
    molecules_tree = Tree(
        number_of_nodes=number_of_molecules,
        tree_type=TreeType(molecules_tree_type),
        tree_name=molecules_tree_name,
    )

    species_tree.to_dataframe().to_csv(
        out_path / f"{species_tree_name}.txt", sep="\t", index=False
    )
    species_tree.generate_papers_number().to_csv(
        out_path / f"{species_tree_name}_papers.txt", sep="\t", index=False
    )

    molecules_tree.to_dataframe().to_csv(
        out_path / f"{molecules_tree_name}.txt", sep="\t", index=False
    )
    molecules_tree.generate_papers_number().to_csv(
        out_path / f"{molecules_tree_name}_papers.txt", sep="\t", index=False
    )

    trees = (species_tree, molecules_tree)

    species_params = ParameterGenerator(
        trees=trees,
        tree_name=species_tree_name,
        mean_log_nu=species_mean_log_nu,
        var_log_nu=species_var_log_nu,
        log_nu=species_log_nu,
        alpha=species_alpha,
        binned_branch_length=species_branch_length,
    )
    molecules_params = ParameterGenerator(
        trees=trees,
        tree_name=molecules_tree_name,
        mean_log_nu=molecules_mean_log_nu,
        var_log_nu=molecules_var_log_nu,
        log_nu=molecules_log_nu,
        alpha=molecules_alpha,
        binned_branch_length=molecules_branch_length,
    )

    pd.concat([species_params.to_dataframe(), molecules_params.to_dataframe()]).to_csv(
        out_path / "acol_input_simulated.txt", sep="\t", index=False
    )

    (out_path / "simulate.sh").write_text(
        _simulate_script(species_tree_name, molecules_tree_name, out)
    )
    (out_path / "simulate.sh").chmod(0o755)

    (out_path / "test_out").mkdir(exist_ok=True)
    (out_path / "infer_lotus.sh").write_text(
        _infer_script(species_tree_name, molecules_tree_name, out)
    )
    (out_path / "infer_lotus.sh").chmod(0o755)

    click.echo(f"Generated scenario in: {out_path.resolve()}")
    click.echo(f"  Run simulation : bash {out_path}/simulate.sh")
    click.echo(f"  Run inference  : bash {out_path}/infer_lotus.sh")


if __name__ == "__main__":
    main()
