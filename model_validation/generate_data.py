import os

import pandas as pd

from cli import parse_args
from params import ParameterGenerator
from tree import Tree, TreeType

TREE_TYPE = {
    "grass": TreeType.grass,
    "star": TreeType.star,
    "balanced": TreeType.balanced,
}


def create_simulate_file(
    species_tree_name: str,
    molecules_tree_name: str,
) -> str:
    return f"""
../../build/acol simulate \\
    --tree_species {species_tree_name}.txt \\
    --tree_molecules {molecules_tree_name}.txt \\
    --num_iterations 10000 \\
    --species_branch_lengths acol_input_simulated.txt \\
    --species_var_log_nu acol_input_simulated.txt \\
    --species_mean_log_nu acol_input_simulated.txt \\
    --species_log_nu acol_input_simulated.txt \\
    --species_alpha acol_input_simulated.txt \\
    --molecules_branch_lengths acol_input_simulated.txt \\
    --molecules_var_log_nu acol_input_simulated.txt \\
    --molecules_mean_log_nu acol_input_simulated.txt \\
    --molecules_log_nu acol_input_simulated.txt \\
    --molecules_alpha acol_input_simulated.txt \\
    --numThreads 1 \\
    --write_Y \\
    --write_Z \\
    --write_joint_log_prob_density \\
    --write_Y_trace \\
    --write_Z_trace \\
    --fixedSeed 42 \\
 """


def create_infer_file(
    species_tree_name: str,
    molecules_tree_name: str,
) -> str:
    return f"""
../../build/acol infer \\
    --out ./test_out/acol \\
    --tree_species {species_tree_name}.txt \\
    --tree_molecules {molecules_tree_name}.txt \\
    --lotus acol_simulated_lotus.tsv \\
    --iterations 20000 \\
    --numThreads 1 \\
    --writeTrace \\
    --writeBurnin \\
    --write_Y_trace \\
    --set_Y acol_simulated_Y.txt \\
    --Y.update false \\
    --write_Z_trace \\
    --Z.update false \\
    --set_species_Z acol_simulated_Z_species.txt \\
    --set_molecules_Z acol_simulated_Z_molecules.txt \\
    --molecules_branch_lengths acol_input_simulated.txt \\
    --molecules_branch_lengths.update false \\
    --molecules_log_nu acol_input_simulated.txt \\
    --molecules_log_nu.update false \\
    --molecules_mean_log_nu acol_input_simulated.txt \\
    --molecules_mean_log_nu.update false \\
    --molecules_var_log_nu acol_input_simulated.txt \\
    --molecules_var_log_nu.update false \\
    --molecules_alpha acol_input_simulated.txt \\
    --molecules_alpha.update false \\
    --species_branch_lengths acol_input_simulated.txt \\
    --species_branch_lengths.update false \\
    --species_var_log_nu acol_input_simulated.txt \\
    --species_var_log_nu.update false \\
    --species_mean_log_nu acol_input_simulated.txt \\
    --species_mean_log_nu.update false \\
    --species_log_nu acol_input_simulated.txt \\
    --species_log_nu.update false \\
    --species_alpha acol_input_simulated.txt \\
    # --species_alpha.update false \\
"""


def main():
    args = parse_args()
    # Check if the output directory exists
    # If not, create it
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    species_tree = Tree(
        number_of_nodes=args.number_of_species,
        tree_type=TREE_TYPE[args.species_tree_type],
        tree_name=args.species_tree_name,
    )

    species_tree.to_dataframe().to_csv(
        os.path.join(args.out, f"{args.species_tree_name}.txt"),
        sep="\t",
        index=False,
    )

    molecules_tree = Tree(
        number_of_nodes=args.number_of_molecules,
        tree_type=TREE_TYPE[args.molecules_tree_type],
        tree_name=args.molecules_tree_name,
    )
    molecules_tree.to_dataframe().to_csv(
        os.path.join(args.out, f"{args.molecules_tree_name}.txt"),
        sep="\t",
        index=False,
    )

    # Generate the parameters for the species tree
    species_params = ParameterGenerator(
        (species_tree, molecules_tree),
        args.species_tree_name,
        args.species_mean_log_nu,
        args.species_var_log_nu,
        args.species_log_nu,
        args.species_alpha,
        args.species_branch_length,
    )

    # Generate the parameters for the molecules tree
    molecules_params = ParameterGenerator(
        (species_tree, molecules_tree),
        args.molecules_tree_name,
        args.molecules_mean_log_nu,
        args.molecules_var_log_nu,
        args.molecules_log_nu,
        args.molecules_alpha,
        args.molecules_branch_length,
    )

    params = pd.concat(
        [
            species_params.to_dataframe(),
            molecules_params.to_dataframe(),
        ]
    )

    params.to_csv(
        os.path.join(args.out, "acol_input_simulated.txt"),
        sep="\t",
        index=False,
    )

    # Write the simulate file
    sim_file = create_simulate_file(
        species_tree_name=args.species_tree_name,
        molecules_tree_name=args.molecules_tree_name,
    )
    with open(os.path.join(args.out, "simulate.sh"), "w") as f:
        f.write(sim_file)

    # now we create the diretory "test_out" for the output of inference
    if not os.path.exists(os.path.join(args.out, "test_out")):
        os.makedirs(os.path.join(args.out, "test_out"))

    # Write the infer file
    inf_file = create_infer_file(
        species_tree_name=args.species_tree_name,
        molecules_tree_name=args.molecules_tree_name,
    )
    with open(os.path.join(args.out, "infer.sh"), "w") as f:
        f.write(inf_file)


if __name__ == "__main__":
    main()
