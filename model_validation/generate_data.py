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


if __name__ == "__main__":
    main()
