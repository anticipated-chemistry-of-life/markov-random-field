import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="A simple data generator to test our model."
    )
    parser.add_argument(
        "--out",
        type=str,
        help="The directory where the output files will be saved",
        required=True,
    )

    # species tree
    parser.add_argument(
        "--species_tree_type",
        type=str,
        choices=["grass", "star", "balanced"],
        help="The type of the species tree",
        required=True,
    )
    parser.add_argument(
        "--number_of_species",
        type=int,
        help="The number of species in the species tree",
        required=True,
    )

    parser.add_argument(
        "--species_tree_name",
        type=str,
        help="The name of the species tree",
        default="species",
    )

    parser.add_argument(
        "--species_branch_length",
        type=int,
        help="The branch length of the species tree",
        default=10,
    )
    parser.add_argument(
        "--species_mean_log_nu",
        type=float,
        help="The mean log nu of the species tree",
        default=0.0,
    )
    parser.add_argument(
        "--species_var_log_nu",
        type=float,
        help="The variance of log nu of the species tree",
        default=0.2,
    )

    parser.add_argument(
        "--species_log_nu",
        type=float,
        help="The log nu of the species tree",
        default=0.5,
    )
    parser.add_argument(
        "--species_alpha",
        type=float,
        help="The alpha of the species tree",
        default=0.5,
    )

    # now for the molcules tree
    parser.add_argument(
        "--molecules_tree_name",
        type=str,
        help="The name of the molecules tree",
        default="molecules",
    )
    parser.add_argument(
        "--molecules_tree_type",
        type=str,
        choices=["grass", "star", "balanced"],
        help="The type of the molecules tree",
        required=True,
    )
    parser.add_argument(
        "--number_of_molecules",
        type=int,
        help="The number of molecules in the molecules tree",
        required=True,
    )

    parser.add_argument(
        "--molecules_branch_length",
        type=int,
        help="The branch length of the molecules tree",
        default=10,
    )

    parser.add_argument(
        "--molecules_mean_log_nu",
        type=float,
        help="The mean log nu of the molecules tree",
        default=0.0,
    )

    parser.add_argument(
        "--molecules_var_log_nu",
        type=float,
        help="The variance of log nu of the molecules tree",
        default=0.2,
    )

    parser.add_argument(
        "--molecules_log_nu",
        type=float,
        help="The log nu of the molecules tree",
        default=0.5,
    )
    parser.add_argument(
        "--molecules_alpha",
        type=float,
        help="The alpha of the molecules tree",
        default=0.5,
    )

    args = parser.parse_args()

    return args
