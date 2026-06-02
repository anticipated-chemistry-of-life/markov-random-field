import click

from .tree import TreeType


TREE_TYPE_CHOICES = click.Choice([t.value for t in TreeType])


def common_tree_options(prefix: str, default_name: str):
    """Return a decorator that adds all CLI options for one tree dimension."""

    def decorator(f):
        opts = [
            click.option(
                f"--{prefix}_tree_type",
                type=TREE_TYPE_CHOICES,
                required=True,
                help=f"Tree topology for the {prefix} dimension.",
            ),
            click.option(
                f"--number_of_{prefix}",
                type=int,
                required=True,
                help=f"Number of nodes in the {prefix} tree.",
            ),
            click.option(
                f"--{prefix}_tree_name",
                type=str,
                default=default_name,
                show_default=True,
                help=f"Name used as tree-node prefix for {prefix}.",
            ),
            click.option(
                f"--{prefix}_branch_length",
                type=int,
                default=10,
                show_default=True,
                help=f"Binned branch length for the {prefix} tree.",
            ),
            click.option(
                f"--{prefix}_mean_log_nu",
                type=float,
                default=0.0,
                show_default=True,
                help=f"Mean of log-nu for the {prefix} tree.",
            ),
            click.option(
                f"--{prefix}_var_log_nu",
                type=float,
                default=0.2,
                show_default=True,
                help=f"Variance of log-nu for the {prefix} tree.",
            ),
            click.option(
                f"--{prefix}_log_nu",
                type=float,
                default=0.5,
                show_default=True,
                help=f"Initial log-nu value for the {prefix} tree.",
            ),
            click.option(
                f"--{prefix}_alpha",
                type=float,
                default=0.5,
                show_default=True,
                help=f"Alpha (Ising coupling strength) for the {prefix} tree.",
            ),
        ]
        for opt in reversed(opts):
            f = opt(f)
        return f

    return decorator
