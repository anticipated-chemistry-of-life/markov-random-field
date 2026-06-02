import pandas as pd

from .tree import Tree


class ParameterGenerator:
    def __init__(
        self,
        trees: tuple[Tree, Tree],
        tree_name: str,
        mean_log_nu: float,
        var_log_nu: float,
        log_nu: float,
        alpha: float,
        binned_branch_length: int,
    ):
        if len(trees) != 2:
            raise ValueError("Trees must be a tuple of two trees.")
        if trees[0].tree_name == tree_name:
            self._this_tree, self._other_tree = trees[0], trees[1]
        else:
            self._this_tree, self._other_tree = trees[1], trees[0]

        if not (0 < alpha < 1):
            raise ValueError("Alpha must be strictly between 0 and 1.")
        if binned_branch_length < 0:
            raise ValueError("Binned branch length must be non-negative.")

        self._binned_branch_lengths: pd.DataFrame
        self._mean_log_nu: pd.DataFrame
        self._var_log_nu: pd.DataFrame
        self._log_nu: pd.DataFrame
        self._alpha: pd.DataFrame

        self._generate_branch_length_params(binned_branch_length)
        self._generate_mean_log_nu_params(mean_log_nu)
        self._generate_var_log_nu_params(var_log_nu)
        self._generate_log_nu_params(log_nu)
        self._generate_alpha_params(alpha)

    def _generate_branch_length_params(self, binned_branch_length: int) -> None:
        node_names = self._this_tree.get_non_root_node_names_in_cpp_order()
        self._binned_branch_lengths = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_branch_lengths_{name}"
                    for name in node_names
                ],
                "value": [int(binned_branch_length)] * len(node_names),
            }
        )

    def _generate_mean_log_nu_params(self, mean_log_nu: float) -> None:
        self._mean_log_nu = pd.DataFrame(
            {
                "name": [f"{self._this_tree.tree_name}_mean_log_nu"],
                "value": [mean_log_nu],
            }
        )

    def _generate_var_log_nu_params(self, var_log_nu: float) -> None:
        self._var_log_nu = pd.DataFrame(
            {
                "name": [f"{self._this_tree.tree_name}_var_log_nu"],
                "value": [var_log_nu],
            }
        )

    def _generate_log_nu_params(self, log_nu: float) -> None:
        leaf_names = self._other_tree.get_leaf_names_in_cpp_order()
        self._log_nu = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_log_nu_{name}" for name in leaf_names
                ],
                "value": [log_nu] * len(leaf_names),
            }
        )

    def _generate_alpha_params(self, alpha: float) -> None:
        leaf_names = self._other_tree.get_leaf_names_in_cpp_order()
        self._alpha = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_alpha_{name}" for name in leaf_names
                ],
                "value": [alpha] * len(leaf_names),
            }
        )

    def to_dataframe(self) -> pd.DataFrame:
        dfs = [
            self._binned_branch_lengths,
            self._mean_log_nu,
            self._var_log_nu,
            self._log_nu,
            self._alpha,
        ]
        for df in dfs:
            df["value"] = df["value"].astype(str)
        return pd.concat(dfs, axis=0).reset_index(drop=True)
