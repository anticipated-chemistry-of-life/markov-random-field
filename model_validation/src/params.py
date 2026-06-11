import numpy as np
import pandas as pd

from .tree import Tree


class ParameterGenerator:
    def __init__(
        self,
        trees: tuple[Tree, Tree],
        tree_name: str,
        mean_log_nu: float | None,
        var_log_nu: float | None,
        log_nu: float | None,
        alpha: float | None,
        binned_branch_length: int | None,
    ):
        if len(trees) != 2:
            raise ValueError("Trees must be a tuple of two trees.")
        if trees[0].tree_name == tree_name:
            self._this_tree, self._other_tree = trees[0], trees[1]
        else:
            self._this_tree, self._other_tree = trees[1], trees[0]

        if alpha is not None and not (0 < alpha < 1):
            raise ValueError("Alpha must be strictly between 0 and 1.")
        if binned_branch_length is not None and binned_branch_length < 0:
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

    def _generate_branch_length_params(self, binned_branch_length: int | None) -> None:
        node_names = self._this_tree.get_non_root_node_names_in_cpp_order()
        n = len(node_names)
        if binned_branch_length is not None:
            values = [int(binned_branch_length)] * n
        else:
            values = np.random.randint(0, 20, size=n).tolist()
        self._binned_branch_lengths = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_branch_lengths_{name}"
                    for name in node_names
                ],
                "value": values,
            }
        )

    def _generate_mean_log_nu_params(self, mean_log_nu: float | None) -> None:
        value = (
            mean_log_nu if mean_log_nu is not None else float(np.random.normal(0, 0.5))
        )
        self._mean_log_nu = pd.DataFrame(
            {
                "name": [f"{self._this_tree.tree_name}_mean_log_nu"],
                "value": [value],
            }
        )

    def _generate_var_log_nu_params(self, var_log_nu: float | None) -> None:
        value = (
            var_log_nu
            if var_log_nu is not None
            else float(np.random.uniform(0.05, 0.5))
        )
        self._var_log_nu = pd.DataFrame(
            {
                "name": [f"{self._this_tree.tree_name}_var_log_nu"],
                "value": [value],
            }
        )

    def _generate_log_nu_params(self, log_nu: float | None) -> None:
        leaf_names = self._other_tree.get_leaf_names_in_cpp_order()
        n = len(leaf_names)
        if log_nu is not None:
            values = [log_nu] * n
        else:
            values = np.random.normal(0, 1, size=n).tolist()
        self._log_nu = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_log_nu_{name}" for name in leaf_names
                ],
                "value": values,
            }
        )

    def _generate_alpha_params(self, alpha: float | None) -> None:
        leaf_names = self._other_tree.get_leaf_names_in_cpp_order()
        n = len(leaf_names)
        if alpha is not None:
            values = [alpha] * n
        else:
            values = np.random.beta(2, 2, size=n).tolist()
        self._alpha = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_alpha_{name}" for name in leaf_names
                ],
                "value": values,
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
