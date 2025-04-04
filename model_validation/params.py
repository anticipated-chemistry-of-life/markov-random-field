from typing import Tuple

import pandas as pd

from tree import Tree


class ParameterGenerator:
    def __init__(
        self,
        trees: Tuple[Tree, Tree],
        tree_name: str,
        mean_log_nu: float,
        var_log_nu: float,
        log_nu: float,
        alpha: float,
        binned_branch_length: int,
    ):
        if trees is None:
            raise ValueError("Trees cannot be None.")

        self._this_tree = None
        self._other_tree = None
        if len(trees) != 2:
            raise ValueError("Trees must be a tuple of two trees.")

        if trees[0].tree_name == tree_name:
            self._this_tree = trees[0]
            self._other_tree = trees[1]
        else:
            self._this_tree = trees[1]
            self._other_tree = trees[0]

        if alpha <= 0 or alpha >= 1:
            raise ValueError("Alpha must be between 0 and 1.")

        if binned_branch_length < 0:
            raise ValueError("Binned branch length must be greater than 0.")

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

    def _generate_branch_length_params(self, binned_branch_length) -> None:
        length_of_dataframe = (
            self._this_tree.number_of_leaves()
            + self._this_tree.number_of_internal_nodes()
        )
        self._binned_branch_lengths = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_branch_lengths_{i+1}"
                    for i in range(length_of_dataframe)
                ],
                "value": [int(binned_branch_length)] * length_of_dataframe,
            }
        )

    def _generate_mean_log_nu_params(self, mean_log_nu) -> None:
        self._mean_log_nu = pd.DataFrame(
            {
                "name": f"{self._this_tree.tree_name}_mean_log_nu",
                "value": [mean_log_nu],
            }
        )

    def _generate_var_log_nu_params(self, var_log_nu) -> None:
        self._var_log_nu = pd.DataFrame(
            {
                "name": f"{self._this_tree.tree_name}_var_log_nu",
                "value": [var_log_nu],
            }
        )

    def _generate_log_nu_params(self, log_nu) -> None:
        length_of_dataframe = self._other_tree.number_of_leaves()

        self._log_nu = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_log_nu_{i+1}"
                    for i in range(length_of_dataframe)
                ],
                "value": [log_nu] * length_of_dataframe,
            }
        )

    def _generate_alpha_params(self, alpha) -> None:
        length_of_dataframe = self._other_tree.number_of_leaves()

        self._alpha = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_alpha_{i+1}"
                    for i in range(length_of_dataframe)
                ],
                "value": [alpha] * length_of_dataframe,
            }
        )

    def to_dataframe(self) -> pd.DataFrame:
        # since some values are integers and some are floats, we need to convert them to strings
        self._binned_branch_lengths["value"] = self._binned_branch_lengths[
            "value"
        ].astype(str)
        self._mean_log_nu["value"] = self._mean_log_nu["value"].astype(str)
        self._var_log_nu["value"] = self._var_log_nu["value"].astype(str)
        self._log_nu["value"] = self._log_nu["value"].astype(str)
        self._alpha["value"] = self._alpha["value"].astype(str)
        df = pd.concat(
            [
                self._binned_branch_lengths,
                self._mean_log_nu,
                self._var_log_nu,
                self._log_nu,
                self._alpha,
            ],
            axis=0,
        )
        df = df.reset_index(drop=True)
        return df
