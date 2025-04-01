from enum import Enum
from typing import Tuple

import networkx as nx
import numpy as np
import pandas as pd


class DataGenerator:
    def __init__(self):
        pass

    def _generate_():
        pass


class TreeType(Enum):
    grass = 1
    star = 2
    balanced = 3


class Tree:
    def __init__(self, number_of_nodes: int, tree_type: TreeType, tree_name: str):
        self._tree_name: str = tree_name
        self._generated: bool = False
        self._tree_type: str = tree_type
        self._graph: nx.Graph
        self._branch_length: float = 0.2

        if tree_type is TreeType.grass:
            if number_of_nodes < 2:
                raise ValueError("Grass tree must have at least 2 nodes.")
            if number_of_nodes % 2 != 0:
                raise ValueError("Grass tree must have an even number of nodes.")

            self._number_of_leaves = int(number_of_nodes / 2)
            self._number_of_roots = int(number_of_nodes / 2)
            self._number_of_internal_nodes = 0

            self._graph = nx.Graph()
            for i in range(self._number_of_leaves):
                self._graph.add_edge(f"root_{i}", f"leaf_{i}")

        elif tree_type is TreeType.star:
            if number_of_nodes < 3:
                raise ValueError("Star tree must have at least 3 nodes.")

            self._number_of_leaves = number_of_nodes - 1
            self._number_of_roots = 1
            self._number_of_internal_nodes = 0

            self._graph = nx.star_graph(self._number_of_leaves)

        elif tree_type is TreeType.balanced:
            if number_of_nodes < 3:
                raise ValueError("Balanced tree must have at least 3 nodes.")

            # the number of nodes must be odd and it can only be 3 or 7 or 15
            # this means that the number of nodes must be 2^k - 1
            if np.log2(number_of_nodes + 1) % 2 != 0:
                raise ValueError("Balanced tree must have 2^k - 1 nodes.")

            self._number_of_roots = 1
            self._number_of_leaves = int((number_of_nodes + 1) / 2)
            self._number_of_internal_nodes = (
                number_of_nodes - self._number_of_leaves - self._number_of_roots
            )

            self._graph = nx.balanced_tree(2, int(np.log2(number_of_nodes + 1)) - 1)

        else:
            raise ValueError("Tree type must be grass, star or balanced.")

        mapping = {i: f"{self.tree_name}_{i}" for i in self._graph.nodes()}
        self._graph = nx.relabel_nodes(self._graph, mapping)

    def tree_type(self) -> str:
        return self._tree_type.name

    def number_of_leaves(self) -> int:
        return self._number_of_leaves

    def number_of_roots(self) -> int:
        return self._number_of_roots

    def number_of_internal_nodes(self) -> int:
        return self._number_of_internal_nodes

    @property
    def tree_name(self) -> str:
        return self._tree_name

    def generated(self) -> bool:
        return self._generated

    def to_dataframe(self) -> pd.DataFrame:
        df = nx.to_pandas_edgelist(self._graph)
        df = df[["target", "source"]]
        df["length"] = self._branch_length
        df.columns = ["child", "parent", "length"]

        return df

    def get_graph(self) -> nx.Graph:
        if not self._generated:
            raise ValueError("Graph has not been generated yet.")
        return self._graph


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
        ) + 1
        self._binned_branch_lengths = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_branch_length_{i+1}"
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
        length_of_dataframe = self._other_tree.number_of_leaves() + 1

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
        length_of_dataframe = self._other_tree.number_of_leaves() + 1

        self._alpha = pd.DataFrame(
            {
                "name": [
                    f"{self._this_tree.tree_name}_alpha_{i+1}"
                    for i in range(length_of_dataframe)
                ],
                "value": [alpha] * length_of_dataframe,
            }
        )

    def params_to_dataframe(self) -> pd.DataFrame:
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
