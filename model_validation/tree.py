from enum import Enum

import networkx as nx
import numpy as np
import pandas as pd


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
            self._create_grass_tree(number_of_nodes)
            self._generated = True

        elif tree_type is TreeType.star:
            self._create_star_tree(number_of_nodes)
            self._generated = True

        elif tree_type is TreeType.balanced:
            self._create_balanced_tree(number_of_nodes)
            self._generated = True

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

    def _create_grass_tree(self, number_of_nodes: int) -> None:
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

    def _create_star_tree(self, number_of_nodes: int) -> None:
        if number_of_nodes < 3:
            raise ValueError("Star tree must have at least 3 nodes.")

        self._number_of_leaves = number_of_nodes - 1
        self._number_of_roots = 1
        self._number_of_internal_nodes = 0

        self._graph = nx.star_graph(self._number_of_leaves)

    def _create_balanced_tree(self, number_of_nodes: int) -> None:
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
