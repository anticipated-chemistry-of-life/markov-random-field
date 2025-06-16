from enum import Enum
from typing import Union

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
        self._graph: Union[nx.Graph, nx.DiGraph]
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

        self._graph = nx.DiGraph()
        for i in range(self._number_of_leaves):
            self._graph.add_edge(f"root_{i}", f"leaf_{i}")

    def _create_star_tree(self, number_of_nodes: int) -> None:
        if number_of_nodes < 3:
            raise ValueError("Star tree must have at least 3 nodes.")

        self._number_of_leaves = number_of_nodes - 1
        self._number_of_roots = 1
        self._number_of_internal_nodes = 0

        self._graph = nx.star_graph(self._number_of_leaves)

    def isPowerofTwo(self, n):
        if n <= 0:
            return False

        # Calculate log base 2 of n
        logValue = int(np.log2(n))

        # Check if log2(n) is an integer
        # and 2^(logn) = n
        return pow(2, logValue) == n

    def _create_balanced_tree(self, number_of_nodes: int) -> None:
        if number_of_nodes < 3:
            raise ValueError("Balanced tree must have at least 3 nodes.")

        # the number of nodes must be odd and it can only be 3 or 7 or 15
        # this means that the number of nodes must be 2^k - 1

        if not self.isPowerofTwo(number_of_nodes + 1):
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

    def generate_papers_number(self) -> pd.DataFrame:
        """Sample from a Poisson distribution the number of papers for each node if the node is a leaf."""
        if not self._generated:
            raise ValueError("Graph has not been generated yet.")

        papers = []
        nodes = []
        if isinstance(self._graph, nx.DiGraph):
            condition = (
                lambda node: self._graph.in_degree(node) == 1
                and self._graph.out_degree(node) == 0
            )
        else:
            condition = lambda node: self._graph.degree(node) == 1

        for node in self._graph.nodes():
            if condition(node):
                # Leaf node, sample number of papers
                papers.append(np.random.poisson(2))
                nodes.append(node)
        df = pd.DataFrame({self.tree_name: nodes, "number_of_papers": papers})

        # drop all nodes with 0 papers
        df = df[df["number_of_papers"] > 0].reset_index(drop=True)
        return df
