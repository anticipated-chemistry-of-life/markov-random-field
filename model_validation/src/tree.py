from enum import Enum
from typing import Union

import networkx as nx
import numpy as np
import pandas as pd


class TreeType(Enum):
    grass = "grass"
    star = "star"
    balanced = "balanced"


class Tree:
    def __init__(self, number_of_nodes: int, tree_type: TreeType, tree_name: str):
        self._tree_name: str = tree_name
        self._generated: bool = False
        self._tree_type: TreeType = tree_type
        self._graph: Union[nx.Graph, nx.DiGraph]
        self._branch_length: float = 0.2
        self._number_of_leaves: int = 0
        self._number_of_roots: int = 0
        self._number_of_internal_nodes: int = 0
        self._leaf_names_cache: list[str] | None = None
        self._non_root_node_names_cache: list[str] | None = None

        match tree_type:
            case TreeType.grass:
                self._create_grass_tree(number_of_nodes)
            case TreeType.star:
                self._create_star_tree(number_of_nodes)
            case TreeType.balanced:
                self._create_balanced_tree(number_of_nodes)

        self._generated = True
        mapping = {i: f"{self.tree_name}_{i}" for i in self._graph.nodes()}
        self._graph = nx.relabel_nodes(self._graph, mapping)

    @property
    def tree_type(self) -> TreeType:
        return self._tree_type

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

        self._number_of_leaves = number_of_nodes // 2
        self._number_of_roots = number_of_nodes // 2
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

    @staticmethod
    def _is_power_of_two(n: int) -> bool:
        return n > 0 and (n & (n - 1)) == 0

    def _create_balanced_tree(self, number_of_nodes: int) -> None:
        if number_of_nodes < 3:
            raise ValueError("Balanced tree must have at least 3 nodes.")
        if not Tree._is_power_of_two(number_of_nodes + 1):
            raise ValueError("Balanced tree must have 2^k - 1 nodes.")

        self._number_of_roots = 1
        self._number_of_leaves = (number_of_nodes + 1) // 2
        self._number_of_internal_nodes = (
            number_of_nodes - self._number_of_leaves - self._number_of_roots
        )

        self._graph = nx.balanced_tree(2, int(np.log2(number_of_nodes + 1)) - 1)

    def to_dataframe(self) -> pd.DataFrame:
        df = nx.to_pandas_edgelist(self._graph)
        df = df[["target", "source"]]
        df["length"] = self._branch_length
        df.columns = pd.Index(["child", "parent", "length"])
        return df

    def get_graph(self) -> Union[nx.Graph, nx.DiGraph]:
        if not self._generated:
            raise ValueError("Graph has not been generated yet.")
        return self._graph

    def _compute_node_ordering(
        self,
    ) -> tuple[list[str], set[str], list[str], list[str]]:
        """Simulate C++ node-order when reading the tree file.

        Returns (all_nodes_ordered, roots_set, leaf_names, non_root_node_names).
        """
        df = self.to_dataframe()
        node_order: list[str] = []
        node_set: set[str] = set()
        is_child: set[str] = set()
        is_parent: set[str] = set()

        for _, row in df.iterrows():
            child, parent = str(row["child"]), str(row["parent"])
            is_child.add(child)
            is_parent.add(parent)
            if parent not in node_set:
                node_order.append(parent)
                node_set.add(parent)
            if child not in node_set:
                node_order.append(child)
                node_set.add(child)

        roots = is_parent - is_child
        leaves = [n for n in node_order if n in is_child and n not in is_parent]
        non_root_nodes = [n for n in node_order if n not in roots]
        return node_order, roots, leaves, non_root_nodes

    def get_leaf_names_in_cpp_order(self) -> list[str]:
        """Leaf names in the order the C++ inference binary sees them."""
        if self._leaf_names_cache is None:
            _, _, leaves, _ = self._compute_node_ordering()
            self._leaf_names_cache = leaves
        return self._leaf_names_cache

    def get_non_root_node_names_in_cpp_order(self) -> list[str]:
        """Names of all non-root nodes (leaves + internal non-roots) in C++ order.

        This is the set used for branch-length parameters.
        """
        if self._non_root_node_names_cache is None:
            _, _, _, non_roots = self._compute_node_ordering()
            self._non_root_node_names_cache = non_roots
        return self._non_root_node_names_cache

    def generate_papers_number(self) -> pd.DataFrame:
        """Sample Poisson-distributed paper counts for each leaf node."""
        if not self._generated:
            raise ValueError("Graph has not been generated yet.")

        if isinstance(self._graph, nx.DiGraph):
            is_leaf = lambda node: (
                self._graph.in_degree(node) == 1 and self._graph.out_degree(node) == 0
            )
        else:
            is_leaf = lambda node: self._graph.degree(node) == 1

        nodes = [n for n in self._graph.nodes() if is_leaf(n)]
        papers = [np.random.poisson(4) for _ in nodes]
        df = pd.DataFrame({self.tree_name: nodes, "number_of_papers": papers})
        return df[df["number_of_papers"] > 0].reset_index(drop=True)
