import numpy as np
import networkx as nx
import gritic.treetools as treetools


class RouteTree:
    def __init__(self, tree, major_cn, minor_cn, wgd_status):
        self.main_tree = tree

        self.major_cn = major_cn
        self.minor_cn = minor_cn

        major_tree, minor_tree = treetools.split_tree(tree)
        self.major_tree = major_tree
        self.minor_tree = minor_tree

        self.wgd_status = wgd_status
        self.node_attributes = treetools.get_node_attributes(
            self.main_tree, self.wgd_status
        )

        self.non_phased_node_order = self.get_node_order()

        self.timeable_nodes = self.get_timeable_nodes()
        self.wgd_nodes = self.get_wgd_nodes()
        self.sum_constraint_matrix = self.get_sum_constraint_matrix()
        self.wgd_constraint_matrix = self.get_wgd_constraint_matrix()

        self.timing_matrix = self.get_timing_matrix()
        self.unphased_timing_matrix = self.get_unphased_timing_matrix()

    def get_timeable_nodes(self):
        timeable_nodes = []
        for node in self.non_phased_node_order:
            if (
                len(self.main_tree.out_edges(node)) == 2
                and not self.main_tree.nodes[node]["WGD_Symbol"]
            ):
                timeable_nodes.append(node)
        return timeable_nodes

    def get_non_terminal_nodes(self):
        non_terminal_nodes = []
        for node in self.non_phased_node_order:
            if len(self.main_tree.out_edges(node)) > 0:
                non_terminal_nodes.append(node)
        return non_terminal_nodes

    def get_terminal_nodes(self):
        terminal_nodes = []
        for node in self.non_phased_node_order:
            if len(self.main_tree.out_edges(node)) == 0:
                terminal_nodes.append(node)
        return terminal_nodes

    def get_wgd_nodes(self):
        wgd_nodes = []
        for node in self.non_phased_node_order:
            if self.main_tree.nodes[node]["WGD_Symbol"]:
                wgd_nodes.append(node)
        return wgd_nodes

    # all nodes
    def get_node_order(self, allele=None):
        if allele is None:
            return list(nx.dfs_preorder_nodes(self.main_tree))
        if allele == "Major":
            return list(nx.dfs_preorder_nodes(self.major_tree))
        if allele == "Minor":
            return list(nx.dfs_preorder_nodes(self.minor_tree))
        raise ValueError(f"allele must be None, Major or Minor, not {allele}")

    # an mxn matrix
    # m = number of timing nodes
    # n is the number of multiplicity states
    def get_timing_matrix(self):

        n_mults = self.major_cn * 2 + self.minor_cn

        mult_offset = 1

        timing_matrix = np.zeros((len(self.non_phased_node_order), n_mults))

        for i, node in enumerate(self.non_phased_node_order):
            final_mult = self.node_attributes[node]["Multiplicity"]
            timing_matrix[i, final_mult - mult_offset] = 1

            if node in self.major_tree.nodes():
                timing_matrix[i, final_mult - mult_offset + self.major_cn] = 1
            if node in self.minor_tree.nodes():
                timing_matrix[i, final_mult - mult_offset + 2 * self.major_cn] = 1
        return timing_matrix

    def get_unphased_timing_matrix(self):
        n_mults = self.major_cn
        mult_offset = 1

        timing_matrix = np.zeros((len(self.non_phased_node_order), n_mults))

        for i, node in enumerate(self.non_phased_node_order):
            final_mult = self.node_attributes[node]["Multiplicity"]
            timing_matrix[i, final_mult - mult_offset] = 1

        return timing_matrix

    # an mxn matrix
    # m = number of paths through the tree
    # n number of timing nodes
    def get_sum_constraint_matrix(self):
        possible_paths = treetools.get_possible_paths(self.main_tree)
        sum_constraint_matrix = np.zeros(
            (len(possible_paths), len(self.non_phased_node_order))
        )

        for row, path in enumerate(possible_paths):
            for node in path:
                sum_constraint_matrix[row, self.non_phased_node_order.index(node)] = 1

        return sum_constraint_matrix

    def get_wgd_constraint_matrix(self):

        wgd_paths = treetools.get_wgd_paths(self.main_tree)
        if len(wgd_paths) == 0:
            return None

        constraint_matrix = np.zeros((len(wgd_paths), len(self.non_phased_node_order)))

        for row, path in enumerate(wgd_paths):
            for node in path:
                constraint_matrix[row, self.non_phased_node_order.index(node)] = 1
        return constraint_matrix

    def get_combined_constraints(self, wgd_timing):
        combined_constraint_matrix = self.sum_constraint_matrix
        combined_constraints_sum = np.ones(combined_constraint_matrix.shape[0])

        if self.wgd_constraint_matrix is not None:
            combined_constraint_matrix = np.vstack(
                [combined_constraint_matrix, self.wgd_constraint_matrix]
            )
            wgd_constraints_sum = (
                np.ones(self.wgd_constraint_matrix.shape[0]) * wgd_timing
            )
            combined_constraints_sum = np.concatenate(
                [combined_constraints_sum, wgd_constraints_sum]
            )

        return combined_constraint_matrix, combined_constraints_sum

    def get_n_events(self, node_timing, wgd_timing):
        n_events_store = []
        if wgd_timing is None:
            # if no wgd the number of events is the number of gains
            n_events = len(
                [
                    node
                    for node in self.main_tree.nodes
                    if len(list(self.main_tree.successors(node))) == 2
                ]
            )
            return n_events, np.nan, np.nan

        pre_wgd_losses = 0
        plotting_tree = treetools.convert_to_plotting_tree(
            self.main_tree, wgd_timing, node_timing, self.non_phased_node_order
        )
        post_wgd_losses = sum(
            [
                int(plotting_tree.nodes[node]["Loss_Symbol"])
                for node in plotting_tree.nodes
            ]
        )

        n_wgds = sum(
            [
                int(self.main_tree.nodes[node]["WGD_Symbol"])
                for node in self.main_tree.nodes
            ]
        )
        n_gains = len(
            [
                node
                for node in self.main_tree.nodes
                if len(list(self.main_tree.successors(node))) == 2
            ]
        )
        # each wgd is only one event
        n_events = n_gains - (n_wgds - 1) + post_wgd_losses
        # add extra loss not accounted for
        if self.minor_cn == 0:
            n_events += 1
            pre_wgd_losses = 1
        return n_events, pre_wgd_losses, post_wgd_losses
