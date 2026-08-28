import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import networkx as nx

from gritic import gritictimer, treetools
from gritic.intervaltools import IntervalSpec


def make_plot_tree():
    tree = nx.DiGraph()
    nodes = (
        ('root', False, False, False),
        ('wgd', True, False, False),
        ('gain', False, False, False),
        ('loss', False, True, True),
        ('wgd-left', False, False, True),
        ('wgd-right', False, False, True),
        ('gain-copy', False, False, True),
    )
    for node, wgd_symbol, loss_symbol, terminal_node in nodes:
        tree.add_node(
            node,
            Allele='Major',
            WGD_Symbol=wgd_symbol,
            Loss_Symbol=loss_symbol,
            Terminal_Node=terminal_node,
            Label=node,
        )
    tree.add_edges_from([
        ('root', 'wgd'),
        ('root', 'gain'),
        ('wgd', 'wgd-left'),
        ('wgd', 'wgd-right'),
        ('gain', 'loss'),
        ('gain', 'gain-copy'),
    ])
    return tree


class TreePlotBranchTest(unittest.TestCase):
    def test_convert_to_nx_tree_rejects_noncanonical_allele(self):
        for allele in ('major', 'Other', None):
            with self.subTest(allele=allele):
                with self.assertRaisesRegex(ValueError, 'Major or Minor'):
                    treetools.convert_to_nx_tree(
                        treetools.TreeNode(),
                        allele,
                    )

    def test_plot_tree_uses_stable_loss_wgd_terminal_and_gain_colors(self):
        tree = make_plot_tree()
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.dict(os.environ, {
                'MPLBACKEND': 'Agg',
                'MPLCONFIGDIR': directory,
                'XDG_CACHE_HOME': directory,
            }):
                import matplotlib.pyplot as plt

                with (
                    mock.patch.object(treetools.nx, 'draw') as draw,
                    mock.patch.object(plt, 'savefig') as savefig,
                    mock.patch.object(plt, 'show') as show,
                ):
                    output_path = Path(directory) / 'tree.pdf'
                    treetools.plot_tree(tree, 'Tree', output_path)

        draw.assert_called_once()
        self.assertEqual(draw.call_args.kwargs['node_color'], [
            '#509BCE',
            '#E8BF5E',
            '#509BCE',
            '#4d4d4d',
            '#F7483B',
            '#F7483B',
            '#F7483B',
        ])
        savefig.assert_called_once_with(output_path)
        show.assert_not_called()

    def test_plot_tree_without_output_path_uses_show_and_closes_figure(self):
        tree = make_plot_tree()
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.dict(os.environ, {
                'MPLBACKEND': 'Agg',
                'MPLCONFIGDIR': directory,
                'XDG_CACHE_HOME': directory,
            }):
                import matplotlib.pyplot as plt

                with (
                    mock.patch.object(treetools.nx, 'draw'),
                    mock.patch.object(plt, 'show') as show,
                    mock.patch.object(
                        plt,
                        'close',
                        wraps=plt.close,
                    ) as close,
                ):
                    treetools.plot_tree(tree, 'Tree')

        show.assert_called_once_with()
        close.assert_called_once()


class RouteClassifierPlotTest(unittest.TestCase):
    @staticmethod
    def fake_route(short_id):
        tree = nx.DiGraph()
        tree.add_node(
            short_id,
            Allele='Major',
            WGD_Symbol=False,
            Terminal_Node=True,
        )
        return SimpleNamespace(
            short_id=short_id,
            route_tree=SimpleNamespace(main_tree=tree),
        )

    def test_plot_trees_rejects_classifier_before_fit(self):
        classifier = gritictimer.RouteClassifier(2, 1, False, 'No_WGD')
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "routes haven't been fit"):
                classifier.plot_trees(
                    Path(directory) / 'plots',
                    'Segment',
                    {},
                )

    def test_plot_trees_writes_each_route_with_probability_and_best_fit_title(self):
        classifier = gritictimer.RouteClassifier.__new__(
            gritictimer.RouteClassifier
        )
        first = self.fake_route('alpha')
        second = self.fake_route('beta')
        classifier.routes = {'route-a': first, 'route-b': second}
        classifier.route_probabilities = {'route-a': 0.25, 'route-b': 0.75}
        interval = IntervalSpec(0.8, 'equal-tailed')

        def labels(route, _wgd_info, gain_interval):
            self.assertIs(gain_interval, interval)
            return {route.short_id: f'label-{route.short_id}'}

        with tempfile.TemporaryDirectory() as directory:
            output_directory = Path(directory) / 'plots'
            with mock.patch.object(
                classifier,
                'get_timing_tree_labels',
                side_effect=labels,
            ) as get_labels, mock.patch.object(
                gritictimer.treetools,
                'plot_tree',
            ) as plot_tree:
                classifier.plot_trees(
                    output_directory,
                    'Segment 1',
                    {'WGD_Timing': 0.5},
                    gain_interval=interval,
                )

        self.assertEqual(get_labels.call_count, 2)
        self.assertEqual(plot_tree.call_count, 2)
        first_call, second_call = plot_tree.call_args_list
        self.assertEqual(
            first_call.args[1],
            'Segment 1\nRoute alpha (Probability = 0.25)',
        )
        self.assertEqual(
            second_call.args[1],
            'Segment 1\nRoute beta (Probability = 0.75) - (Best Fit)',
        )
        self.assertEqual(
            first_call.kwargs['output_path'],
            f'{output_directory}/route_alpha.pdf',
        )
        self.assertEqual(
            second_call.kwargs['output_path'],
            f'{output_directory}/route_beta.pdf',
        )
        self.assertEqual(first_call.args[0].nodes['alpha']['Label'], 'label-alpha')
        self.assertEqual(second_call.args[0].nodes['beta']['Label'], 'label-beta')
        self.assertNotIn('Label', first.route_tree.main_tree.nodes['alpha'])
        self.assertNotIn('Label', second.route_tree.main_tree.nodes['beta'])


if __name__ == '__main__':
    unittest.main()
