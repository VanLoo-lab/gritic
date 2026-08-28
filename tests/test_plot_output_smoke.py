import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import networkx as nx

from gritic import treetools
from tests._plot_cache import MPL_CONFIG_DIRECTORY, XDG_CACHE_DIRECTORY


class PlotOutputSmokeTest(unittest.TestCase):
    def test_route_tree_plot_writes_a_nonempty_pdf_in_headless_mode(self):
        tree = nx.DiGraph()
        tree.add_node(
            0,
            Allele='Major',
            WGD_Symbol=False,
            Terminal_Node=False,
            Label='gain',
        )
        tree.add_node(
            1,
            Allele='Major',
            WGD_Symbol=False,
            Terminal_Node=True,
            Label='copy',
        )
        tree.add_node(
            2,
            Allele='Major',
            WGD_Symbol=False,
            Terminal_Node=True,
            Label='copy',
        )
        tree.add_edges_from([(0, 1), (0, 2)])

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / 'tree.pdf'
            with mock.patch.dict(
                os.environ,
                {
                    'MPLBACKEND': 'Agg',
                    'MPLCONFIGDIR': str(MPL_CONFIG_DIRECTORY),
                    'XDG_CACHE_HOME': str(XDG_CACHE_DIRECTORY),
                },
            ):
                treetools.plot_tree(tree, 'GRITIC tree', output_path)

            self.assertTrue(output_path.is_file())
            self.assertGreater(output_path.stat().st_size, 1_000)
            self.assertEqual(output_path.read_bytes()[:5], b'%PDF-')


if __name__ == '__main__':
    unittest.main()
