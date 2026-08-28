import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pandas as pd

from gritic import gritictimer
from gritic.tableschemas import GAIN_TIMING_TABLE_COLUMNS


class GainTimingNodeDtypeTest(unittest.TestCase):
    def test_empty_timing_table_does_not_promote_written_nodes_to_float(self):
        classifier = gritictimer.RouteClassifier.__new__(
            gritictimer.RouteClassifier
        )
        classifier._get_output_routes = mock.Mock(
            return_value=[
                (
                    SimpleNamespace(
                        short_id='route-without-gains',
                        route_tree=SimpleNamespace(timeable_nodes=[]),
                    ),
                    1.0,
                )
            ]
        )
        empty_timing_table = classifier.get_timing_table()
        empty_timing_table['Segment_ID'] = 'segment-without-gains'
        empty_timing_table['Sample_ID'] = 'sample'

        timing_table = pd.DataFrame(
            {
                'Sample_ID': ['sample', 'sample'],
                'Segment_ID': ['segment-with-gains', 'segment-with-gains'],
                'Route': ['route', 'route'],
                'Node': [0, 2],
                'Node_Phasing': ['Major', 'Major'],
                'Timing': [0.1, 0.9],
                'Timing_CI_Low': [0.0, 0.8],
                'Timing_CI_High': [0.2, 1.0],
            }
        )

        with tempfile.TemporaryDirectory() as directory:
            table_path = Path(directory) / 'gain_timing_table.tsv'
            gritictimer._initialize_table(
                table_path,
                GAIN_TIMING_TABLE_COLUMNS,
            )
            gritictimer.write_gain_timing_table(timing_table, table_path)
            gritictimer.write_gain_timing_table(
                empty_timing_table,
                table_path,
            )

            written_table = pd.read_csv(table_path, sep='\t')
            node_column = [
                line.split('\t')[3]
                for line in table_path.read_text().splitlines()[1:]
            ]

        self.assertEqual(node_column, ['0', '2'])
        self.assertTrue(pd.api.types.is_integer_dtype(written_table['Node']))
        self.assertEqual(written_table['Node'].tolist(), [0, 2])


if __name__ == '__main__':
    unittest.main()
