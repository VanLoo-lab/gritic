import unittest

from gritic import tableschemas


class TableSchemaContractTest(unittest.TestCase):
    def test_key_columns_are_ordered_prefixes_of_their_tables(self):
        self.assertEqual(
            tableschemas.ROUTE_KEY_COLUMNS,
            ['Sample_ID', 'Segment_ID', 'Route'],
        )
        self.assertEqual(
            tableschemas.SEGMENT_KEY_COLUMNS,
            ['Sample_ID', 'Segment_ID'],
        )
        self.assertEqual(
            tableschemas.ROUTE_TABLE_COLUMNS[:3],
            tableschemas.ROUTE_KEY_COLUMNS,
        )
        self.assertEqual(
            tableschemas.GAIN_TIMING_TABLE_COLUMNS[:3],
            tableschemas.ROUTE_KEY_COLUMNS,
        )
        self.assertEqual(
            tableschemas.SEGMENT_METADATA_COLUMNS[:2],
            tableschemas.SEGMENT_KEY_COLUMNS,
        )

    def test_segment_metadata_is_shared_verbatim_by_every_segment_level_output(self):
        metadata = tableschemas.SEGMENT_METADATA_COLUMNS
        for name in (
            'GAIN_DRAW_COLUMNS',
            'ROUTE_DRAW_COLUMNS',
            'SUMMARY_COLUMNS',
        ):
            with self.subTest(schema=name):
                schema = getattr(tableschemas, name)
                self.assertEqual(schema[:len(metadata)], metadata)

        route_table = tableschemas.ROUTE_TABLE_COLUMNS
        self.assertEqual(route_table[:3], tableschemas.ROUTE_KEY_COLUMNS)
        self.assertEqual(route_table[3:11], metadata[2:-1])
        self.assertEqual(route_table[-4], 'WGD_Status')

    def test_specialized_schema_tails_have_stable_semantic_order(self):
        self.assertEqual(tableschemas.GAIN_TIMING_TABLE_COLUMNS, [
            'Sample_ID',
            'Segment_ID',
            'Route',
            'Node',
            'Node_Phasing',
            'Timing',
            'Timing_CI_Low',
            'Timing_CI_High',
        ])
        self.assertEqual(
            tableschemas.NODE_PHASING_LABELS,
            ('Major', 'Minor'),
        )
        self.assertEqual(tableschemas.SUMMARY_VALUE_COLUMNS, [
            'Gain_Index',
            'Proportion',
            'Timing_Median',
            'Timing_Low_CI',
            'Timing_High_CI',
            'Pre_WGD_Probability',
            'Post_WGD_Probability',
            'WGD_Timing_Median',
            'WGD_Timing_Low_CI',
            'WGD_Timing_High_CI',
        ])
        self.assertEqual(
            tableschemas.SUMMARY_COLUMNS,
            tableschemas.SEGMENT_METADATA_COLUMNS
            + tableschemas.SUMMARY_VALUE_COLUMNS,
        )

    def test_draw_schemas_distinguish_route_and_gain_specific_fields(self):
        metadata_length = len(tableschemas.SEGMENT_METADATA_COLUMNS)
        self.assertEqual(tableschemas.ROUTE_DRAW_COLUMNS[metadata_length:], [
            'Posterior_Sample_Index',
            'Route',
            'WGD_Timing',
        ])
        self.assertEqual(tableschemas.GAIN_DRAW_COLUMNS[metadata_length:], [
            'Posterior_Sample_Index',
            'Route',
            'Node',
            'Node_Phasing',
            'Gain_Timing',
            'WGD_Timing',
            'Gain_Index',
        ])

    def test_every_concrete_schema_has_unique_column_names(self):
        for name in (
            'SEGMENT_METADATA_COLUMNS',
            'ROUTE_TABLE_COLUMNS',
            'GAIN_TIMING_TABLE_COLUMNS',
            'GAIN_DRAW_COLUMNS',
            'ROUTE_DRAW_COLUMNS',
            'SUMMARY_COLUMNS',
        ):
            with self.subTest(schema=name):
                schema = getattr(tableschemas, name)
                self.assertEqual(
                    len(schema),
                    len(set(schema)),
                    f'{name} contains duplicate columns',
                )


if __name__ == '__main__':
    unittest.main()
