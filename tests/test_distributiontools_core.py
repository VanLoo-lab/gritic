import unittest

from gritic.distributiontools import get_ids_with_maximum_overlap


class MaximumIntervalOverlapTest(unittest.TestCase):
    def test_empty_interval_store_is_rejected(self):
        with self.assertRaisesRegex(
            ValueError,
            'At least one timing interval is required',
        ):
            get_ids_with_maximum_overlap({}, {})

    def test_single_interval_is_its_own_maximum(self):
        result = get_ids_with_maximum_overlap(
            {'only': (0.25, 0.75)},
            {'only': 17},
        )

        self.assertEqual(result, (['only'], 17, 0.25))

    def test_maximum_is_weighted_by_segment_width_not_segment_count(self):
        intervals = {
            'long-lived': (0.0, 10.0),
            'first-short': (1.0, 2.0),
            'second-short': (3.0, 4.0),
        }
        widths = {
            'long-lived': 1,
            'first-short': 5,
            'second-short': 6,
        }

        overlapping, overlap_width, timing = get_ids_with_maximum_overlap(
            intervals,
            widths,
        )

        self.assertEqual(overlapping, ['long-lived', 'second-short'])
        self.assertEqual(overlap_width, 7)
        self.assertEqual(timing, 3.0)

    def test_closed_interval_endpoints_count_as_overlapping(self):
        overlapping, overlap_width, timing = get_ids_with_maximum_overlap(
            {
                'ending': (0.0, 1.0),
                'starting': (1.0, 2.0),
            },
            {'ending': 2, 'starting': 3},
        )

        self.assertEqual(overlapping, ['ending', 'starting'])
        self.assertEqual(overlap_width, 5)
        self.assertEqual(timing, 1.0)

    def test_disjoint_intervals_select_the_heaviest_segment(self):
        result = get_ids_with_maximum_overlap(
            {
                'early': (-2.0, -1.0),
                'middle': (0.0, 1.0),
                'late': (2.0, 3.0),
            },
            {'early': 2, 'middle': 11, 'late': 7},
        )

        self.assertEqual(result, (['middle'], 11, 0.0))

    def test_nested_intervals_accumulate_all_active_widths(self):
        result = get_ids_with_maximum_overlap(
            {
                'outer': (0.0, 8.0),
                'middle': (1.0, 7.0),
                'inner': (2.0, 6.0),
            },
            {'outer': 2.5, 'middle': 3.5, 'inner': 4.0},
        )

        self.assertEqual(result, (['outer', 'middle', 'inner'], 10.0, 2.0))

    def test_equal_maxima_keep_the_first_encountered_overlap(self):
        result = get_ids_with_maximum_overlap(
            {
                'persistent': (0.0, 5.0),
                'first': (1.0, 2.0),
                'second': (3.0, 4.0),
            },
            {'persistent': 1, 'first': 3, 'second': 3},
        )

        self.assertEqual(result, (['persistent', 'first'], 4, 1.0))

    def test_missing_width_for_an_encountered_segment_is_not_silently_ignored(self):
        with self.assertRaises(KeyError):
            get_ids_with_maximum_overlap(
                {'a': (0.0, 2.0), 'b': (1.0, 3.0)},
                {'a': 1},
            )


if __name__ == '__main__':
    unittest.main()
