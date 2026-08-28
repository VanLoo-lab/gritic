import unittest
from types import SimpleNamespace
import numpy as np

from gritic import gritictimer


class MirrorRouteGuardTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.classifier = gritictimer.RouteClassifier(
            4,
            4,
            False,
            'No_WGD',
        )
        cls.source = next(
            route
            for route in cls.classifier.routes.values()
            if route.mirror_route_id != route.route_id
        )
        cls.target = cls.classifier.routes[cls.source.mirror_route_id]

    def test_validate_mirror_rejects_route_without_balanced_components(self):
        unbalanced = next(iter(gritictimer.RouteClassifier(
            2,
            1,
            False,
            'No_WGD',
        ).routes.values()))

        with self.assertRaisesRegex(
            ValueError,
            'Only balanced two-component routes have mirrors',
        ):
            unbalanced._validate_mirror_route(self.source)

    def test_validate_mirror_rejects_wrong_ordered_route(self):
        wrong_route = SimpleNamespace(
            route_id='not-the-expected-mirror',
            mirror_route_id=self.source.route_id,
        )

        with self.assertRaisesRegex(
            ValueError,
            'does not match the expected ordered mirror route ID',
        ):
            self.source._validate_mirror_route(wrong_route)

    def test_validate_mirror_requires_reciprocal_relationship(self):
        nonreciprocal = SimpleNamespace(
            route_id=self.source.mirror_route_id,
            mirror_route_id='not-the-source',
        )

        with self.assertRaisesRegex(
            ValueError,
            'mirror relationships must be reciprocal',
        ):
            self.source._validate_mirror_route(nonreciprocal)

    def test_mirror_node_map_rejects_nonisomorphic_route_structures(self):
        target_tree = self.target.route_tree.main_tree.copy()
        target_tree.remove_node(next(iter(target_tree)))

        with self.assertRaisesRegex(
            ValueError,
            'must be structurally isomorphic',
        ):
            gritictimer._get_mirror_node_map(
                self.source.route_tree.main_tree,
                target_tree,
            )

    def test_transform_mirror_mult_store_requires_balanced_positive_cn(self):
        unbalanced = next(iter(gritictimer.RouteClassifier(
            2,
            1,
            False,
            'No_WGD',
        ).routes.values()))

        with self.assertRaisesRegex(
            ValueError,
            'requires balanced positive copy number',
        ):
            unbalanced._transform_mirror_mult_store(
                np.empty((2, 6)),
                n_subclones=0,
            )

    def test_transform_mirror_mult_store_validates_rank_and_column_count(self):
        expected_columns = 3 * self.target.major_cn + 2
        invalid_stores = (
            np.empty(expected_columns),
            np.empty((3, expected_columns - 1)),
            np.empty((3, expected_columns + 1)),
        )

        for source_store in invalid_stores:
            with self.subTest(shape=source_store.shape):
                with self.assertRaisesRegex(
                    ValueError,
                    'unexpected number of columns',
                ):
                    self.target._transform_mirror_mult_store(
                        source_store,
                        n_subclones=2,
                    )

    def test_transform_mirror_timing_store_validates_source_node_rows(self):
        expected_rows = len(self.source.route_tree.non_phased_node_order)
        invalid_stores = (
            np.empty(3),
            np.empty((expected_rows - 1, 3)),
            np.empty((expected_rows + 1, 3)),
        )

        for source_store in invalid_stores:
            with self.subTest(shape=source_store.shape):
                with self.assertRaisesRegex(
                    ValueError,
                    'rows must match the source route node order',
                ):
                    self.target._transform_mirror_timing_store(
                        source_store,
                        self.source,
                    )

    def test_reuse_requires_every_fitted_mirror_store(self):
        fitted_values = {
            'log_evidence': -1.0,
            'node_timing': np.empty((1, 1)),
            'wgd_timing_store': np.empty(1),
            'mult_store': np.empty((1, 12)),
            'n_events_store': {},
            'density': 1.0,
            'density_high': 1.0,
        }
        original_values = {
            name: getattr(self.source, name, None)
            for name in fitted_values
        }
        try:
            for name, value in fitted_values.items():
                setattr(self.source, name, value)
            self.source.mult_store = None

            with self.assertRaisesRegex(
                ValueError,
                r'fully fitted before reuse; missing: mult_store',
            ):
                self.target._reuse_fitted_unphased_mirror(
                    self.source,
                    n_subclones=0,
                )
        finally:
            for name, value in original_values.items():
                setattr(self.source, name, value)


class MirrorSerializationAndFitGuardTest(unittest.TestCase):
    def test_fit_rejects_balanced_route_whose_mirror_is_missing(self):
        classifier = gritictimer.RouteClassifier(
            4,
            4,
            False,
            'No_WGD',
        )
        source = next(
            route
            for route in classifier.routes.values()
            if route.mirror_route_id != route.route_id
        )
        classifier.routes = {source.route_id: source}

        with self.assertRaisesRegex(
            ValueError,
            'must have its mirror in the route classifier',
        ):
            gritictimer.fit_route_classifiers(
                [(classifier, None, None)],
                None,
            )

    def test_shared_fit_rejects_wgd_status_and_route_id_mismatches(self):
        for mismatch in ('wgd-status', 'route-ids'):
            with self.subTest(mismatch=mismatch):
                first = gritictimer.RouteClassifier(
                    2,
                    1,
                    False,
                    'No_WGD',
                )
                second = gritictimer.RouteClassifier(
                    2,
                    1,
                    False,
                    'No_WGD',
                )
                if mismatch == 'wgd-status':
                    second.wgd_status = True
                else:
                    second.routes = dict(second.routes)
                    second.routes.pop(next(iter(second.routes)))

                with self.assertRaisesRegex(
                    ValueError,
                    'matching copy-number state, WGD status, and route IDs',
                ):
                    gritictimer.fit_route_classifiers(
                        [(first, None, None), (second, None, None)],
                        None,
                    )


if __name__ == '__main__':
    unittest.main()
