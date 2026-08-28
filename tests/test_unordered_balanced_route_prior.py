import unittest

import numpy as np
import pandas as pd

from gritic import gritictimer, posteriortablegen


class UnorderedBalancedRoutePriorTest(unittest.TestCase):
    @staticmethod
    def make_classifier(
        major_cn,
        minor_cn,
        *,
        unordered_balanced_route_prior=False,
    ):
        return gritictimer.RouteClassifier(
            major_cn,
            minor_cn,
            False,
            'No_WGD',
            unordered_balanced_route_prior=unordered_balanced_route_prior,
        )

    @staticmethod
    def fit_with_route_evidence(classifier, evidence_by_route_id):
        classifier._prepare_fit()
        for route_id, route in classifier.routes.items():
            route.log_evidence = np.log(evidence_by_route_id[route_id])
        classifier._finalize_fit()

    def test_default_keeps_equal_balanced_route_probabilities_uniform(self):
        classifier = self.make_classifier(4, 4)
        self.assertFalse(classifier.unordered_balanced_route_prior)
        self.assertEqual(len(classifier.routes), 4)

        self.fit_with_route_evidence(
            classifier,
            {route_id: 1.0 for route_id in classifier.routes},
        )

        for probability in classifier.route_probabilities.values():
            self.assertAlmostEqual(probability, 1 / 4)

    def test_unordered_prior_halves_each_off_diagonal_orientation(self):
        classifier = self.make_classifier(
            4,
            4,
            unordered_balanced_route_prior=True,
        )
        self.fit_with_route_evidence(
            classifier,
            {route_id: 1.0 for route_id in classifier.routes},
        )

        diagonal_routes = [
            route
            for route in classifier.routes.values()
            if route.mirror_route_id == route.route_id
        ]
        off_diagonal_routes = [
            route
            for route in classifier.routes.values()
            if route.mirror_route_id != route.route_id
        ]
        self.assertEqual(len(diagonal_routes), 2)
        self.assertEqual(len(off_diagonal_routes), 2)

        for route in diagonal_routes:
            self.assertEqual(classifier.get_route_prior_weight(route), 1.0)
            self.assertAlmostEqual(
                classifier.route_probabilities[route.route_id],
                1 / 3,
            )
        for route in off_diagonal_routes:
            self.assertEqual(classifier.get_route_prior_weight(route), 0.5)
            self.assertAlmostEqual(
                classifier.route_probabilities[route.route_id],
                1 / 6,
            )

    def test_unordered_prior_multiplies_unequal_marginal_evidence(self):
        classifier = self.make_classifier(
            4,
            4,
            unordered_balanced_route_prior=True,
        )
        evidence_by_route_id = {
            route_id: float(index)
            for index, route_id in enumerate(
                sorted(classifier.routes),
                start=1,
            )
        }
        self.fit_with_route_evidence(classifier, evidence_by_route_id)

        weighted_evidence = {
            route_id: evidence_by_route_id[route_id]
            * classifier.get_route_prior_weight(route)
            for route_id, route in classifier.routes.items()
        }
        normalizer = sum(weighted_evidence.values())
        for route_id, evidence in weighted_evidence.items():
            self.assertAlmostEqual(
                classifier.route_probabilities[route_id],
                evidence / normalizer,
            )

    def test_penalized_probability_inherits_prior_without_reapplying_it(self):
        classifier = self.make_classifier(
            4,
            4,
            unordered_balanced_route_prior=True,
        )
        self.fit_with_route_evidence(
            classifier,
            {route_id: 1.0 for route_id in classifier.routes},
        )
        route_ids = list(classifier.route_probabilities)
        route_table = pd.DataFrame({
            'Sample_ID': ['sample'] * len(route_ids),
            'Segment_ID': ['segment'] * len(route_ids),
            'Route': route_ids,
            'Probability': [
                classifier.route_probabilities[route_id]
                for route_id in route_ids
            ],
            'Average_N_Events': [2.0] * len(route_ids),
        })

        penalized_table = posteriortablegen.add_penalized_probability(
            route_table
        )

        np.testing.assert_allclose(
            penalized_table['Penalized_Probability'],
            penalized_table['Probability'],
        )

    def test_wgd_decorations_are_part_of_balanced_route_identity(self):
        classifier = gritictimer.RouteClassifier(
            3,
            3,
            True,
            'Default',
            unordered_balanced_route_prior=True,
        )
        diagonal_routes = [
            route
            for route in classifier.routes.values()
            if route.mirror_route_id == route.route_id
        ]
        off_diagonal_routes = [
            route
            for route in classifier.routes.values()
            if route.mirror_route_id != route.route_id
        ]

        self.assertEqual(len(classifier.routes), 9)
        self.assertEqual(len(diagonal_routes), 3)
        self.assertEqual(len(off_diagonal_routes), 6)
        self.assertEqual(
            sum(
                classifier.get_route_prior_weight(route)
                for route in classifier.routes.values()
            ),
            6.0,
        )
        for route in off_diagonal_routes:
            self.assertIn(route.mirror_route_id, classifier.routes)
            self.assertEqual(classifier.get_route_prior_weight(route), 0.5)

    def test_prior_mode_does_not_change_the_ordered_route_set(self):
        default_classifier = self.make_classifier(4, 4)
        unordered_prior_classifier = self.make_classifier(
            4,
            4,
            unordered_balanced_route_prior=True,
        )

        self.assertEqual(
            set(unordered_prior_classifier.routes),
            set(default_classifier.routes),
        )

    def test_unbalanced_and_single_component_routes_are_unaffected(self):
        for major_cn, minor_cn in ((4, 3), (4, 0)):
            with self.subTest(major_cn=major_cn, minor_cn=minor_cn):
                default_classifier = self.make_classifier(
                    major_cn,
                    minor_cn,
                )
                unordered_prior_classifier = self.make_classifier(
                    major_cn,
                    minor_cn,
                    unordered_balanced_route_prior=True,
                )
                evidence_by_route_id = {
                    route_id: float(index)
                    for index, route_id in enumerate(
                        sorted(default_classifier.routes),
                        start=1,
                    )
                }
                self.fit_with_route_evidence(
                    default_classifier,
                    evidence_by_route_id,
                )
                self.fit_with_route_evidence(
                    unordered_prior_classifier,
                    evidence_by_route_id,
                )

                for route in unordered_prior_classifier.routes.values():
                    self.assertEqual(
                        unordered_prior_classifier.get_route_prior_weight(
                            route
                        ),
                        1.0,
                    )
                self.assertEqual(
                    unordered_prior_classifier.route_probabilities,
                    default_classifier.route_probabilities,
                )

    def test_unordered_prior_option_requires_a_boolean(self):
        for invalid_value in (None, 0, 1, 'true', object()):
            with self.subTest(invalid_value=invalid_value):
                with self.assertRaisesRegex(
                    ValueError,
                    'unordered_balanced_route_prior must be a boolean',
                ):
                    self.make_classifier(
                        4,
                        4,
                        unordered_balanced_route_prior=invalid_value,
                    )

        classifier = self.make_classifier(
            4,
            4,
            unordered_balanced_route_prior=np.bool_(True),
        )
        self.assertIs(classifier.unordered_balanced_route_prior, True)


if __name__ == '__main__':
    unittest.main()
