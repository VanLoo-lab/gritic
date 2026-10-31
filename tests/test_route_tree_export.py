"""The topology export is an inference output, independent of any renderer."""

import numpy as np
import pandas as pd
import pytest

from gritic import gritictimer, sampletools
from gritic.tableschemas import ROUTE_TREE_TABLE_COLUMNS


def test_conditioned_wgd_and_exact_edges_are_exported():
    classifier = gritictimer.RouteClassifier(3, 1, True, 'Default')
    route = next(route for route in classifier.routes.values()
                 if route.route_tree.wgd_nodes)
    classifier.route_probabilities = {route.route_id: 1.0}
    route.wgd_timing_store = np.array([.35, .4, .45, .5, .55])
    route.node_timing = np.ones((len(route.route_tree.non_phased_node_order), 5))
    for index, node in enumerate(route.route_tree.non_phased_node_order):
        if node in route.route_tree.wgd_nodes:
            route.node_timing[index] = route.wgd_timing_store
        elif node in route.route_tree.timeable_nodes:
            route.node_timing[index] = np.array([.65, .7, .75, .8, .85])
    exported = classifier.get_route_tree_table()
    assert list(exported) == ROUTE_TREE_TABLE_COLUMNS[2:]
    assert set(exported['Route_WGD_Source']) == {'Route_Posterior'}
    np.testing.assert_allclose(exported['Route_WGD_Timing'], .45)
    tree = route.route_tree.main_tree
    for record in exported.itertuples(index=False):
        parents = list(tree.predecessors(record.Node))
        assert (pd.isna(record.Parent_Node) if not parents else record.Parent_Node == parents[0])
        assert record.Node_Phasing == route.route_tree.node_attributes[record.Node]['Phasing']
        if record.WGD_Symbol:
            assert record.Timing == record.Route_WGD_Timing
        elif tree.out_degree(record.Node) == 0:
            assert (record.Timing, record.Timing_CI_Low, record.Timing_CI_High) == (1, 1, 1)


@pytest.mark.parametrize('minor_cn', [0, 1])
def test_uniform_no_gain_wgd_explicitly_uses_sample_posterior(minor_cn):
    classifier = gritictimer.RouteClassifier(1, minor_cn, True, 'Default')
    route = next(iter(classifier.routes.values()))
    classifier.route_probabilities = {route.route_id: 1.0}
    assert route.is_uniform_no_gain
    assert route.wgd_timing_store is None
    info = dict(WGD_Timing=.5, WGD_Timing_CI_Low=.3, WGD_Timing_CI_High=.7)
    exported = classifier.get_route_tree_table(wgd_info=info)
    assert set(exported['Route_WGD_Source']) == {'Sample_Posterior'}
    assert set(exported['Timing']) == {1.0}
    assert set(exported['Route_WGD_Timing']) == {.5}
    assert exported['Parent_Node'].isna().all()


def _make_sample(wgd):
    states = [(2, 1), (3, 1), (1, 1), (1, 0)] if wgd else [(1, 1), (2, 1), (1, 0)]
    starts = [0, 10000, 11000, 12000][:len(states)]
    ends = [10000, 11000, 12000, 13000][:len(states)]
    copy_number = pd.DataFrame({
        'Chromosome': ['1'] * len(states), 'Segment_Start': starts, 'Segment_End': ends,
        'Major_CN': [state[0] for state in states], 'Minor_CN': [state[1] for state in states],
    })
    positions = [start + offset for start in starts for offset in range(10, 22)]
    mutations = pd.DataFrame({
        'Chromosome': ['1'] * len(positions), 'Position': positions,
        'Tumor_Ref_Count': [20] * len(positions), 'Tumor_Alt_Count': [5] * len(positions),
    })
    return sampletools.Sample(mutations, copy_number, None, sample_id='TREE_WGD' if wgd else 'TREE_NO_WGD',
                              purity=.8, merge_cn=False)


@pytest.mark.parametrize('wgd', [False, True])
def test_process_sample_writes_complete_current_topology(tmp_path, wgd):
    sample = _make_sample(wgd)
    directory = tmp_path / sample.sample_id
    gritictimer.process_sample(sample, sample_dir=directory, wgd_count=int(wgd), random_seed=2031)
    topology = pd.read_csv(directory / f'{sample.sample_id}_route_tree_table.tsv', sep='\t')
    routes = pd.read_csv(directory / f'{sample.sample_id}_route_table.tsv', sep='\t')
    assert list(topology) == ROUTE_TREE_TABLE_COLUMNS
    assert not topology.empty
    assert set(zip(topology.Segment_ID, topology.Route)) == set(zip(routes.Segment_ID, routes.Route))
    for (segment_id, route_id), tree in topology.groupby(['Segment_ID', 'Route']):
        metadata = routes.loc[(routes.Segment_ID == segment_id) & (routes.Route == route_id)].iloc[0]
        parent_ids = set(tree.Parent_Node.dropna())
        leaves = tree.loc[~tree.Node.isin(parent_ids)]
        assert (leaves.Node_Phasing == 'Major').sum() == metadata.Major_CN
        assert (leaves.Node_Phasing == 'Minor').sum() == metadata.Minor_CN
        assert (leaves.Timing == 1).all()
        medians = tree.set_index('Node').Timing
        for row in tree.itertuples(index=False):
            if not pd.isna(row.Parent_Node):
                assert medians[row.Parent_Node] <= row.Timing
        if wgd:
            assert np.isfinite(tree.Route_WGD_Timing).all()
            expected_source = 'Sample_Posterior' if metadata.Major_CN == 1 else 'Route_Posterior'
            assert set(tree.Route_WGD_Source) == {expected_source}
        else:
            assert tree.Route_WGD_Timing.isna().all()
            assert tree.Route_WGD_Source.isna().all()
    if wgd:
        assert 'Sample_Posterior' in set(topology.Route_WGD_Source)
        assert topology.WGD_Symbol.any()
