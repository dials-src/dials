from __future__ import absolute_import, division, print_function

import pytest

import libtbx
from cctbx import sgtbx

from dials.algorithms.symmetry.cosym._generate_test_data import generate_test_data
from dials.algorithms.symmetry.cosym import phil_scope
from dials.algorithms.symmetry.cosym import analyse_datasets


@pytest.mark.parametrize(
    ("space_group", "dimensions"), [("P2", None), ("P3", None), ("I23", libtbx.Auto)]
)
def test_cosym_analyse_datasets(space_group, dimensions, run_in_tmpdir):
    import matplotlib

    matplotlib.use("Agg")

    datasets, expected_reindexing_ops = generate_test_data(
        space_group=sgtbx.space_group_info(symbol=space_group).group(),
        unit_cell_volume=10000,
        d_min=1.5,
        map_to_p1=True,
        sample_size=20,
    )
    expected_space_group = sgtbx.space_group_info(symbol=space_group).group()

    params = phil_scope.extract()
    params.cluster.n_clusters = len(expected_reindexing_ops)
    params.dimensions = dimensions

    result = analyse_datasets(datasets, params)
    d = result.as_dict()
    assert d["subgroup_scores"][0]["likelihood"] > 0.89
    assert (
        sgtbx.space_group(d["subgroup_scores"][0]["patterson_group"])
        == sgtbx.space_group_info(space_group).group().build_derived_patterson_group()
    )

    space_groups = {}
    reindexing_ops = {}
    for dataset_id in result.reindexing_ops.iterkeys():
        if 0 in result.reindexing_ops[dataset_id]:
            cb_op = result.reindexing_ops[dataset_id][0]
            reindexing_ops.setdefault(cb_op, set())
            reindexing_ops[cb_op].add(dataset_id)
        if dataset_id in result.space_groups:
            space_groups.setdefault(result.space_groups[dataset_id], set())
            space_groups[result.space_groups[dataset_id]].add(dataset_id)

    assert len(reindexing_ops) == len(expected_reindexing_ops)
    assert sorted(reindexing_ops.keys()) == sorted(expected_reindexing_ops.keys())

    for ridx_set in reindexing_ops.values():
        for expected_set in expected_reindexing_ops.values():
            assert (len(ridx_set.symmetric_difference(expected_set)) == 0) or (
                len(ridx_set.intersection(expected_set)) == 0
            )
