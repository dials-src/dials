from __future__ import annotations

import collections
import os
import pathlib
import shutil
import subprocess

import pytest

from cctbx import uctbx
from dxtbx.model import ExperimentList
from dxtbx.serialize import load

from dials.array_family import flex


def unit_cells_are_similar(
    uc1, uc2, relative_length_tolerance=0.01, absolute_angle_tolerance=1
):
    # see also uctbx.cpp unit_cell::is_similar_to()
    l1 = uc1.parameters()
    l2 = uc2.parameters()
    for i in range(3):
        if abs(min(l1[i], l2[i]) / max(l1[i], l2[i]) - 1) > relative_length_tolerance:
            return False
    for i in range(3, 6):
        if abs(l1[i] - l2[i]) > absolute_angle_tolerance:
            if abs(l1[i] - (180 - l2[i])) > absolute_angle_tolerance:
                return False
    return True


_indexing_result = collections.namedtuple(
    "indexing", ["indexed_reflections", "experiments", "rmsds"]
)


def run_indexing(
    reflections,
    experiment,
    working_directory,
    extra_args,
    expected_unit_cell,
    expected_rmsds,
    expected_hall_symbol,
    n_expected_lattices=1,
    relative_length_tolerance=0.005,
    absolute_angle_tolerance=0.5,
):
    commands = [shutil.which("dials.index")]
    if isinstance(reflections, list):
        commands.extend(reflections)
    else:
        commands.append(reflections)
    if isinstance(experiment, list):
        commands.extend(experiment)
    else:
        commands.append(experiment)
    commands.extend(extra_args)

    result = subprocess.run(commands, cwd=working_directory, capture_output=True)
    assert not result.returncode and not result.stderr

    out_expts = working_directory / "indexed.expt"
    out_refls = working_directory / "indexed.refl"
    assert out_expts.is_file()
    assert out_refls.is_file()

    experiments_list = load.experiment_list(out_expts, check_format=False)
    assert len(experiments_list.crystals()) == n_expected_lattices
    indexed_reflections = flex.reflection_table.from_file(out_refls)
    indexed_reflections.assert_experiment_identifiers_are_consistent(experiments_list)
    rmsds = None

    for i, experiment in enumerate(experiments_list):
        assert unit_cells_are_similar(
            experiment.crystal.get_unit_cell(),
            expected_unit_cell,
            relative_length_tolerance=relative_length_tolerance,
            absolute_angle_tolerance=absolute_angle_tolerance,
        ), (
            experiment.crystal.get_unit_cell().parameters(),
            expected_unit_cell.parameters(),
        )
        sg = experiment.crystal.get_space_group()
        assert sg.type().hall_symbol() == expected_hall_symbol, (
            sg.type().hall_symbol(),
            expected_hall_symbol,
        )
        reflections = indexed_reflections.select(indexed_reflections["id"] == i)
        mi = reflections["miller_index"]
        assert (mi != (0, 0, 0)).count(False) == 0
        reflections = reflections.select(mi != (0, 0, 0))
        reflections = reflections.select(
            reflections.get_flags(reflections.flags.used_in_refinement)
        )
        assert len(reflections) > 0
        obs_x, obs_y, obs_z = reflections["xyzobs.mm.value"].parts()
        calc_x, calc_y, calc_z = reflections["xyzcal.mm"].parts()
        rmsd_x = flex.mean(flex.pow2(obs_x - calc_x)) ** 0.5
        rmsd_y = flex.mean(flex.pow2(obs_y - calc_y)) ** 0.5
        rmsd_z = flex.mean(flex.pow2(obs_z - calc_z)) ** 0.5
        rmsds = (rmsd_x, rmsd_y, rmsd_z)
        for actual, expected in zip(rmsds, expected_rmsds):
            assert actual <= expected, f"{rmsds} {expected_rmsds}"
        assert experiment.identifier != ""
        expt = ExperimentList()
        expt.append(experiment)
        reflections.assert_experiment_identifiers_are_consistent(expt)

    return _indexing_result(indexed_reflections, experiments_list, rmsds)


def test_index_i04_weak_data_fft3d(dials_data, tmp_path):
    # thaumatin
    data_dir = dials_data("i04_weak_data")
    pickle_path = data_dir / "full.pickle"
    sequence_path = data_dir / "experiments_import.json"
    extra_args = [
        "bin_size_fraction=0.25",
        "image_range=1,20",
        "image_range=250,270",
        "image_range=520,540",
    ]
    expected_unit_cell = uctbx.unit_cell((57.7, 57.7, 149.8, 90, 90, 90))
    expected_rmsds = (0.05, 0.04, 0.0005)
    expected_hall_symbol = " P 1"

    run_indexing(
        pickle_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_index_trypsin_two_lattice_P212121(dials_data, tmp_path):
    data_dir = dials_data("semisynthetic_multilattice", pathlib=True)
    strong_path = data_dir / "ag_strong_1_50.refl"
    sequence_path = data_dir / "ag_imported_1_50.expt"
    extra_args = [
        "indexing.method=real_space_grid_search",
        "reflections_per_degree=10",
        "n_macro_cycles=5",
        "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
        "known_symmetry.space_group=P212121",
        "image_range=0,10",
        "beam.fix=all",
        "detector.fix=all",
        "max_cell=70",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.31, 0.49, 0.002)
    expected_hall_symbol = " P 2ac 2ab"
    n_expected_lattices = 1

    run_indexing(
        strong_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=n_expected_lattices,
        relative_length_tolerance=0.02,
        absolute_angle_tolerance=1,
    )


def test_index_i04_weak_data_fft1d(dials_data, tmp_path):
    # thaumatin
    data_dir = dials_data("i04_weak_data")
    pickle_path = data_dir / "full.pickle"
    sequence_path = data_dir / "experiments_import.json"
    extra_args = [
        "n_macro_cycles=2",
        "indexing.method=fft1d",
        "bin_size_fraction=0.25",
        "image_range=1,20",
        "image_range=250,270",
        "image_range=520,540",
    ]
    expected_unit_cell = uctbx.unit_cell((57.7, 57.7, 149.9, 90, 90, 90))
    expected_rmsds = (0.06, 0.05, 0.0005)
    expected_hall_symbol = " P 1"

    run_indexing(
        pickle_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_index_trypsin_index_assignment_local(dials_data, tmp_path):
    data_dir = dials_data("semisynthetic_multilattice", pathlib=True)
    reflections_path = data_dir / "ag_strong_1_50.refl"
    sequence_path = data_dir / "ag_imported_1_50.expt"
    extra_args = [
        "indexing.method=real_space_grid_search",
        "d_min_start=3",
        "n_macro_cycles=3",
        "known_symmetry.unit_cell=54.3,58.3,66.5,90,90,90",
        "known_symmetry.space_group=P212121",
        "image_range=0,10",
        "beam.fix=all",
        "detector.fix=all",
        "max_lattices=2",
        "index_assignment.method=local",
        "nearest_neighbours=50",
    ]
    expected_unit_cell = uctbx.unit_cell((54.3, 58.3, 66.5, 90, 90, 90))
    expected_rmsds = (0.32, 0.41, 0.004)
    expected_hall_symbol = " P 2ac 2ab"
    n_expected_lattices = 2

    run_indexing(
        reflections_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=n_expected_lattices,
        relative_length_tolerance=0.02,
        absolute_angle_tolerance=1,
    )


def test_index_peak_search_clean(dials_data, tmp_path):
    # test indexing from single image of i04_weak_data
    data_dir = dials_data("i04_weak_data")
    pickle_path = data_dir / "first_image.pickle"
    sequence_path = data_dir / "experiments_import.json"
    extra_args = [
        "indexing.method=fft3d",
        "known_symmetry.space_group=P4",
        "known_symmetry.unit_cell=57.8,57.8,150,90,90,90",
        "peak_search=clean",
        "min_samples=15",
        "n_macro_cycles=4",
        "reciprocal_space_grid.d_min=4",
    ]

    expected_unit_cell = uctbx.unit_cell((57.8, 57.8, 150, 90, 90, 90))
    expected_rmsds = (0.06, 0.07, 0.003)
    expected_hall_symbol = " P 4"

    run_indexing(
        pickle_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.mark.parametrize("specify_unit_cell", [False, True])
def test_index_imosflm_tutorial(dials_data, tmp_path, specify_unit_cell):
    # test on spots derived from imosflm tutorial data:
    # http://www.ccp4.ac.uk/courses/BCA2005/tutorials/dataproc-tutorial.html
    data_dir = dials_data("indexing_test_data")
    reflections_path = data_dir / "imosflm_hg_mar-strong.pickle"
    sequence_path = data_dir / "imosflm_hg_mar-experiments.json"

    unit_cell = uctbx.unit_cell((58.373, 58.373, 155.939, 90, 90, 120))
    hall_symbol = '-R 3 2"'

    extra_args = [
        "bin_size_fraction=0.25",
        'known_symmetry.space_group="Hall: %s"' % hall_symbol.replace('"', '\\"'),
    ]
    if specify_unit_cell:
        extra_args.append(
            'known_symmetry.unit_cell="{} {} {} {} {} {}"'.format(
                *unit_cell.parameters()
            )
        )

    expected_unit_cell = unit_cell
    expected_hall_symbol = hall_symbol
    expected_rmsds = (0.08, 0.11, 0.004)

    run_indexing(
        reflections_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.fixture(scope="session")
def insulin_spotfinding(dials_data, tmp_path_factory):
    """Return experiment and reflection files for 2 images of the insulin dataset"""

    data_dir = dials_data("insulin", pathlib=True)
    tmp_path = tmp_path_factory.mktemp("insulin")

    command = [shutil.which("dials.import")]
    for i, image_path in enumerate(("insulin_1_001.img", "insulin_1_045.img")):
        command.append(data_dir / image_path)

    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    experiment = tmp_path / "imported.expt"
    assert experiment.is_file()

    command = [shutil.which("dials.find_spots"), "nproc=1", experiment]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    reflections = tmp_path / "strong.refl"
    assert reflections.is_file()

    return experiment, reflections


@pytest.mark.parametrize("method", ["fft3d", "fft1d", "real_space_grid_search"])
def test_index_insulin_multi_sequence(insulin_spotfinding, tmp_path, method):
    experiment, reflections = insulin_spotfinding
    expected_unit_cell = uctbx.unit_cell(
        (78.163, 78.163, 78.163, 90.000, 90.000, 90.000)
    )
    expected_hall_symbol = " I 2 2 3"
    expected_rmsds = (0.05, 0.06, 0.01)
    extra_args = [
        'known_symmetry.unit_cell="{} {} {} {} {} {}"'.format(
            *expected_unit_cell.parameters()
        ),
        f'known_symmetry.space_group="Hall: {expected_hall_symbol}"',
        f"indexing.method={method}",
        "treat_single_image_as_still=False",
    ]
    run_indexing(
        reflections,
        experiment,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.fixture(scope="session")
def insulin_spotfinding_stills(dials_data, tmp_path_factory):
    """Return experiment and reflection files for 1 image of the insulin
    dataset treated as still image"""

    data_dir = dials_data("insulin", pathlib=True)
    tmp_path = tmp_path_factory.mktemp("insulin")

    command = [
        shutil.which("dials.import"),
        "convert_sequences_to_stills=True",
        data_dir / "insulin_1_001.img",
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    experiment = tmp_path / "imported.expt"
    assert experiment.is_file()

    command = [shutil.which("dials.find_spots"), "nproc=1", experiment]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    reflections = tmp_path / "strong.refl"
    assert reflections.is_file()

    return experiment, reflections


@pytest.mark.parametrize("method", ["fft3d", "fft1d", "real_space_grid_search"])
def test_index_insulin_force_stills(insulin_spotfinding_stills, tmp_path, method):
    experiment, reflections = insulin_spotfinding_stills
    expected_unit_cell = uctbx.unit_cell(
        (78.092, 78.092, 78.092, 90.000, 90.000, 90.000)
    )
    expected_hall_symbol = " I 2 2 3"
    expected_rmsds = (0.05, 0.06, 0.01)

    extra_args = [
        "stills.indexer=stills",
        'known_symmetry.unit_cell="{} {} {} {} {} {}"'.format(
            *expected_unit_cell.parameters()
        ),
        f'known_symmetry.space_group="Hall: {expected_hall_symbol}"',
        f"indexing.method={method}",
    ]

    run_indexing(
        reflections,
        experiment,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_multiple_experiments(dials_data: pathlib.Path, tmp_path):
    # Test indexing 4 lysozyme still shots in a single dials.index job
    #   - the first image doesn't index
    #   - the last three images do index
    data_dir = dials_data("indexing_test_data")
    reflections_path = data_dir / "i24_lyso_still-strong.pickle"
    experiments_json = data_dir / "i24_lyso_still-imported.json"

    expected_unit_cell = uctbx.unit_cell((38.06, 78.78, 78.91, 90, 90, 90))
    expected_hall_symbol = " P 1"
    expected_rmsds = (0.1, 0.07, 0.0)

    extra_args = [
        "stills.indexer=sequences",
        "joint_indexing=False",
        "outlier.algorithm=sauter_poon",
    ]

    run_indexing(
        reflections_path,
        experiments_json,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=3,
        relative_length_tolerance=0.01,
    )


def test_index_4rotation(dials_regression: pathlib.Path, tmp_path):
    # 1440 images of 1° rotation each
    data_dir = dials_regression / "indexing_test_data" / "4rotation"
    pickle_path = data_dir / "strong.pickle"
    sequence_path = data_dir / "experiments.json"
    extra_args = [
        "max_refine=10",
        "reflections_per_degree=50",
        "known_symmetry.space_group=R3",
        "n_macro_cycles=3",
    ]
    expected_unit_cell = uctbx.unit_cell((48.397, 48.397, 284.767, 90, 90, 120))
    expected_rmsds = (0.06, 0.08, 0.22)
    expected_hall_symbol = " R 3"

    result = run_indexing(
        pickle_path,
        sequence_path,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 276800, len(result.indexed_reflections)


def test_index_small_molecule_multi_sequence_4(dials_data, tmp_path):
    # test for small molecule multi-sequence indexing, 4 sequences with different values
    # of goniometer.fixed_rotation()
    data_dir = dials_data("indexing_test_data")
    reflections_paths = [
        (data_dir / f"multi_sweep-SWEEP{i + 1}-strong.pickle") for i in range(4)
    ]
    sequence_paths = [
        data_dir / f"multi_sweep-SWEEP{i + 1}-experiments.json" for i in range(4)
    ]
    extra_args = ["known_symmetry.space_group=I4", "filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell((7.310, 7.310, 6.820, 90.000, 90.000, 90.000))
    expected_rmsds = (0.10, 0.7, 0.5)
    expected_hall_symbol = " I 4"

    result = run_indexing(
        reflections_paths,
        sequence_paths,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 1250, len(result.indexed_reflections)


def test_index_small_molecule_multi_sequence_3(
    dials_regression: pathlib.Path, tmp_path
):
    # test for small molecule multi-sequence indexing, 3 sequences with different values
    # of goniometer setting rotation (i.e. phi scans)
    data_dir = dials_regression / "dials-191"
    print(data_dir)
    pickle_paths = [
        sorted(data_dir.glob(f"*_SWEEP{i + 1}_strong.pickle"))[0] for i in range(3)
    ]
    sequence_paths = [
        sorted(data_dir.glob(f"*_SWEEP{i + 1}_experiments.json"))[0] for i in range(3)
    ]
    extra_args = ["filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell(
        (9.440, 15.313, 17.126, 90.073, 90.106, 79.248)
    )
    expected_rmsds = (0.32, 0.34, 0.005)
    expected_hall_symbol = " P 1"

    result = run_indexing(
        pickle_paths,
        sequence_paths,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 12000, len(result.indexed_reflections)
    # expect at least indexed 2000 reflections per experiment
    for i in range(3):
        assert (result.indexed_reflections["id"] == i).count(True) > 2000
    n_indexed_run1 = result.indexed_reflections.get_flags(
        result.indexed_reflections.flags.indexed
    ).count(True)
    # reindex with known orientations
    result = run_indexing(
        tmp_path / "indexed.refl",
        tmp_path / "indexed.expt",
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert (
        result.indexed_reflections.get_flags(
            result.indexed_reflections.flags.indexed
        ).count(True)
        > n_indexed_run1
    )


def test_index_small_molecule_ice_max_cell(dials_regression: pathlib.Path, tmp_path):
    # test for small molecule indexing: presence of ice rings makes max-cell
    # estimation tricky
    data_dir = os.path.join(dials_regression, "indexing_test_data", "MXSW-904")
    pickle_path = os.path.join(data_dir, "1_SWEEP1_strong.pickle")
    experiments = os.path.join(data_dir, "1_SWEEP1_experiments.json")
    extra_args = ["filter_ice=False"]
    expected_unit_cell = uctbx.unit_cell((11.72, 11.72, 11.74, 109.08, 109.24, 108.99))
    expected_rmsds = (0.06, 0.05, 0.04)
    expected_hall_symbol = " P 1"

    result = run_indexing(
        pickle_path,
        experiments,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )
    assert len(result.indexed_reflections) > 1300, len(result.indexed_reflections)


def test_refinement_failure_on_max_lattices_a15(dials_data, tmp_path):
    """Problem: Sometimes there is enough data to index, but not enough to
    refine. If this happens in the (N>1)th crystal of max_lattices, then
    all existing solutions are also dropped."""
    lpe4_expt = (
        dials_data("indexing_test_data", pathlib=True)
        / "lattice_failure-lpe4-2-a15.expt"
    )
    lpe4_pickle = (
        dials_data("indexing_test_data", pathlib=True)
        / "lattice_failure-lpe4-2-a15_strong.pickle"
    )

    result = subprocess.run(
        [
            shutil.which("dials.index"),
            lpe4_pickle,
            lpe4_expt,
            "max_lattices=3",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    experiments_list = load.experiment_list(
        tmp_path / "indexed.expt", check_format=False
    )
    assert len(experiments_list) == 2

    # now try to reindex with existing model
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            tmp_path / "indexed.expt",
            lpe4_pickle,
            "max_lattices=2",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    experiments_list = load.experiment_list(
        tmp_path / "indexed.expt", check_format=False
    )
    assert len(experiments_list) == 2


@pytest.mark.parametrize(
    "indexer_type,expected_n_lattices",
    (("sequences", 3), ("stills", 1)),
)
def test_indexers_dont_lose_reflections(
    dials_data, tmp_path, indexer_type, expected_n_lattices
):
    refl = flex.reflection_table.from_file(
        dials_data("cunir_serial_processed") / "strong_1.refl"
    )
    expts = load.experiment_list(
        dials_data("cunir_serial_processed") / "imported_with_ref_5.expt",
        check_format=False,
    )[0:1]
    refl["imageset_id"] = flex.int(refl.size(), 0)  # needed for centroid_px_to_mm
    refl.centroid_px_to_mm(expts)
    refl.map_centroids_to_reciprocal_space(expts)

    from cctbx import sgtbx, uctbx

    from dials.algorithms.indexing.indexer import Indexer
    from dials.command_line.ssx_index import phil_scope as ssx_index_phil_scope

    params = ssx_index_phil_scope.extract()
    params.indexing.known_symmetry.space_group = sgtbx.space_group_info("P 21 3")
    params.indexing.known_symmetry.unit_cell = uctbx.unit_cell(
        (96.4, 96.4, 96.4, 90, 90, 90)
    )
    params.indexing.multiple_lattice_search.max_lattices = 5
    params.indexing.method = "real_space_grid_search"
    params.indexing.stills.indexer = indexer_type
    if indexer_type == "sequences":
        params.refinement.reflections.outlier.algorithm = "null"
    idxr = Indexer.from_parameters(refl, expts, params=params)
    idxr.index()
    assert len(idxr.refined_experiments) == expected_n_lattices
    assert (
        idxr.refined_reflections.size() + idxr.unindexed_reflections.size()
    ) == refl.size()
    refined = idxr.refined_reflections
    assert (refined.get_flags(refined.flags.indexed)).count(True) == refined.size()
    assert (refined["miller_index"] == (0, 0, 0)).count(True) == 0
    unindexed = idxr.unindexed_reflections
    assert unindexed.get_flags(unindexed.flags.indexed).count(True) == 0
    assert (unindexed["miller_index"] != (0, 0, 0)).count(True) == 0


def test_index_multi_lattice_multi_sweep(dials_data, tmp_path):
    loc = dials_data("semisynthetic_multilattice", pathlib=True)
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            loc / "ag_imported_1_50.expt",
            loc / "ag_strong_1_50.refl",
            loc / "bh_imported_1_50.expt",
            loc / "bh_strong_1_50.refl",
            "max_lattices=2",
            "joint_indexing=False",
            "n_macro_cycles=2",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    expts = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    assert len(expts) == 4
    assert len(expts.crystals()) == 4
    refls = flex.reflection_table.from_file(tmp_path / "indexed.refl")
    refls.assert_experiment_identifiers_are_consistent(expts)
    assert set(refls["imageset_id"]) == {0, 1}
    n_indexed_first = refls.get_flags(refls.flags.indexed).count(True)

    # now try to reindex with existing model
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            tmp_path / "indexed.expt",
            tmp_path / "indexed.refl",
            "max_lattices=2",
            "joint_indexing=False",
            "n_macro_cycles=2",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    expts = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    assert len(expts) == 4
    assert len(expts.crystals()) == 4
    refls = flex.reflection_table.from_file(tmp_path / "indexed.refl")
    refls.assert_experiment_identifiers_are_consistent(expts)
    assert set(refls["imageset_id"]) == {0, 1}
    assert refls.get_flags(refls.flags.indexed).count(True) >= n_indexed_first


def test_stills_indexer_multi_lattice_bug_MosaicSauter2014(dials_data, tmp_path):
    """Problem: In stills_indexer, before calling the refine function, the
    experiment list contains a list of dxtbx crystal models (that are not
    MosaicSauter2014 models). The conversion to MosaicSauter2014 is made
    during the refine step when functions from nave_parameters is called.
    If the experiment list contains more than 1 experiment, for eg.
    multiple lattices, only the first crystal gets assigned mosaicity. In
    actuality, all crystal models should be assigned mosaicity. This test
    only compares whether or not all crystal models have been assigned a
    MosaicSauter2014 model."""

    import dxtbx.model
    from dxtbx.model import Crystal
    from dxtbx.model.experiment_list import (
        Experiment,
        ExperimentList,
        ExperimentListFactory,
    )

    from dials.algorithms.indexing.stills_indexer import StillsIndexer
    from dials.array_family import flex
    from dials.command_line.stills_process import (
        phil_scope as stills_process_phil_scope,
    )

    data_dir = dials_data("iterative_cspad_refinement", pathlib=True)
    experiment_data = data_dir / "cspad_refined_experiments_step6_level2_300.json"
    reflection_data = data_dir / "cspad_reflections_step7_300.pickle"

    refl = flex.reflection_table.from_file(reflection_data)
    explist = ExperimentListFactory.from_json_file(experiment_data, check_format=False)[
        0:2
    ]
    reflist = refl.select(refl["id"] < 2)  # Only use the first 2 for convenience
    # Construct crystal models that don't have mosaicity. These A,B,C values are the same
    # as read in from the dials_regression folder
    # Crystal-0
    cs0 = Crystal(explist[0].crystal)
    exp0 = Experiment(
        imageset=explist[0].imageset,
        beam=explist[0].beam,
        detector=explist[0].detector,
        goniometer=None,
        scan=None,
        crystal=cs0,
    )

    # Crystal-1
    cs1 = Crystal(explist[1].crystal)
    exp1 = Experiment(
        imageset=explist[1].imageset,
        beam=explist[1].beam,
        detector=explist[1].detector,
        goniometer=None,
        scan=None,
        crystal=cs1,
    )
    # Construct a new experiment_list that will be passed on for refinement
    unrefined_explist = ExperimentList([exp0, exp1])
    # Get default params from stills_process and construct StillsIndexer, then run refinement
    params = stills_process_phil_scope.extract()
    SI = StillsIndexer(reflist, unrefined_explist, params=params)
    refined_explist, new_reflist = SI.refine(unrefined_explist, reflist)
    # Now check whether the models have mosaicity after stills_indexer refinement
    # Also check that mosaicity values are within expected limits
    for ii, crys in enumerate(refined_explist.crystals()):
        assert isinstance(crys, dxtbx.model.MosaicCrystalSauter2014)
        if ii == 0:
            assert crys.get_domain_size_ang() == pytest.approx(2242.0, rel=0.1)
        if ii == 1:
            assert crys.get_domain_size_ang() == pytest.approx(2689.0, rel=0.1)


def test_pink_indexer(
    dials_data,
    tmp_path,
):
    data_dir = dials_data("cunir_serial_processed", pathlib=True)
    expt_file = data_dir / "imported_with_ref_5.expt"
    refl_file = data_dir / "strong_5.refl"

    command = [shutil.which("dials.split_experiments"), expt_file, refl_file]
    result = subprocess.run(command, cwd=tmp_path)
    assert not result.returncode and not result.stderr

    command = [shutil.which("dials.combine_experiments")]
    for i in range(5):
        command.append(f"split_{i}.expt")
        command.append(f"split_{i}.refl")
    result = subprocess.run(command, cwd=tmp_path)
    assert not result.returncode and not result.stderr

    extra_args = [
        "joint_indexing=False",
        "indexing.method=pink_indexer",
        "min_lattices=5",
        "percent_bandwidth=2",
        'known_symmetry.space_group="P 21 3"',
        "known_symmetry.unit_cell=96.410,96.410,96.410,90.0,90.0,90.0",
    ]

    expected_unit_cell = uctbx.unit_cell((96.41, 96.41, 96.41, 90, 90, 90))
    expected_rmsds = (0.200, 0.200, 0.000)
    expected_hall_symbol = " P 2ac 2ab 3"

    run_indexing(
        "combined.expt",
        "combined.refl",
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=5,
    )


def test_ffbidx(
    dials_data,
    tmp_path,
):
    try:
        import ffbidx  # noqa: F401
    except ModuleNotFoundError:
        pytest.skip("ffbidx not installed")
    try:
        ffbidx.Indexer()
    except RuntimeError:
        pytest.skip("ffbidx installed but not functional on this system")

    data_dir = dials_data("cunir_serial_processed", pathlib=True)
    expt_file = data_dir / "imported_with_ref_5.expt"
    refl_file = data_dir / "strong_5.refl"

    command = [shutil.which("dials.split_experiments"), expt_file, refl_file]
    result = subprocess.run(command, cwd=tmp_path)
    assert not result.returncode and not result.stderr

    command = [shutil.which("dials.combine_experiments")]
    for i in range(5):
        command.append(f"split_{i}.expt")
        command.append(f"split_{i}.refl")
    result = subprocess.run(command, cwd=tmp_path)
    assert not result.returncode and not result.stderr

    extra_args = [
        "joint_indexing=False",
        "indexing.method=ffbidx",
        'known_symmetry.space_group="P 21 3"',
        "known_symmetry.unit_cell=96.410,96.410,96.410,90.0,90.0,90.0",
    ]

    expected_unit_cell = uctbx.unit_cell((96.41, 96.41, 96.41, 90, 90, 90))
    expected_rmsds = (0.200, 0.200, 0.000)
    expected_hall_symbol = " P 2ac 2ab 3"

    run_indexing(
        "combined.expt",
        "combined.refl",
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
        n_expected_lattices=5,
    )


@pytest.mark.parametrize(
    "indexer_type,fix_cell",
    (("sequences", False), ("stills", True)),
)
def test_index_ED_still_low_res_spot_match(
    dials_data, tmp_path, indexer_type, fix_cell
):
    # test indexing from a single simulated lysozyme ED still

    image_path = (
        dials_data("image_examples", pathlib=True) / "simtbx_FormatSMVJHSim_001.img"
    )

    command = [shutil.which("dials.import"), image_path]
    result = subprocess.run(command, cwd=tmp_path)
    assert not result.returncode and not result.stderr

    experiment = tmp_path / "imported.expt"
    assert experiment.is_file()

    command = [shutil.which("dials.find_spots"), "nproc=1", experiment]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True)
    assert not result.returncode and not result.stderr

    reflections = tmp_path / "strong.refl"

    extra_args = [
        "indexing.method=low_res_spot_match",
        "known_symmetry.space_group=P43212",
        "known_symmetry.unit_cell=78.84,78.84,38.29,90,90,90",
        "stills.indexer=" + indexer_type,
        "n_macro_cycles=2",
        "detector.fix_list=Dist",
    ]
    if fix_cell:
        extra_args += ["crystal.fix=cell"]

    expected_unit_cell = uctbx.unit_cell((78.84, 78.84, 38.29, 90, 90, 90))
    expected_rmsds = (0.0065, 0.0065, 0.000)
    expected_hall_symbol = " P 4nw 2abw"

    run_indexing(
        reflections,
        experiment,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


@pytest.mark.parametrize(
    "cell_params",
    [
        (44.47, 52.85, 62.23, 115.14, 101.72, 90.01),
        (52.85, 62.23, 44.47, 101.72, 90.01, 115.14),
    ],
)
def test_unconventional_P1_cell(dials_data, tmp_path, cell_params):
    """
    Indexing in P1 should succeed even if the cell parameters are provided in
    a non-conventional setting
    """
    data_dir = dials_data("mpro_x0305_processed", pathlib=True)
    experiment = data_dir / "imported.expt"
    reflections = data_dir / "strong.refl"

    cell_params_str = ",".join([str(x) for x in cell_params])
    extra_args = [
        "indexing.method=fft3d",
        "known_symmetry.space_group=P1",
        "known_symmetry.unit_cell=" + cell_params_str,
    ]
    expected_unit_cell = uctbx.unit_cell(cell_params)
    expected_rmsds = (1, 1, 1)
    expected_hall_symbol = " P 1"

    run_indexing(
        reflections,
        experiment,
        tmp_path,
        extra_args,
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_real_space_grid_search_no_unit_cell(dials_data, tmp_path):
    data_dir = dials_data("i04_weak_data")
    experiments_json = data_dir / "experiments_import.json"
    pickle_path = data_dir / "full.pickle"
    commands = [
        shutil.which("dials.index"),
        experiments_json,
        pickle_path,
        "indexing.method=real_space_grid_search",
    ]
    result = subprocess.run(commands, cwd=tmp_path, capture_output=True)
    assert result.stderr
    assert (
        result.stderr.strip()
        == b"Target unit cell must be provided for real_space_grid_search"
    )


def test_index_known_orientation(dials_data, tmp_path):
    data_dir = dials_data("vmxi_proteinase_k_sweeps", pathlib=True)
    experiments_json = data_dir / "experiments_0.json"
    reflections = data_dir / "reflections_0.pickle"

    expected_unit_cell = uctbx.unit_cell((68.395, 68.395, 104, 90, 90, 90))
    expected_rmsds = (0.013, 0.012, 0.008)
    expected_hall_symbol = " P 4"

    run_indexing(
        reflections,
        experiments_json,
        tmp_path,
        [],
        expected_unit_cell,
        expected_rmsds,
        expected_hall_symbol,
    )


def test_all_expt_ids_have_expts(dials_data, tmp_path):
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            dials_data("vmxi_thaumatin_grid_index", pathlib=True) / "split_07602.expt",
            dials_data("vmxi_thaumatin_grid_index", pathlib=True) / "split_07602.refl",
            "stills.indexer=sequences",
            "indexing.method=real_space_grid_search",
            "space_group=P4",
            "unit_cell=58,58,150,90,90,90",
            "max_lattices=8",
            "beam.fix=all",
            "detector.fix=all",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "indexed.refl").is_file()

    refl = flex.reflection_table.from_file(tmp_path / "indexed.refl")
    expt = ExperimentList.from_file(tmp_path / "indexed.expt", check_format=False)
    assert (refl["id"] != -1).count(True) == refl.get_flags(refl.flags.indexed).count(
        True
    )
    refl.assert_experiment_identifiers_are_consistent(expt)

    assert flex.max(refl["id"]) + 1 == len(expt)


def test_multi_lattice_multi_sweep_joint(dials_data, tmp_path):
    # this test data is not really multi-lattice, but we can force it to find multiple
    # lattices by setting a very low minimum_angular_separation=0.001
    # A test to demonstrate a fix for https://github.com/dials/dials/issues/1821

    # first check that all is well if we don't find the extra lattice
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            dials_data("l_cysteine_dials_output", pathlib=True) / "indexed.expt",
            dials_data("l_cysteine_dials_output", pathlib=True) / "indexed.refl",
            "max_lattices=2",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "indexed.refl").is_file()

    expts = ExperimentList.from_file(tmp_path / "indexed.expt", check_format=False)
    refls = flex.reflection_table.from_file(tmp_path / "indexed.refl")
    assert len(expts) == 4
    assert len(expts.crystals()) == 1
    refls.assert_experiment_identifiers_are_consistent(expts)

    # now force it to find a second shared lattice
    result = subprocess.run(
        [
            shutil.which("dials.index"),
            dials_data("l_cysteine_dials_output", pathlib=True) / "indexed.expt",
            dials_data("l_cysteine_dials_output", pathlib=True) / "indexed.refl",
            "max_lattices=2",
            "minimum_angular_separation=0.001",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "indexed.refl").is_file()

    expts = ExperimentList.from_file(tmp_path / "indexed.expt", check_format=False)
    refls = flex.reflection_table.from_file(tmp_path / "indexed.refl")
    assert len(expts) == 8
    assert len(expts.crystals()) == 2
    refls.assert_experiment_identifiers_are_consistent(expts)
