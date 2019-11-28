from __future__ import absolute_import, division, print_function

import os
import json
import procrunner
import pytest
from cctbx import sgtbx, uctbx
import scitbx.matrix
from dxtbx.serialize import load
from dxtbx.model import Crystal, Experiment, ExperimentList

from dials.array_family import flex
from dials.algorithms.symmetry.cosym._generate_test_data import (
    generate_experiments_reflections,
)
from dials.command_line.symmetry import map_to_minimum_cell


def test_symmetry_laue_only(dials_data, tmpdir):
    """Simple test to check that dials.symmetry completes"""

    result = procrunner.run(
        [
            "dials.symmetry",
            dials_data("l_cysteine_dials_output") / "20_integrated_experiments.json",
            dials_data("l_cysteine_dials_output") / "20_integrated.pickle",
            dials_data("l_cysteine_dials_output") / "25_integrated_experiments.json",
            dials_data("l_cysteine_dials_output") / "25_integrated.pickle",
            "systematic_absences.check=False",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check()
    assert tmpdir.join("symmetrized.expt").check()
    exps = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert str(exps[0].crystal.get_space_group().info()) == "P 2 2 2"


def test_symmetry_basis_changes_for_C2(tmpdir):
    """Test the correctness of change of basis operations in dials.symmetry

    Supply the unit cell of beta-lactamase, which triggers a change of
    basis from input to minimum during symmetry analysis."""
    os.chdir(tmpdir.strpath)
    unit_cell = (53.173, 61.245, 69.292, 90.0, 93.04675, 90.0)
    space_group = sgtbx.space_group_info("C 2").group()
    experiments, reflections, _ = generate_experiments_reflections(
        space_group=space_group,
        unit_cell=unit_cell,
        sample_size=1,
        map_to_minimum=False,
    )
    experiments.as_json("tmp.expt")
    expt_file = tmpdir.join("tmp.expt").strpath
    joint_table = flex.reflection_table()
    for r in reflections:
        joint_table.extend(r)
    joint_table.as_file("tmp.refl")
    refl_file = tmpdir.join("tmp.refl").strpath

    command = ["dials.symmetry", expt_file, refl_file, "json=symmetry.json"]
    result = procrunner.run(command, working_directory=tmpdir.strpath)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check(file=1)
    assert tmpdir.join("symmetrized.expt").check(file=1)

    expts = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    for v, expected in zip(expts[0].crystal.get_unit_cell().parameters(), unit_cell):
        assert v == pytest.approx(expected)

    # Using the change of basis ops from the json output we should be able to
    # reindex the input experiments to match the output experiments
    with tmpdir.join("symmetry.json").open() as f:
        d = json.load(f)
        cs = experiments[0].crystal.get_crystal_symmetry()
        cb_op_inp_min = sgtbx.change_of_basis_op(str(d["cb_op_inp_min"][0]))
        cb_op_min_best = sgtbx.change_of_basis_op(str(d["subgroup_scores"][0]["cb_op"]))
        assert cs.change_basis(cb_op_min_best * cb_op_inp_min).is_similar_symmetry(
            expts[0].crystal.get_crystal_symmetry()
        )


@pytest.mark.parametrize("option", ["", "exclude_images=0:1500:1800"])
def test_symmetry_with_absences(dials_data, tmpdir, option):
    """Simple test to check that dials.symmetry, with absences, completes"""

    cmd = [
        "dials.symmetry",
        dials_data("l_cysteine_dials_output") / "20_integrated_experiments.json",
        dials_data("l_cysteine_dials_output") / "20_integrated.pickle",
        dials_data("l_cysteine_dials_output") / "25_integrated_experiments.json",
        dials_data("l_cysteine_dials_output") / "25_integrated.pickle",
    ]
    if option:
        cmd.append(option)

    result = procrunner.run(cmd, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check()
    assert tmpdir.join("symmetrized.expt").check()
    expts = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert str(expts[0].crystal.get_space_group().info()) == "P 21 21 21"


def test_symmetry_with_laue_group_override(dials_data, tmpdir):
    """Simple test to check that dials.symmetry, with overridden laue group, completes"""

    result = procrunner.run(
        [
            "dials.symmetry",
            "laue_group=P121",
            "change_of_basis_op=-b,-a,-c",
            dials_data("l_cysteine_dials_output") / "20_integrated_experiments.json",
            dials_data("l_cysteine_dials_output") / "20_integrated.pickle",
            dials_data("l_cysteine_dials_output") / "25_integrated_experiments.json",
            dials_data("l_cysteine_dials_output") / "25_integrated.pickle",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("symmetrized.refl").check()
    assert tmpdir.join("symmetrized.expt").check()
    expts = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert str(expts[0].crystal.get_space_group().info()) == "P 1 21 1"
    # Verify that the unit cell has been reindexed correctly
    assert expts[0].crystal.get_unit_cell().parameters() == pytest.approx(
        (8.21578444269, 5.4815363434, 12.1457047712, 90.0, 90.0, 90.0)
    )


def test_symmetry_absences_only(dials_data, tmpdir):
    """Test the command line script with real data. Proteinase K in P41"""
    location = dials_data("vmxi_proteinase_k_sweeps")

    command = ["dials.symmetry", "laue_group=None"]
    command.append(location.join("experiments_0.json").strpath)
    command.append(location.join("reflections_0.pickle").strpath)

    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    assert tmpdir.join("dials-symmetry.html").check()
    assert tmpdir.join("symmetrized.expt").check()
    exps = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert str(exps[0].crystal.get_space_group().info()) == "P 41"

    # Now try with a d_min
    command += ["d_min=4.0"]
    result = procrunner.run(command, working_directory=tmpdir)
    assert not result.returncode and not result.stderr
    exps = load.experiment_list(
        tmpdir.join("symmetrized.expt").strpath, check_format=False
    )
    assert str(exps[0].crystal.get_space_group().info()) == "P 41"


def test_map_to_minimum_cell():
    # Input and expected output
    input_ucs = [
        (39.7413, 183.767, 140.649, 90, 90, 90),
        (40.16, 142.899, 92.4167, 90, 102.48, 90),
        (180.613, 40.1558, 142.737, 90, 90.0174, 90),
    ]
    input_sgs = ["C 2 2 21", "P 1 2 1", "C 1 2 1"]
    input_hkl = [
        [(1, -75, -71), (1, -73, -70), (1, -71, -69)],
        [(14, -37, -36), (-2, -35, -46), (-3, -34, -47)],
        [(-31, -5, -3), (-25, -3, -3), (-42, -8, -2)],
    ]
    expected_ucs = [
        (39.7413, 94.00755450320204, 140.649, 90.0, 90.0, 77.79717980856927),
        (40.16, 92.46399390642911, 142.899, 90.0, 90.0, 77.3882749092846),
        (
            40.1558,
            92.51154528306184,
            142.73699999999997,
            89.9830147351441,
            90.0,
            77.46527404307477,
        ),
    ]
    expected_output_hkl = [
        [(-1, 37, -71), (-1, 36, -70), (-1, 35, -69)],
        [(-14, 22, 37), (2, 48, 35), (3, 50, 34)],
        [(-5, 13, -3), (-3, 11, -3), (-8, 17, -2)],
    ]

    # Setup the input experiments and reflection tables
    expts = ExperimentList()
    reflections = []
    for uc, sg, hkl in zip(input_ucs, input_sgs, input_hkl):
        uc = uctbx.unit_cell(uc)
        sg = sgtbx.space_group_info(sg).group()
        B = scitbx.matrix.sqr(uc.fractionalization_matrix()).transpose()
        expts.append(Experiment(crystal=Crystal(B, space_group=sg, reciprocal=True)))
        refl = flex.reflection_table()
        refl["miller_index"] = flex.miller_index(hkl)
        reflections.append(refl)

    # Actually run the method we are testing
    expts_min, reflections_min, cb_ops = map_to_minimum_cell(
        expts, reflections, max_delta=5
    )

    # Verify that the unit cells have been transformed as expected
    for expt, uc in zip(expts, expected_ucs):
        assert expt.crystal.get_unit_cell().parameters() == pytest.approx(uc, abs=4e-2)

    # Verify that the cb_ops map the input unit cells to the expected minimum unit cells
    for input_uc, expected_uc, cb_op in zip(input_ucs, expected_ucs, cb_ops):
        assert (
            uctbx.unit_cell(input_uc)
            .change_basis(cb_op)
            .is_similar_to(uctbx.unit_cell(expected_uc))
        )

    # Space group should be set to P1
    assert [expt.crystal.get_space_group().type().number() for expt in expts_min] == [
        1,
        1,
        1,
    ]

    # Verify that the reflections have been reindexed as expected
    # Because the exact choice of minimum cell can be platform-dependent,
    # compare the magnitude, but not the sign of the output hkl values
    for refl, expected_hkl in zip(reflections, expected_output_hkl):
        for hkl, e_hkl in zip(refl["miller_index"], expected_hkl):
            assert [abs(h) for h in hkl] == [abs(eh) for eh in e_hkl]


def test_map_to_minimum_cell_1037():
    # See https://github.com/dials/dials/issues/1037

    # Input and expected output
    input_ucs = [
        (
            4.805202948916906,
            12.808064769657364,
            16.544899201125446,
            106.45808502003258,
            90.0065567098825,
            100.77735674275475,
        ),
        (
            4.808011343212577,
            12.821894835790472,
            16.557339561965573,
            106.48431244651402,
            90.0252848479048,
            100.77252933676507,
        ),
        (
            4.8096632137789985,
            12.815648858527567,
            16.55931712239122,
            106.48990701341536,
            90.01703141314147,
            100.80397887485773,
        ),
        (
            4.807294085194974,
            12.822386757910516,
            16.560411742466663,
            106.43185845358086,
            90.02067929544215,
            100.79522302759383,
        ),
    ]
    input_sgs = ["P1"] * 4
    input_hkl = [[], [], [], []]
    expected_ucs = [
        (4.8052, 12.8081, 16.5449, 106.458, 90.0066, 100.777),
        (4.80801, 12.8219, 16.5573, 106.484, 90.0253, 100.773),
        (4.80966, 12.8156, 16.5593, 106.49, 90.017, 100.804),
        (4.80729, 12.8224, 16.5604, 106.432, 90.0207, 100.795),
    ]
    expected_output_hkl = [[], [], [], []]

    # Setup the input experiments and reflection tables
    expts = ExperimentList()
    reflections = []
    for uc, sg, hkl in zip(input_ucs, input_sgs, input_hkl):
        uc = uctbx.unit_cell(uc)
        sg = sgtbx.space_group_info(sg).group()
        B = scitbx.matrix.sqr(uc.fractionalization_matrix()).transpose()
        expts.append(Experiment(crystal=Crystal(B, space_group=sg, reciprocal=True)))
        refl = flex.reflection_table()
        refl["miller_index"] = flex.miller_index(hkl)
        reflections.append(refl)

    # Actually run the method we are testing
    expts_min, reflections_min, cb_ops = map_to_minimum_cell(
        expts, reflections, max_delta=5
    )

    # Verify that the unit cells have been transformed as expected
    for expt, uc in zip(expts, expected_ucs):
        assert expt.crystal.get_unit_cell().parameters() == pytest.approx(uc, abs=4e-2)

    # Verify that the cb_ops map the input unit cells to the expected minimum unit cells
    for input_uc, expected_uc, cb_op in zip(input_ucs, expected_ucs, cb_ops):
        assert (
            uctbx.unit_cell(input_uc)
            .change_basis(cb_op)
            .is_similar_to(uctbx.unit_cell(expected_uc))
        )

    # Space group should be set to P1
    assert [expt.crystal.get_space_group().type().number() for expt in expts_min] == [
        1,
        1,
        1,
        1,
    ]

    # Verify that the reflections have been reindexed as expected
    # Because the exact choice of minimum cell can be platform-dependent,
    # compare the magnitude, but not the sign of the output hkl values
    for refl, expected_hkl in zip(reflections, expected_output_hkl):
        for hkl, e_hkl in zip(refl["miller_index"], expected_hkl):
            assert [abs(h) for h in hkl] == [abs(eh) for eh in e_hkl]
