from __future__ import annotations

import json
import os.path
import pathlib
import shutil
import subprocess

import pytest

from dxtbx.serialize import load

from dials.command_line.ssx_index import run


def test_ssx_index_reference_geometry(dials_data, tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = ssx / "imported_with_ref_5.expt"
    refls = ssx / "strong_5.refl"
    pathlib.Path.mkdir(tmp_path / "nuggets")
    result = subprocess.run(
        [
            shutil.which("dials.ssx_index"),
            expts,
            refls,
            "output.nuggets=nuggets",
            "min_spots=72",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "dials.ssx_index.html").is_file()
    experiments = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    assert len(experiments) == 3  # only 3 out of the 5 get indexed
    for i in [0, 1, 2, 4]:
        assert tmp_path.joinpath(
            f"nuggets/nugget_index_merlin0047_1700{i}.cbf.json"
        ).is_file()
    filtered_json = tmp_path.joinpath("nuggets/nugget_index_filtered_images.json")
    assert filtered_json.is_file()
    with open(filtered_json) as f:
        data = json.load(f)
    assert data["filtered_images"] == [4]


@pytest.mark.parametrize("indexer", ["stills", "sequences"])
def test_ssx_index_no_reference_geometry(dials_data, tmp_path, indexer):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = ssx / "imported_no_ref_5.expt"
    refls = ssx / "strong_5.refl"

    args = [
        shutil.which("dials.ssx_index"),
        expts,
        refls,
        f"stills.indexer={indexer}",
        "-vv",
    ]
    if indexer == "sequences":
        args.append("refinement.reflections.outlier.algorithm=null")
    result = subprocess.run(
        args,
        cwd=tmp_path,
        capture_output=True,
    )

    assert not result.returncode and not result.stderr
    assert (tmp_path / "indexed.refl").is_file()
    assert (tmp_path / "indexed.expt").is_file()
    assert (tmp_path / "dials.ssx_index.html").is_file()
    experiments = load.experiment_list(tmp_path / "indexed.expt", check_format=False)
    if indexer == "stills":
        assert (
            len(experiments) == 3
        )  # only 3 out of the 5 get indexed if no reference geometry
    elif indexer == "sequences":
        assert (
            len(experiments) == 5
        )  # all 5 get indexed, albeit some with questionably high rmsds.


def test_ssx_index_bad_input(dials_data, run_in_tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = str(ssx / "imported_no_ref_5.expt")
    refls = str(ssx / "strong_1.refl")

    run([expts, refls])
    assert os.path.exists("indexed.refl")
    assert os.path.exists("indexed.expt")
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert len(experiments) == 0


def test_ssx_index_input_unit_cell(dials_data, run_in_tmp_path):
    ssx = dials_data("cunir_serial_processed", pathlib=True)
    expts = str(ssx / "imported_with_ref_5.expt")
    refls = str(ssx / "strong_5.refl")

    # invoke the run function
    run(
        [
            expts,
            refls,
            "max_lattices=2",
            "unit_cell=96.4,96.4,96.4,90,90,90",
            "method=fft1d+real_space_grid_search",
        ]
    )

    assert os.path.exists("indexed.refl")
    assert os.path.exists("indexed.expt")
    assert os.path.exists("dials.ssx_index.html")
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert (
        len(experiments) == 5
    )  # only 4 out of the 5 images gets indexed, one has two lattices

    # now if just fft1d method:
    # invoke the run function
    run([expts, refls, "method=fft1d", "unit_cell=96.4,96.4,96.4,90,90,90"])
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert (
        len(experiments) == 2
    )  # only 2 out of the 5 images get indexed without real space grid search

    # test we can run the pink_indexer method through ssx_index also
    run([expts, refls, "method=pink_indexer", "unit_cell=96.4,96.4,96.4,90,90,90"])
    experiments = load.experiment_list("indexed.expt", check_format=False)
    assert len(experiments) == 5
