from __future__ import annotations

import shutil
import subprocess

import iotbx.mtz
from cctbx import uctbx
from libtbx.test_utils import approx_equal


def test(dials_data, tmp_path):
    g = sorted(f for f in dials_data("x4wide", pathlib=True).glob("*"))
    assert len(g) == 90

    commands = [
        [shutil.which("dials.import")] + g,
        [shutil.which("dials.slice_sequence"), "imported.expt", "image_range=80,90"],
        [shutil.which("dials.find_spots"), "imported_80_90.expt", "nproc=1"],
        [
            shutil.which("dials.index"),
            "imported_80_90.expt",
            "strong.refl",
            "space_group=P41212",
        ],
        [
            shutil.which("dials.refine"),
            "indexed.expt",
            "indexed.refl",
            "scan_varying=True",
        ],
        [shutil.which("dials.integrate"), "refined.expt", "indexed.refl", "nproc=1"],
        [
            shutil.which("dials.export"),
            "refined.expt",
            "integrated.refl",
            "partiality_threshold=0.99",
        ],
    ]

    for cmd in commands:
        # print cmd
        result = subprocess.run(cmd, cwd=tmp_path, capture_output=True)
        assert not result.returncode and not result.stderr

    integrated_mtz = tmp_path / "integrated.mtz"
    assert integrated_mtz.is_file()

    mtz_object = iotbx.mtz.object(file_name=str(integrated_mtz))
    assert mtz_object.column_labels()[:14] == [
        "H",
        "K",
        "L",
        "M_ISYM",
        "BATCH",
        "IPR",
        "SIGIPR",
        "I",
        "SIGI",
        "BG",
        "SIGBG",
        "FRACTIONCALC",
        "XDET",
        "YDET",
    ]

    assert len(mtz_object.batches()) == 11
    batch = mtz_object.batches()[0]
    expected_unit_cell = uctbx.unit_cell((42.5787, 42.5787, 40.2983, 90, 90, 90))
    assert expected_unit_cell.is_similar_to(uctbx.unit_cell(list(batch.cell()))), (
        expected_unit_cell.parameters(),
        list(batch.cell()),
    )
    assert mtz_object.space_group().type().lookup_symbol() == "P 41 21 2"
    assert approx_equal(mtz_object.n_reflections(), 7446, eps=2e3)
