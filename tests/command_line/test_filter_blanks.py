from __future__ import annotations

from dxtbx.model.experiment_list import ExperimentList

from dials.array_family import flex
from dials.command_line import filter_blanks


def test_strong(dials_data, capsys, run_in_tmp_path):
    expts_file = dials_data("insulin_processed", pathlib=True) / "imported.expt"
    refl_file = dials_data("insulin_processed", pathlib=True) / "strong.refl"
    refl = flex.reflection_table.from_file(refl_file)
    z = refl["xyzobs.px.value"].parts()[2]
    refl_subset = refl.select(((z > 10) & (z < 20)) | ((z > 30) & (z < 40)))
    refl_subset.as_file(run_in_tmp_path / "mod.refl")

    filter_blanks.run([str(expts_file), str(run_in_tmp_path / "mod.refl")])

    # FIXME add reading if output experiment file

    el = ExperimentList.from_file(str(run_in_tmp_path / "not_blank.expt"))
    assert len(el) == 2
    assert el[0].scan.get_array_range() == (11, 20)
    assert el[1].scan.get_array_range() == (29, 40)
