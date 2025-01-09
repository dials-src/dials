from __future__ import annotations

import procrunner


def test_indexed_as_integrated(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed", pathlib=True)
    refl = data_dir / "indexed.refl"

    command = [
        "dev.dials.pixel_histogram",
        refl,
    ]
    result = procrunner.run(command, working_directory=tmp_path)
    assert not result.returncode and not result.stderr
    assert (tmp_path / "histogram.json").is_file()
