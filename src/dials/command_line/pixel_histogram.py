# LIBTBX_SET_DISPATCHER_NAME dev.dials.pixel_histogram
from __future__ import annotations

import json
import sys

from dials.array_family import flex


def pixel_histogram(spot_filename: str, maximum: int = 0x10000) -> None:
    """Compute a histogram of pixel values up to the given maximum value"""

    spots = flex.reflection_table.from_file(spot_filename)
    spots = spots.select(spots.get_flags(spots.flags.indexed))

    if "shoebox" not in spots:
        sys.exit(f"no shoeboxes in {spot_filename}")

    histogram = flex.histogram(
        flex.double(), data_min=0, data_max=maximum, n_slots=maximum
    )

    for j in range(len(spots)):
        box = spots[j]["shoebox"].data
        histogram.update(
            flex.histogram(
                box.as_double().as_1d(), data_min=0, data_max=maximum, n_slots=maximum
            )
        )

    pixels = flex.ceil(histogram.slot_centers()).iround()
    values = histogram.slots()
    sel = values > 0
    with open("histogram.json", "w") as f:
        json.dump(
            {"pixels": list(pixels.select(sel)), "counts": list(values.select(sel))}, f
        )


if __name__ == "__main__":
    pixel_histogram(sys.argv[1])
