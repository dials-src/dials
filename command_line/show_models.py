from __future__ import division

def run(args):


  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks
  from libtbx.utils import Abort

  parser = OptionParser(
    read_experiments=True,
    read_datablocks=True,
    read_datablocks_from_images=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  datablocks = flatten_datablocks(params.input.datablock)

  if experiments is not None:
    for detector in experiments.detectors():
      print detector
    for beam in experiments.beams():
      print beam
    for scan in experiments.scans():
      print scan
    for goniometer in experiments.goniometers():
      print goniometer
    for crystal in experiments.crystals():
      crystal.show(show_scan_varying=True)
  if datablocks is not None:
    for datablock in datablocks:
      if datablock.format_class() is not None:
        print 'Format: %s' %datablock.format_class()
      imagesets = datablock.extract_imagesets()
      for imageset in imagesets:
        try: print imageset.get_template()
        except Exception: pass
        print imageset.get_detector()
        print imageset.get_beam()
        if imageset.get_scan() is not None:
          print imageset.get_scan()
        if imageset.get_goniometer() is not None:
          print imageset.get_goniometer()
  if experiments is None and datablocks is None:
    raise Abort('No experiments or datablocks specified')
  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
