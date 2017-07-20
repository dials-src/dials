from __future__ import absolute_import, division
from dials.array_family import flex # import dependency


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "image_examples/SACLA_MPCCD_Cheetah")

  def run(self):
    self.test_sacla_h5()

  def test_sacla_h5(self):
    from os.path import join, exists
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    assert exists(join(self.path, 'run266702-0-subset.h5'))

    geometry_path = join(self.path, 'refined_experiments_level1.json')
    assert exists(geometry_path)

    f = open("process.phil", 'w')
    f.write("""
      input.reference_geometry=%s
      indexing {
        known_symmetry {
          space_group = P43212
          unit_cell = 78.9 78.9 38.1 90 90 90
        }
        method = fft1d
        refinement_protocol.d_min_start = 2.2
      }

      spotfinder {
        filter.min_spot_size = 2
        threshold {
          xds {
            gain = 5.46 # from dials.estimate_gain run266702-0-subset.h5 max_images=4
            global_threshold = 50
          }
        }
      }

      refinement {
        parameterisation {
          beam.fix = all
          detector.fix_list = Dist,Tau1
          auto_reduction {
            action = fix
            min_nref_per_parameter = 1
          }
          crystal {
            unit_cell {
              restraints {
                tie_to_target {
                  values = 78.9,78.9,38.1,90,90,90
                  sigmas = 1,1,1,0,0,0
                }
              }
            }
          }
        }
      }
      integration {
        integrator = stills
        profile.fitting = False
        background {
          algorithm = simple
          simple {
            model.algorithm = linear2d
            outlier.algorithm = tukey
          }
        }
      }
      profile {
        gaussian_rs {
          min_spots.overall = 0
        }
      }
      """%geometry_path)
    f.close()

    # Call dials.stills_process
    result = easy_run.fully_buffered([
      'dials.stills_process',
      join(self.path, 'run266702-0-subset.h5'),
      'process.phil',
    ]).raise_if_errors()
    result.show_stdout()

    import cPickle as pickle
    # Frame 1 no longer indexing after cctbx r25607 which made wavelengths be on a per-image basis
    #for result, n_refls in zip(["idx-run266702-0-subset_00000_integrated.pickle",
    #                            "idx-run266702-0-subset_00001_integrated.pickle"],
    #                            [range(109,114), range(80,85)]): # large ranges to handle platform-specific differences
    #for result, n_refls in zip(["idx-run266702-0-subset_00000_integrated.pickle"],
    #                            [range(109,114)]): # large ranges to handle platform-specific differences
    # dxtbx r25668 and 25669 flip X axis in the SACLA format class and changed indexing results.
    #for result, n_refls in zip(["idx-run266702-0-subset_00001_integrated.pickle"],
    #                            [range(90,96)]): # large ranges to handle platform-specific differences
    # 02/12/17 Handle change to stills_process refining after indexing plus new spotfinding params
    #for result, n_refls in zip(["idx-run266702-0-subset_00000_integrated.pickle",
    #                            "idx-run266702-0-subset_00001_integrated.pickle",
    #                            "idx-run266702-0-subset_00003_integrated.pickle"],
    #                            [range(75,90),range(220,230),range(285,295)]): # large ranges to handle platform-specific differences
    # 02/14/17 Further changes to stills_process: resetting rejected reflections before re-refinement
    #for result, n_refls in zip(["idx-run266702-0-subset_00000_integrated.pickle",
    #                            "idx-run266702-0-subset_00001_integrated.pickle",
    #                            "idx-run266702-0-subset_00003_integrated.pickle"],
    #                            [range(80,95),range(225,235),range(235,245)]): # large ranges to handle platform-specific differences
    # 02/21/17 Changes to stills_process: refine during indexing instead of after. Also used refined metrology from Rahel
    #for result, n_refls in zip(["idx-run266702-0-subset_00001_integrated.pickle",
    #                            "idx-run266702-0-subset_00003_integrated.pickle"],
    #                            [range(600,610),range(505,520)]): # large ranges to handle platform-specific differences
    # 04/25/17 Changes after reverting sign_error_27Feb2014_through_15Feb2017 in xfel/mono_simulation/max_like.py
    #for result, n_refls in zip(["idx-run266702-0-subset_00001_integrated.pickle",
    #                            "idx-run266702-0-subset_00003_integrated.pickle"],
    #                            [range(565,580),range(495,510)]): # large ranges to handle platform-specific differences
    # 09/20/17 Changes to still indexer: refine candidate basis vectors in target symmetry if supplied
    for result, n_refls in zip(["idx-run266702-0-subset_00001_integrated.pickle",
                                "idx-run266702-0-subset_00003_integrated.pickle"],
                                [range(100,115),range(155,165)]): # large ranges to handle platform-specific differences
      table = pickle.load(open(result, 'rb'))
      assert len(table) in n_refls, len(table)
      assert 'id' in table
      assert (table['id'] == 0).count(False) == 0
    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
