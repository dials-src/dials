#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import os
import pytest

@pytest.fixture(scope="session")
def dials_regression():
  '''Return the absolute path to the dials_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    import dials_regression as dr
  except ImportError:
    pytest.skip("dials_regression required for this test")
  return os.path.dirname(dr.__file__)

@pytest.fixture(scope="session")
def xia2_regression():
  '''Return the absolute path to the xia2_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    import xia2_regression as xr
  except ImportError:
    pytest.skip("xia2_regression required for this test")
  return os.path.dirname(xr.__file__)

@pytest.fixture(scope="session")
def xia2_regression_build():
  '''Return the absolute path to the xia2_regression directory within the build
     path as a string. Skip the test if xia2_regression is not installed.'''
  try:
    x2rpath = os.path.join(os.environ.get('LIBTBX_BUILD'), 'xia2_regression')
  except AttributeError:
    x2rpath = ''
  if not os.path.exists(x2rpath):
    pytest.skip("xia2_regression required for this test")
  if 'test_data' not in os.listdir(x2rpath):
    pytest.skip("xia2_regression files need to be downloaded for this test. Run xia2_regression.fetch_test_data")
  return x2rpath

from libtbx.test_utils.pytest import libtbx_collector
pytest_collect_file = libtbx_collector()
