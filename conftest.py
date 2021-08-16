#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#


import multiprocessing
import os
import sys
import warnings

import pytest

# https://stackoverflow.com/a/40846742
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

if sys.version_info[:2] == (3, 7) and sys.platform == "darwin":
    multiprocessing.set_start_method("forkserver")

collect_ignore = []


def pytest_configure(config):
    if not config.pluginmanager.hasplugin("dials_data"):

        @pytest.fixture(scope="session")
        def dials_data():
            pytest.skip("This test requires the dials_data package to be installed")

        globals()["dials_data"] = dials_data


@pytest.fixture(scope="session")
def dials_regression():
    """Return the absolute path to the dials_regression module as a string.
    Skip the test if dials_regression is not installed."""

    if "DIALS_REGRESSION" in os.environ:
        return os.environ["DIALS_REGRESSION"]

    try:
        import dials_regression as dr

        return os.path.dirname(dr.__file__)
    except ImportError:
        pass  # dials_regression not configured
    try:
        import socket

        reference_copy = "/dls/science/groups/scisoft/DIALS/repositories/git-reference/dials_regression"
        if (
            os.name == "posix"
            and socket.gethostname().endswith(".diamond.ac.uk")
            and os.path.exists(reference_copy)
        ):
            return reference_copy
    except ImportError:
        pass  # Cannot tell whether in DLS network or not
    pytest.skip("dials_regression required for this test")


@pytest.fixture
def run_in_tmpdir(tmpdir):
    """Shortcut to create a temporary directory and then run the test inside
    this directory."""
    cwd = os.getcwd()
    tmpdir.chdir()
    yield tmpdir
    os.chdir(cwd)
