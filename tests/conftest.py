import os
import sys

import pytest

from sourmash_tst_utils import TempDirectory, RunnerContext
sys.stdout = sys.stderr

@pytest.fixture
def runtmp():
    with TempDirectory() as location:
        yield RunnerContext(location)

# Set environment variable PYTEST_RUNNING
def pytest_configure(config):
    os.environ["PYTEST_RUNNING"] = "1"

def pytest_unconfigure(config):
    del os.environ["PYTEST_RUNNING"]
