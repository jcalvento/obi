import os

from tests.utils import results_dir


def pytest_sessionstart(session):
    to_create = results_dir()
    if not os.path.isdir(to_create):
        os.mkdir(to_create)
