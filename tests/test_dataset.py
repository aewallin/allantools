from allantools.dataset import Dataset
from allantools import noise
import pytest


@pytest.fixture
def dataset():
    return Dataset(noise.white(10))


def test_no_function_in_allantools(dataset):
    with pytest.raises(AttributeError):
        dataset.compute("nosuchfunction")


def test_blacklisted_function(dataset):
    with pytest.raises(RuntimeError):
        dataset.compute("calc_mtotdev_phase")


def test_compute_functions(dataset):
    dataset.compute("adev")
    dataset.compute("oadev")
    dataset.compute("mdev")
    dataset.compute("hdev")
    dataset.compute("ohdev")
    dataset.compute("tdev")
    dataset.compute("totdev")
    dataset.compute("mtotdev")
    dataset.compute("ttotdev")
    dataset.compute("htotdev")
    dataset.compute("theo1")
    dataset.compute("mtie")
    dataset.compute("tierms")
