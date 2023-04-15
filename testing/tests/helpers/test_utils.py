import tempfile
from PIL import Image
import numpy
from pathlib import Path
import pytest

from pymol import cmd


def datafile(filename: str) -> Path:
    """
    Return path to filename, the current directory and the data
    directory are searched for filename.
    :filename: filename to search for
    """
    if Path(filename).exists():
        return filename
    PYMOL_TESTING_ROOT = 2
    # pymol-testing root
    pymol_test_dir = Path(__file__).parents[PYMOL_TESTING_ROOT]
    return pymol_test_dir.joinpath('data', filename)


class mktemp(object):
    """
    Context manager which returns a temporary filename and
    deletes the file in the end, if it exists.
    """

    def __init__(self, suffix: str = ''):
        self.file = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
        self.file.close()  # Allow other processes to access the file
        self.filename = self.file.name

    def __enter__(self):
        return self.filename

    def __exit__(self, exc_type, exc_value, traceback):
        if Path(self.filename).exists():
            Path(self.filename).unlink()


def ambientOnly(cmd) -> None:
    """
    Set up a scene with only ambient lighting.
    """
    cmd.set('ambient', 1)
    cmd.set('antialias', 0)
    cmd.set('light_count', 1)
    cmd.set('depth_cue', 0)
    # needed for open-source
    cmd.set('reflect', 0)
    cmd.set('direct', 0)


def get_imagearray(cmd, **kwargs) -> numpy.array:
    """
    Return the image as a numpy array.
    :kwargs: keyword arguments passed to cmd.png
    """

    with mktemp('.png') as filename:
        cmd.png(filename, **kwargs)
        img = Image.open(filename)

        return numpy.array(img.getdata(),
                           numpy.uint8).reshape((img.size[1], img.size[0], -1))


def _imageHasColor(cmd, color, img, delta=0) -> bool:
    """
    Return True if the image contains the given color.
    :color: color to search for
    :img: image as numpy array
    :delta: maximum difference between color and image color
    """
    if isinstance(color, str):
        color = [int(v*255) for v in cmd.get_color_tuple(color)]
    else:
        color = list(color)
    dim = img.shape[-1]
    if dim == len(color) + 1:
        dim -= 1
        img = img[..., :dim]
    diff = abs(img.reshape((-1, dim)) - color)
    return (diff - delta <= 0).prod(1).sum()


def imageHasColor(cmd, color: str, delta: float = 0) -> bool:
    """
    Return True if the image contains the given color.
    :color: color to search for
    :delta: maximum difference between color and image color
    """
    img = get_imagearray(cmd)
    return _imageHasColor(cmd, color, img, delta)


def assert_in_names_undo(cmd, name: str) -> None:
    """
    Assert that the object is in the names list and undoing the
    command removes it from the list.
    :name: name to search for
    """
    assert name in cmd.get_names()
    cmd.undo()
    assert name not in cmd.get_names()
    cmd.redo()
    assert name in cmd.get_names()


def compatible_with(version: str) -> bool:
    def tupleize_version(str_: str):
        return tuple(int(x) for x in str_.partition('.')[0::2] if x.isdigit())

    PYMOL_VERSION = cmd.get_version()
    PYMOL_VERSION_TUPLE = tupleize_version(PYMOL_VERSION[0])

    if isinstance(version, int):
        return version <= PYMOL_VERSION[2]
    elif isinstance(version, float):
        return version <= PYMOL_VERSION[1]
    else:
        return tupleize_version(version) <= PYMOL_VERSION_TUPLE


def requires_version(version, reason=None):
    def decorator(test_func):
        return pytest.mark.skipif(
            not compatible_with(version),
            reason=reason or f"Requires PyMOL {version}"
        )(test_func)
    return decorator
