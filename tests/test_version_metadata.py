from importlib.metadata import version

import compactlib


def test_package_version_matches_metadata():
    assert compactlib.__version__ == version("compact-dia-library")
