from utils.pipeline_version import version


def test_version_is_a_string():
    assert isinstance(version(), str)
