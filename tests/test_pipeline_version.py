from utils.pipeline_version import version


def test_version_returns_non_empty_string():
    result = version()
    assert isinstance(result, str)
    assert result != ""
