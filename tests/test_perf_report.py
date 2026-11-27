"""Unit tests for the SLURM perf-report parsers (no sacct needed)."""
from scripts.perf_report import hms_to_s, size_to_bytes


def test_hms_to_s():
    assert hms_to_s("00:14:13") == 853
    assert hms_to_s("03:43:54") == 3 * 3600 + 43 * 60 + 54
    assert hms_to_s("1-02:03:04") == 86400 + 2 * 3600 + 3 * 60 + 4
    assert hms_to_s("") == 0
    assert hms_to_s("Unknown") == 0


def test_size_to_bytes():
    assert size_to_bytes("1024K") == 1024 * 1024
    assert size_to_bytes("512M") == 512 * 1024 ** 2
    assert size_to_bytes("1.5G") == int(1.5 * 1024 ** 3)
    assert size_to_bytes("0") == 0
    assert size_to_bytes("") == 0
