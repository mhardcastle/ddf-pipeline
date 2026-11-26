"""Test that scripts/repackage_dataset.sh strips stale table locks and keeps data.

The archived datasets carry `table.lock` files (one per table/sub-table) that make
casacore deadlock on Lustre. The repackager must remove them and leave the
StandardStMan/TiledStMan data files (`table.f0`, ...) untouched.
"""
import os
import subprocess
import tarfile

SCRIPT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts",
                                      "repackage_dataset.sh"))


def _make_source(path):
    """Build a tiny Sim.tgz: Sim/x.MS/{table.lock, table.f0, table.dat}."""
    src = path / "src" / "Sim" / "x.MS"
    src.mkdir(parents=True)
    (src / "table.lock").write_bytes(b"")
    (src / "table.f0").write_bytes(b"data-f0")
    (src / "table.dat").write_bytes(b"data-dat")
    tgz = path / "Sim.tgz"
    with tarfile.open(tgz, "w:gz") as tf:
        tf.add(path / "src" / "Sim", arcname="Sim")
    return tgz


def test_repackage_removes_locks_and_keeps_data(tmp_path):
    src = _make_source(tmp_path)
    out = tmp_path / "Sim_clean.tgz"

    result = subprocess.run(["bash", SCRIPT, str(src), str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "removed=1" in result.stdout
    assert "remaining=0" in result.stdout

    with tarfile.open(out, "r:gz") as tf:
        names = set(tf.getnames())

    assert "Sim/x.MS/table.lock" not in names
    assert "Sim/x.MS/table.f0" in names
    assert "Sim/x.MS/table.dat" in names


def test_repackage_prints_sha256(tmp_path):
    src = _make_source(tmp_path)
    out = tmp_path / "Sim_clean.tgz"
    result = subprocess.run(["bash", SCRIPT, str(src), str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "sha256 " in result.stdout
