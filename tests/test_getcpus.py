from utils.getcpus import getcpus, get_physical_cpus


def test_getcpus_counts_pbs_nodefile(tmp_path, monkeypatch):
    nodefile = tmp_path / "nodes"
    nodefile.write_text("node1\nnode2\nnode3\n")
    monkeypatch.setenv("PBS_NODEFILE", str(nodefile))
    assert getcpus() == 3


def test_getcpus_falls_back_to_affinity(monkeypatch):
    monkeypatch.delenv("PBS_NODEFILE", raising=False)
    n = getcpus()
    assert n == get_physical_cpus()
    assert n >= 1
