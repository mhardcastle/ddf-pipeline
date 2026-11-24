# What runs where: fast CI on GitHub, smoke tests on the HPC

Fast, always-green checks (unit tests, `compileall`, `mkdocs build --strict`) run on
GitHub Actions on every push and pull request. The smoke tests — running `pipeline.py`
on the simulated dataset under MPI — run on Jean Zay and Adastra through the `hpc`
agent, because they need the sibling DDFacet/KillMS checkouts, system casacore, and a
cluster node that GitHub-hosted runners cannot provide.

Considered options: running the smoke tests on GitHub larger/self-hosted runners
(deferred — the HPC path already exists and is where the data lives).
