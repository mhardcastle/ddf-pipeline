# Datasets and data access

The DDF Pipeline experiments use three datasets plus a shared set of catalogs. All
downloads are public shares on the CNRS sdrive. Keep code and compressed data on `WORK`;
decompress into `SCRATCH` for each experiment.

> **Restart caveat.** The pipeline's restart mechanism is unreliable. If a run is
> interrupted, clear the execution directory and decompress a fresh dataset before
> re-running.

## Catalogue

### LoTSS L234026 "big" — realistic, several GB
- Measurement Sets: <https://sdrive.cnrs.fr/s/azWk38oXb6BwrRS>
- List `mslist.txt`: <https://sdrive.cnrs.fr/s/YjPXwZJ62FwX6Jo>
- List `big-mslist.txt`: <https://sdrive.cnrs.fr/s/NRHNdrLGZs6o7a4>
- Configuration `lotss_mpi.cfg`: <https://sdrive.cnrs.fr/s/omnAgeZL5cPHAAn>

### Simulated LOFAR "small" — the smoke test
- Measurement Sets and lists `Sim.tgz`: <https://sdrive.cnrs.fr/s/PNmRiW3zc7XmEwN>
- Configuration `Simul-Hackaton.cfg`: <https://sdrive.cnrs.fr/s/qgptXzGa6BCE56J>

### Simulated SKA Low AAstar
- Measurement Sets: on Adastra, `$WORK_cad14221_SHARED` (no download)
- Configuration `skalow_aastar_mini.cfg`: <https://sdrive.cnrs.fr/s/Rn6sdE5n85K8SjY>

### Catalogs
- `CATALOGS.zip`: <https://sdrive.cnrs.fr/s/awWCRsaoeax8t4K>

## Downloading

Download from the link, or with `wget` against the DAV endpoint
(`/public.php/dav/files/<id>`):

```bash
wget -O mslist.txt             https://sdrive.cnrs.fr/public.php/dav/files/YjPXwZJ62FwX6Jo
wget -O big-mslist.txt         https://sdrive.cnrs.fr/public.php/dav/files/NRHNdrLGZs6o7a4
wget -O lotss_mpi.cfg          https://sdrive.cnrs.fr/public.php/dav/files/omnAgeZL5cPHAAn
wget -O Sim.tgz                https://sdrive.cnrs.fr/public.php/dav/files/PNmRiW3zc7XmEwN
wget -O Simul-Hackaton.cfg     https://sdrive.cnrs.fr/public.php/dav/files/qgptXzGa6BCE56J
wget -O skalow_aastar_mini.cfg https://sdrive.cnrs.fr/public.php/dav/files/Rn6sdE5n85K8SjY
wget -O CATALOGS.zip           https://sdrive.cnrs.fr/public.php/dav/files/awWCRsaoeax8t4K
```

The **big** Measurement Sets are a folder share (several GB); download them from the web
link on `WORK`. The **SKA Low** Measurement Sets are already on Adastra.

## First run (simulated "small")

1. **Install** the software on Jean Zay with `install_uv_jean_zay.sh` (loads the modules,
   clones DDFacet + KillMS, `uv sync --extra mpi-support`).
2. **Install the data** with `install_data_jean_zay.sh` — it downloads and SHA256-verifies
   `CATALOGS.zip`, `Sim.tgz` and `Simul-Hackaton.cfg`.
3. **Prepare** — code and compressed data on `WORK`; decompress to `SCRATCH` per experiment.
   If interrupted, clear the directory and re-decompress.
4. **Run** — `pipeline.py <config.cfg>`, or on a cluster `sbatch slurm/<platform>.slurm`
   from the repo root.
5. **Verify** — the run succeeded if `image_full_ampphase_di_m.NS.app.restored.fits` exists.

> The big dataset is several GB and realistic; expect it to take much longer than the
> simulated one.

## Scripts

- `install_uv_<platform>.sh` — platform environment (modules + `uv sync`); the MPI/compiler
  modules differ per cluster.
- `install_data_jean_zay.sh` — the simulated "small" dataset + catalogs, SHA256-pinned.
- `install_data.sh <big|skalow|catalogs>` — the other datasets and the shared catalogs,
  SHA256-pinned. Fill the `*_SHA256` pins from a login node before first use.
