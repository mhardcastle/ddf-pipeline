# CI runs on GitHub Actions, not GitLab

ddf-pipeline's canonical repository is on GitHub; the GitLab repository is a read-only
mirror. We put CI in GitHub Actions (`.github/workflows/ci.yml`) so there is a single
source of truth and no duplicated CI config.

Considered options: GitLab CI (rejected — it would duplicate config on a mirror with no
other role). The heavy smoke job (needs ~40 CPUs) will later run on a GitHub
larger/self-hosted runner, still within Actions.
