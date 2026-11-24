#!/bin/bash
# Scaling + reproducibility matrix on the simulated dataset.
#
# For each node count in NODES, run the pipeline twice (pair a/b) sequentially,
# collect the io/mem monitoring traces, diff the pair's science products, and
# append job timings to scaling.csv.
#
# Env overrides (defaults shown):
#   SIM_PRISTINE=$SCRATCH/ddf-sim-dataset/Sim     pristine simulated dataset
#   SCALING_ROOT=$SCRATCH/ddf-sim-dataset/scaling  per-run dirs + scaling.csv
#   NODES="3 6 12"                                 add 1 to redo the 1-node pair
#   TIME="00:40:00"
#   CFG=$REPO_DIR/examples/Simul-Hackaton.cfg
#   LAUNCHER=$REPO_DIR/slurm/jean-zay.slurm      (simulated launch; --nodes overrides it)
#   MONITOR=$HOME/ddf-monitor                      where mem*/iomonitor* are moved
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/.." && pwd)"
SIM_PRISTINE="${SIM_PRISTINE:-${SCRATCH:-$HOME/scratch}/ddf-sim-dataset/Sim}"
SCALING_ROOT="${SCALING_ROOT:-${SCRATCH:-$HOME/scratch}/ddf-sim-dataset/scaling}"
NODES="${NODES:-3 6 12}"
TIME="${TIME:-00:40:00}"
CFG="${CFG:-$REPO_DIR/examples/Simul-Hackaton.cfg}"
LAUNCHER="${LAUNCHER:-$REPO_DIR/slurm/jean-zay.slurm}"
MONITOR="${MONITOR:-$HOME/ddf-monitor}"
EXCLUDE="${EXCLUDE:-mem,iomonitor*,logs/*,*.log,*.ddfcache/*,table.*,*.hash}"
CSV="$SCALING_ROOT/scaling.csv"

[ -d "$SIM_PRISTINE" ] || { echo "error: pristine sim data not found: $SIM_PRISTINE" >&2; exit 1; }
[ -f "$CFG" ] || { echo "error: cfg not found: $CFG" >&2; exit 1; }
mkdir -p "$SCALING_ROOT" "$MONITOR"
[ -f "$CSV" ] || echo "nodes,pair,jobid,elapsed,exit" > "$CSV"

for n in $NODES; do
    for pair in a b; do
        run="$SCALING_ROOT/${n}${pair}"
        mon="$MONITOR/${n}${pair}"
        echo "== ${n} node(s), pair ${pair} -> ${run}"
        rm -rf "$run"
        mkdir -p "$run" "$mon"
        cp -r "$SIM_PRISTINE"/. "$run"/

        set +e
        submit=$(DDF_SIM_DIR="$run" DDF_PIPELINE_CFG="$CFG" DDF_SMOKE_CHECK=1 \
                 sbatch --export=ALL --wait --nodes="$n" --time="$TIME" "$LAUNCHER" 2>&1)
        rc=$?
        set -e
        jobid=$(printf '%s\n' "$submit" | grep -oE '[0-9]+' | tail -1)
        elapsed=$(sacct -j "$jobid" --format=Elapsed -Pn 2>/dev/null | head -1 | tr -d ' ')

        mv "$HOME"/mem* "$HOME"/iomonitor* "$mon/" 2>/dev/null || true
        echo "   job ${jobid:-?} elapsed=${elapsed:-?} rc=$rc"
        echo "$n,$pair,${jobid:-?},${elapsed:-?},$rc" >> "$CSV"
    done

    echo "   reproducibility ${n}a vs ${n}b (science only):"
    python3 "$REPO_DIR/scripts/compare_runs.py" manifest "$SCALING_ROOT/${n}a" "$SCALING_ROOT/${n}a.manifest" --exclude "$EXCLUDE" >/dev/null
    python3 "$REPO_DIR/scripts/compare_runs.py" manifest "$SCALING_ROOT/${n}b" "$SCALING_ROOT/${n}b.manifest" --exclude "$EXCLUDE" >/dev/null
    python3 "$REPO_DIR/scripts/compare_runs.py" compare "$SCALING_ROOT/${n}a.manifest" "$SCALING_ROOT/${n}b.manifest" | sed 's/^/     /'
done

echo
echo "== $CSV =="
cat "$CSV"
