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
EXCLUDE="${EXCLUDE:-DDF_monitor*,DDF_register*,DDF_io_monitor*,mem,iomonitor*,logs/*,*.log,*.ddfcache/*,table.*,*.hash,local_*mslist.txt,ClipCal.last,last_killMS.obj,pop.myPickle,BatchCurrentParset.parset,solslist_*.txt}"
CSV="$SCALING_ROOT/scaling.csv"
PY="$REPO_DIR/.venv/bin/python"; [ -x "$PY" ] || PY=python3

[ -d "$SIM_PRISTINE" ] || { echo "error: pristine sim data not found: $SIM_PRISTINE" >&2; exit 1; }
[ -f "$CFG" ] || { echo "error: cfg not found: $CFG" >&2; exit 1; }
mkdir -p "$SCALING_ROOT" "$MONITOR"
# Snapshot the config: a re-install or edit during the matrix must not change a run pair.
CFG_SNAPSHOT="$SCALING_ROOT/$(basename "$CFG")"
cp "$CFG" "$CFG_SNAPSHOT"
CFG="$CFG_SNAPSHOT"
[ -f "$CSV" ] || echo "nodes,pair,jobid,elapsed,exit" > "$CSV"

for n in $NODES; do
    for pair in a b; do
        run="$SCALING_ROOT/${n}${pair}"
        mon="$MONITOR/${n}${pair}"
        echo "== ${n} node(s), pair ${pair} -> ${run}"
        rm -rf "$run"
        mkdir -p "$run" "$mon"
        cp -r "$SIM_PRISTINE"/. "$run"/
        find "$run" -name table.lock -delete   # stale locks from the archive deadlock casacore on Lustre
        # Every MS named in the run's lists must exist (a partial copy or an
        # extraction still in progress surfaces later as "MS does not exist").
        for msf in "$run/mslist.txt" "$run/big-mslist.txt"; do
            [ -f "$msf" ] || continue
            while read -r ms; do
                [ -z "$ms" ] && continue
                if [ ! -d "$run/$ms" ]; then
                    echo "error: $run/$ms missing (partial copy or extraction in progress)" >&2
                    exit 1
                fi
            done < "$msf"
        done

        set +e
        submit=$(DDF_SIM_DIR="$run" DDF_PIPELINE_CFG="$CFG" DDF_SMOKE_CHECK=1 \
                 sbatch --export=ALL --wait --nodes="$n" --time="$TIME" "$LAUNCHER" 2>&1)
        rc=$?
        set -e
        jobid=$(printf '%s\n' "$submit" | grep -oE '[0-9]+' | tail -1)
        elapsed=$(sacct -j "$jobid" --format=Elapsed -Pn 2>/dev/null | head -1 | tr -d ' ')

        mv "$HOME"/mem* "$HOME"/iomonitor* "$HOME"/DDF_* "$mon/" 2>/dev/null || true
        echo "   job ${jobid:-?} elapsed=${elapsed:-?} rc=$rc"
        echo "$n,$pair,${jobid:-?},${elapsed:-?},$rc" >> "$CSV"
    done

    echo "   reproducibility ${n}a vs ${n}b (science only):"
    python3 "$REPO_DIR/scripts/compare_runs.py" manifest "$SCALING_ROOT/${n}a" "$SCALING_ROOT/${n}a.manifest" --exclude "$EXCLUDE" >/dev/null
    python3 "$REPO_DIR/scripts/compare_runs.py" manifest "$SCALING_ROOT/${n}b" "$SCALING_ROOT/${n}b.manifest" --exclude "$EXCLUDE" >/dev/null
    python3 "$REPO_DIR/scripts/compare_runs.py" compare "$SCALING_ROOT/${n}a.manifest" "$SCALING_ROOT/${n}b.manifest" | sed 's/^/     /'
    echo "   statistical ${n}a vs ${n}b (key products):"
    for f in image_full_ampphase_di_m.NS.app.restored.fits image_full_ampphase_di_m.NS.int.residual.fits; do
        if [ -f "$SCALING_ROOT/${n}a/$f" ] && [ -f "$SCALING_ROOT/${n}b/$f" ]; then
            "$PY" "$REPO_DIR/scripts/compare_products.py" "$SCALING_ROOT/${n}a/$f" "$SCALING_ROOT/${n}b/$f" | sed 's/^/     /'
        fi
    done
    sa="$SCALING_ROOT/${n}a/SOLSDIR/0000.p00.MS/killMS.DDS3_full.sols.npz"
    sb="$SCALING_ROOT/${n}b/SOLSDIR/0000.p00.MS/killMS.DDS3_full.sols.npz"
    if [ -f "$sa" ] && [ -f "$sb" ]; then
        echo "   solutions ${n}a vs ${n}b (phase-invariant):"
        "$PY" "$REPO_DIR/scripts/compare_sols.py" "$sa" "$sb" | sed 's/^/     /'
    fi
done

echo
echo "== $CSV =="
cat "$CSV"

echo
echo "== performance (sacct, pipeline step) =="
"$PY" "$REPO_DIR/scripts/perf_report.py" --from-csv "$CSV"
