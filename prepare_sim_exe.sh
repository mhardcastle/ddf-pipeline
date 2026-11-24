#!/bin/bash
# Fetch the smoke-test inputs (catalogs, simulated dataset, configuration file)
# and verify each one against a pinned SHA256 before it is used, so a run is
# reproducible from a clone: an unpinned or tampered download is never extracted.
set -euo pipefail

SCRIPT_DIR="$(dirname "$(realpath -- "${BASH_SOURCE[0]}")")"
cd "$SCRIPT_DIR"

# --- Pinned inputs ------------------------------------------------------------
# Keep the URLs; the SHA256 pins make what they serve immutable.

# Catalogs, extracted to ./CATALOGS ($DDF_PIPELINE_CATALOGS)
CATALOGS_URL="https://sdrive.cnrs.fr/public.php/dav/files/awWCRsaoeax8t4K/CATALOGS.zip"
CATALOGS_SHA256="TODO"

# Simulated dataset, extracted to $SCRATCH/ddf-sim-dataset
SIM_URL="https://sdrive.cnrs.fr/public.php/dav/files/PNmRiW3zc7XmEwN"
SIM_SHA256="TODO"

# Pipeline configuration working with the simulated dataset, installed to examples/
CFG_URL="https://sdrive.cnrs.fr/public.php/dav/files/qgptXzGa6BCE56J"
CFG_SHA256="TODO"

# To fill a pin, download once on a login node and hash the three files:
#   wget -O CATALOGS.zip <CATALOGS_URL>   # and the two others, then
#   sha256sum CATALOGS.zip Sim.tgz Simul-Hackaton.cfg

# Dataset destination: $SCRATCH on the cluster, the repository root elsewhere.
SIM_DIR="${SCRATCH:-$SCRIPT_DIR}/ddf-sim-dataset"

die() {
    echo "error: $*" >&2
    exit 1
}

# verify <file> <sha256> — 0 on match, 1 on mismatch (with both values shown).
verify() {
    local file=$1 pin=$2 actual
    actual=$(sha256sum "$file" | cut -d' ' -f1)
    if [ "$actual" = "$pin" ]; then
        echo "ok   $file"
        return 0
    fi
    echo "FAIL $file" >&2
    echo "     pinned $pin" >&2
    echo "     actual $actual" >&2
    return 1
}

# require_pinned <file> <sha256> — refuse placeholders and malformed pins.
require_pinned() {
    local file=$1 pin=$2
    if [ "$pin" = "TODO" ]; then
        die "$file has no pinned SHA256. Hash it on a login node (sha256sum $file)
and set the matching *_SHA256 at the top of $0."
    fi
    if [ ${#pin} -ne 64 ] || [[ $pin = *[!0-9a-f]* ]]; then
        die "$file: pinned SHA256 is not 64 lowercase hex characters: '$pin'"
    fi
}

# fetch_pinned <url> <dest> <sha256> — reuse <dest> only if it matches the pin,
# otherwise download to <dest>.part, verify, then replace <dest>.
fetch_pinned() {
    local url=$1 dest=$2 pin=$3
    if [ -f "$dest" ] && verify "$dest" "$pin" >/dev/null 2>&1; then
        echo "reuse $dest"
        return 0
    fi

    echo "fetch $url"
    rm -f "$dest.part"
    wget -nv -O "$dest.part" "$url"
    if ! verify "$dest.part" "$pin"; then
        rm -f "$dest.part"
        die "$url: checksum mismatch, $dest was left untouched"
    fi
    mv "$dest.part" "$dest"
}

# Refuse placeholders and malformed pins before touching the network.
require_pinned CATALOGS.zip "$CATALOGS_SHA256"
require_pinned Sim.tgz "$SIM_SHA256"
require_pinned Simul-Hackaton.cfg "$CFG_SHA256"

# Install the catalogs
fetch_pinned "$CATALOGS_URL" CATALOGS.zip "$CATALOGS_SHA256"
unzip -o -q CATALOGS.zip

# Install the simulated dataset
fetch_pinned "$SIM_URL" Sim.tgz "$SIM_SHA256"
mkdir -p "$SIM_DIR"
tar -xvf Sim.tgz --directory "$SIM_DIR"

# Install the ddf pipeline configuration file working with the simulated dataset
fetch_pinned "$CFG_URL" Simul-Hackaton.cfg "$CFG_SHA256"
mv Simul-Hackaton.cfg examples/

echo "smoke inputs ready: CATALOGS/, examples/Simul-Hackaton.cfg, $SIM_DIR"
