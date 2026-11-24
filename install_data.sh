#!/bin/bash
# Fetch datasets and catalogs for DDF Pipeline experiments.
# (The simulated "small" set is installed by install_data_jean_zay.sh.)
# Every download is verified against a pinned SHA256 before it is used.
#
# Usage: ./install_data.sh <big|skalow|catalogs>
set -euo pipefail

SCRIPT_DIR="$(dirname "$(realpath -- "${BASH_SOURCE[0]}")")"
cd "$SCRIPT_DIR"

die() { echo "error: $*" >&2; exit 1; }

verify() {
    local file=$1 pin=$2 actual
    actual=$(sha256sum "$file" | cut -d' ' -f1)
    [ "$actual" = "$pin" ] && { echo "ok   $file"; return 0; }
    echo "FAIL $file" >&2
    echo "     pinned $pin" >&2
    echo "     actual $actual" >&2
    return 1
}

require_pinned() {
    local file=$1 pin=$2
    [ "$pin" = "TODO" ] && die "$file has no pinned SHA256. Hash it on a login node and set the pin at the top of $0."
    { [ ${#pin} -eq 64 ] && [ -z "${pin//[0-9a-f]/}" ]; } || die "$file: pin is not 64 lowercase hex: '$pin'"
}

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

# --- shared catalogs ---------------------------------------------------------
CATALOGS_URL="https://sdrive.cnrs.fr/public.php/dav/files/awWCRsaoeax8t4K"
CATALOGS_SHA256="8a0edbebc67613266251f9871b09c588c1c8102368daf1161045ecc26d4d701d"

# --- LoTSS L234026 "big" (lists + config; the MSs are a folder share) --------
MSLIST_URL="https://sdrive.cnrs.fr/public.php/dav/files/YjPXwZJ62FwX6Jo"
MSLIST_SHA256="TODO"
BIGMSLIST_URL="https://sdrive.cnrs.fr/public.php/dav/files/NRHNdrLGZs6o7a4"
BIGMSLIST_SHA256="TODO"
LOTSS_CFG_URL="https://sdrive.cnrs.fr/public.php/dav/files/omnAgeZL5cPHAAn"
LOTSS_CFG_SHA256="TODO"

# --- Simulated SKA Low AAstar (config only; the MSs are on Adastra) ----------
SKALOW_CFG_URL="https://sdrive.cnrs.fr/public.php/dav/files/Rn6sdE5n85K8SjY"
SKALOW_CFG_SHA256="TODO"

install_catalogs() {
    require_pinned CATALOGS.zip "$CATALOGS_SHA256"
    fetch_pinned "$CATALOGS_URL" CATALOGS.zip "$CATALOGS_SHA256"
    unzip -o -q CATALOGS.zip
}

case "${1:-}" in
  catalogs)
    install_catalogs
    ;;
  big)
    require_pinned mslist.txt "$MSLIST_SHA256"
    require_pinned big-mslist.txt "$BIGMSLIST_SHA256"
    require_pinned lotss_mpi.cfg "$LOTSS_CFG_SHA256"
    fetch_pinned "$MSLIST_URL" mslist.txt "$MSLIST_SHA256"
    fetch_pinned "$BIGMSLIST_URL" big-mslist.txt "$BIGMSLIST_SHA256"
    fetch_pinned "$LOTSS_CFG_URL" examples/lotss_mpi.cfg "$LOTSS_CFG_SHA256"
    install_catalogs
    cat <<'NOTE'
LoTSS "big": mslist.txt, big-mslist.txt and examples/lotss_mpi.cfg installed.
The Measurement Sets are a multi-GB folder share: download them from
  https://sdrive.cnrs.fr/s/azWk38oXb6BwrRS
on WORK, then point the mslist at them.
NOTE
    ;;
  skalow)
    require_pinned skalow_aastar_mini.cfg "$SKALOW_CFG_SHA256"
    fetch_pinned "$SKALOW_CFG_URL" examples/skalow_aastar_mini.cfg "$SKALOW_CFG_SHA256"
    install_catalogs
    cat <<'NOTE'
SKA Low AAstar: examples/skalow_aastar_mini.cfg installed.
The Measurement Sets are on Adastra at $WORK_cad14221_SHARED — no download needed.
NOTE
    ;;
  *)
    die "usage: $0 <big|skalow|catalogs>"
    ;;
esac
