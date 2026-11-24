#!/bin/bash
# Fetch the LoTSS L234026 "big" dataset for a realistic, non-CI run:
# mslist.txt, big-mslist.txt, examples/lotss_mpi.cfg and the shared CATALOGS.
# Every download is SHA256-verified. The Measurement Sets themselves are a
# multi-GB folder share and are NOT downloaded here (see the note at the end).
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

MSLIST_URL="https://sdrive.cnrs.fr/public.php/dav/files/YjPXwZJ62FwX6Jo"
MSLIST_SHA256="8fdadfc4a9e1cf9e85cd621f2a3883cc57560beb0934e2b19ad155c29ce06397"

BIGMSLIST_URL="https://sdrive.cnrs.fr/public.php/dav/files/NRHNdrLGZs6o7a4"
BIGMSLIST_SHA256="33d6b67cdab8ee17f50def90ac5ba936999ff9685d48a7341874912072044a5c"

LOTSS_CFG_URL="https://sdrive.cnrs.fr/public.php/dav/files/omnAgeZL5cPHAAn"
LOTSS_CFG_SHA256="e0d0a16ad0cfd39183ab9d4c70d0cce8fae501e220cb0dd7f960f4bf3ee3ecfe"

CATALOGS_URL="https://sdrive.cnrs.fr/public.php/dav/files/awWCRsaoeax8t4K/CATALOGS.zip"
CATALOGS_SHA256="8a0edbebc67613266251f9871b09c588c1c8102368daf1161045ecc26d4d701d"

require_pinned mslist.txt "$MSLIST_SHA256"
require_pinned big-mslist.txt "$BIGMSLIST_SHA256"
require_pinned lotss_mpi.cfg "$LOTSS_CFG_SHA256"
require_pinned CATALOGS.zip "$CATALOGS_SHA256"

fetch_pinned "$MSLIST_URL" mslist.txt "$MSLIST_SHA256"
fetch_pinned "$BIGMSLIST_URL" big-mslist.txt "$BIGMSLIST_SHA256"
fetch_pinned "$LOTSS_CFG_URL" examples/lotss_mpi.cfg "$LOTSS_CFG_SHA256"
fetch_pinned "$CATALOGS_URL" CATALOGS.zip "$CATALOGS_SHA256"
unzip -o -q CATALOGS.zip

cat <<'NOTE'
LoTSS "big": mslist.txt, big-mslist.txt, examples/lotss_mpi.cfg and CATALOGS/ installed.
The Measurement Sets are a multi-GB folder share: download them from
  https://sdrive.cnrs.fr/s/azWk38oXb6BwrRS
on WORK, then point the mslist at them. This dataset is for a realistic run, not CI.
NOTE
