#!/usr/bin/env bash
_hera_grid_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export TMPDIR="${_hera_grid_root}/.tmp" TMP="${_hera_grid_root}/.tmp" TEMP="${_hera_grid_root}/.tmp"
export MPLCONFIGDIR="${TMPDIR}/matplotlib" XDG_CACHE_HOME="${TMPDIR}/cache"
mkdir -p "$TMPDIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
unset _hera_grid_root
