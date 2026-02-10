#!/usr/bin/env bash
set -euo pipefail
# run_rx_sex_id.sh
#
# Thin wrapper around rx_sex_id.py so it can be called like a normal CLI.
#
# Examples:
#   ./run_rx_sex_id.sh *bam /path/ref.fna results.tsv
#   ./run_rx_sex_id.sh --x-contig NC_079873.1 --y-contig NC_079874.1 *.bam ref.fna out.tsv
#
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SCRIPT_DIR}/rx_sex_id.py"
exec python3 "$PY" "$@"
