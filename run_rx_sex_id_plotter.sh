#!/usr/bin/env bash
set -euo pipefail
# run_rx_sex_id_plotter.sh
#
# Usage:
#   ./run_rx_sex_id_plotter.sh out_prefix Label=/path/file.tsv [Label2=/path/file2.tsv ...]
#
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/rx_sex_id_plotter.R"
Rscript "$R_SCRIPT" --out-prefix "$1" --inputs "${@:2}"
