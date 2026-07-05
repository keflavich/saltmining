#!/bin/bash
# Activate Bhula's XCLASS env and run a fit against one of our extracted spectra.
#
# Usage:
#   bash run_xclass.sh G326.6618+00.5207 [spw_filename_stem]
#
# Activates the python env that has task_myXCLASS importable, then
# runs the per-spw runner inside the target's lineid subdirectory.
set -e

TARGET="${1:?need target name (e.g. G326.6618+00.5207)}"
LABEL="${2:-${TARGET}}"  # subdir label under lineid/

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TGT_DIR="${HERE}/${LABEL}"
[ -d "${TGT_DIR}" ] || { echo "no dir ${TGT_DIR}"; exit 1; }

source /orange/adamginsburg/miniconda3/bin/activate /blue/adamginsburg/adamginsburg/miniconda3/envs/xclass

# Make libgfortran.so.3 discoverable for myNewXCLASS.exe (Bhula's setup
# symlinks /apps/libgfortran/3.0.0/lib/libgfortran.so.3 into his work dir;
# we just prepend the dir to LD_LIBRARY_PATH).
export LD_LIBRARY_PATH="/apps/libgfortran/3.0.0/lib:${LD_LIBRARY_PATH}"

cd "${TGT_DIR}"
python "${HERE}/xclass_runner.py" "${LABEL}"
