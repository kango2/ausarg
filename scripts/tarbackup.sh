#!/bin/bash
#PBS -P xl04
#PBS -N tar
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=2:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+gdata/xl04+gdata/te53
#PBS -j oe

set -ex
set -o pipefail
set -u

usage() {
  echo "Usage: qsub -v FOLDER_PATH=/path/to/folder[,OUTDIR=/path/to/output] tar_folder.sh"
  exit 1
}

[[ "${1:-}" == "--help" ]] && usage

if [[ -z "${FOLDER_PATH:-}" ]]; then
  echo "[ERROR] FOLDER_PATH not set."
  usage
fi

if [[ ! -d "${FOLDER_PATH}" ]]; then
  echo "[ERROR] Input folder ${FOLDER_PATH} not found."
  exit 1
fi

OUTDIR="${OUTDIR:-$(pwd)}"
mkdir -p "${OUTDIR}"

BASENAME=$(basename "${FOLDER_PATH}")
DATE=$(date +%Y%m%d)
TARFILE="${OUTDIR}/${BASENAME}_${DATE}.tar.gz"

if [[ -f "${OUTDIR}/${BASENAME}_${DATE}.done" ]]; then
  echo "[SKIP] Archive already created: ${TARFILE}"
  exit 0
fi

touch "${OUTDIR}/${BASENAME}_${DATE}.running"

echo "[$(date)] 📦 Archiving ${FOLDER_PATH}..."
tar -czf "${TARFILE}" -C "$(dirname "${FOLDER_PATH}")" "${BASENAME}"
echo "[$(date)] ✅ Archive created at ${TARFILE}"

touch "${OUTDIR}/${BASENAME}_${DATE}.done"
rm -f "${OUTDIR}/${BASENAME}_${DATE}.running"
