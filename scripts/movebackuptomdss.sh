#!/bin/bash
#PBS -P xl04
#PBS -N move_to_mdss
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/if89+gdata/xl04+gdata/te53
#PBS -j oe
#PBS -l wd

set -ex
set -o pipefail
set -u

usage() {
  echo ""
  echo "Usage: qsub -v ARCHIVE=/path/to/file.tar.gz[,MDSS_DIR=/mdss/path] move_to_mdss.sh"
  echo ""
  echo "Required:"
  echo "  ARCHIVE    Full path to the archive file"
  echo "Optional:"
  echo "  MDSS_DIR   Target directory on MDSS (default: /mdss/te53/archive)"
  exit 1
}

[[ "${1:-}" == "--help" ]] && usage

if [[ -z "${ARCHIVE:-}" ]]; then
  echo "[ERROR] ARCHIVE variable is required."
  usage
fi

if [[ ! -f "${ARCHIVE}" ]]; then
  echo "[ERROR] Archive file ${ARCHIVE} does not exist."
  exit 1
fi

MDSS_DIR="${MDSS_DIR:-/mdss/te53/archive}"
ARCHIVE_NAME=$(basename "${ARCHIVE}")
DONE_FLAG="${ARCHIVE}.mdss.done"

if [[ -f "${DONE_FLAG}" ]]; then
  echo "[SKIP] Archive already moved to MDSS: ${ARCHIVE_NAME}"
  exit 0
fi

touch "${ARCHIVE}.mdss.running"

module load mdss/1.2.11

echo "[$(date)] 📤 Moving ${ARCHIVE} to ${MDSS_DIR} on MDSS..."

mdss mkdir -p "${MDSS_DIR}"
mdss put -P xl04 "${ARCHIVE}" "${MDSS_DIR}/"

echo "[$(date)] ✅ Transfer complete: ${MDSS_DIR}/${ARCHIVE_NAME}"

touch "${DONE_FLAG}"
rm -f "${ARCHIVE}.mdss.running"
