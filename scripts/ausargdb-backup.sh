#!/bin/bash
#PBS -N sqlite_backup
#PBS -P xl04
#PBS -q copyq
#PBS -l walltime=10:00:00,ncpus=1,mem=4GB
#PBS -j oe
#PBS -l storage=gdata/xl04+massdata/xl04
#PBS -m ae


set -ex

db_dir=$(dirname "$filepath")
db_file=$(basename "$filepath")

#PBS -o ${db_dir}/logs

# Path to the backup mdss directory
backup_dir="sqliteDB_backup"

cd $db_dir

timestamp=$(date +%Y%m%d%H%M%S)

target="$backup_dir/backup_${timestamp}_${db_file}"

echo "$db_dir, $db_file"

mdss -P xl04 put $db_file $target

if [ $? -eq 0 ]; then
    echo "${db_file} backup created at ${timestamp}: $target" >> "${db_dir}/logs/backlog.log"
else
    echo "${db_file} backup failed at ${timestamp}" >> "${db_dir}/logs/backlog.log"
fi

TIMES=("0400" "1400" "1900")

# Loop through the times and submit jobs
for TIME in "${TIMES[@]}"; do
        qsub -o /g/data/xl04/genomeprojects/database/logs/PBS -a ${TIME} -v filepath="$filepath" /g/data/xl04/ka6418/github/ausarg/scripts/ausargdb-backup.sh
done


