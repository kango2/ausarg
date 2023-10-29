# Tailored for environments that use the PBS scheduler, ensuring that the tasks are efficiently queued and executed within the PBS system.

#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=63GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o md5sum.log

module load python3
cd $PBS_O_WORKDIR
python3 ncig-md5-updater.py
python3 ncig-ont-table-updater.py
python3 ncig-ont-sex-updater.py
python3 ncig-illumina-table-updater.py
python3 ncig-illumina-sex-updater.py
