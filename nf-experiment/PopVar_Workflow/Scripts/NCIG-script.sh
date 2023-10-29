#!/bin/bash

# This script is designed to be run within a PBS (Portable Batch System) job scheduler environment.
# It specifies resource requirements such as maximum wall time, number of CPUs, memory, etc.
# for the job to run on the cluster. Additionally, it loads the Nextflow module and then executes
# a Nextflow workflow with specific parameters.
#
# Key resource directives:
# - Wall time: The maximum time the job can run.
# - ncpus: Number of CPUs allocated for this job.
# - mem: Maximum memory allocated.
# - jobfs: Local disk space required.
# - wd: Start the job in the directory from which the `qsub` command was issued.
# - storage: Specifying data storage resources needed by the job.
# - q: Queue name.
# - P: Project identifier.
#
# After setting up the environment and requirements, the script loads the specified version of Nextflow
# and then runs a Nextflow workflow (main.nf) located in the workflows directory, with provided parameters.

# Authors: Kosar Hooshmand

#PBS -l walltime=48:00:00
#PBS -l ncpus=240
#PBS -l mem=190GB
#PBS -l jobfs=400GB
#PBS -l iointensive=10
#PBS -l wd
#PBS -l storage=scratch/te53+gdata/te53+gdata/xy86+gdata/if89
#PBS -q normal
#PBS -P te53
#PBS -N nextflow

module load nextflow/23.04.1

nextflow -log ./my.log -C ./nextflow_test.config run ./workflows/main-test.nf -profile NCI --postaltjs /g/data/te53/software/bwa/0.7.17-r1198-dirty/bwa-postalt.js --alt -resume <process_Name>

