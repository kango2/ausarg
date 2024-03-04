#!/bin/bash
#PBS -N downloading
#PBS -P xl04
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l mem=2GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -j oe
#PBS -M kirat.alreja@anu.edu.au

export CKAN_API_TOKEN="61beb2cd-aad9-4b96-b9e5-051fa2e6428e"
bash ${script}