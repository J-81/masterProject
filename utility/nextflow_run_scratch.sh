#!/bin/bash

#SBATCH -p nodes
#SBATCH --job-name="nf_wf_manage"
#SBATCH -o slurm_logs/std_output_%j.out
#SBATCH -e slurm_logs/std_error_%j.err
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=23:59:59
#SBATCH -D /home/joribello/mP/masterProject
# #SBATCH --test-only

export TOWER_ACCESS_TOKEN=7e38ad64660cf92608d67a9425562bc55722d54c
export NXF_VER=20.07.1

nextflow run main.nf -c config/default.config -with-tower -w '../workdir/GLDS-104/work' -resume -profile cos_hpc_4_node
