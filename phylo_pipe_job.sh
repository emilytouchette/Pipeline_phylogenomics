#!/bin/bash
#SBATCH --job-name=phylo_pipe_job
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:0:0
#SBATCH --output=phylo_pipe_job.out
#SBATCH --error=phylo_pipe_job.err
#SBATCH --mail-user=touchete@oregonstate.edu
#SBATCH --mail-type=END


#Bash commands to run within job
python phylogenomic_pipeline.py