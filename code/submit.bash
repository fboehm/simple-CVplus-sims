#!/bin/bash
#
#SBATCH --job-name=jkplus
#SBATCH --partition=mulan,main
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

source ~/miniconda3/etc/profile.d/conda.sh
# https://github.com/conda/conda/issues/7980

conda activate barber1
quarto render barber-JK+.qmd
conda deactivate



