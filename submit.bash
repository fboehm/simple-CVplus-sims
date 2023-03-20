#!/bin/bash
#
#SBATCH --job-name=jk+
#SBATCH --output=out.txt
#SBATCH --partition=debug
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

conda activate barber1
quarto render barber-JK+.qmd
conda deactivate



