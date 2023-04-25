#!/bin/bash
#
#SBATCH --job-name=jkplus
#SBATCH --partition=mulan,main
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G


conda run -n barber1 quarto render simple2.qmd --to html


