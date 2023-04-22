#!/bin/bash
#
#SBATCH --job-name=jkplus
#SBATCH --partition=mulan,main
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

quarto render simple2.qmd --to html,gfm


