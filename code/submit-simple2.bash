#!/bin/bash
#
#SBATCH --job-name=simple2
#SBATCH --partition=mulan,main
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/net/mulan/home/fredboe/research/simple-CVplus-sims/cluster_outputs/simple2.out
#SBATCH --error=/net/mulan/home/fredboe/research/simple-CVplus-sims/cluster_outputs/simple2.err

conda run -n barber1 quarto render simple2.qmd --to html


