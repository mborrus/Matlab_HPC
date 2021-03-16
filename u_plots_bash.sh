#!/bin/bash
#SBATCH --job-name=matlab_u_plotter
#SBATCH --output=matlab_u_plotter."%j".out
#SBATCH --error=matlab_u_plotter."%j".err
#SBATCH --partition=normal
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL

module load matlab
matlab -nodisplay < u_plots.m