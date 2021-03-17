#!/bin/bash
#SBATCH --job-name=matlab_u_plotter
#SBATCH --output=matlab_u_plotter."%j".out
#SBATCH --error=matlab_u_plotter."%j".err
#SBATCH --partition=serc
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=ALL

module load matlab
module load system rclone

#matlab -nodisplay < u_plots.m
matlab -nodisplay < RMS_plots.m

rclone copy /home/users/mborrus/Matlab_HPC/plots remote:plot_folder
