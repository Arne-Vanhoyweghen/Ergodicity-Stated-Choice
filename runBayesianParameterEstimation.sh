#!/bin/bash
#SBATCH --partition=HPC
matlab -nodesktop -nojvm -nosplash -r "runBayesianParameterEstimation"
