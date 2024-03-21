#!/bin/bash -e

# When pointed at directory prepared by 'prepRun' will launch matlab wrapper there.
#
# Arguments:
#   SP_WORKDIR: Directory containing decomposed job.
#   SP_MLVER(optional): Version of matlab for workers to use.

cat $0

# Assign args to env vars.
export SP_WORKDIR=$1
export SP_MLVER=$2

# Worker ID will be slurm array task ID unless otherwise set here.
export SP_WORKER_ID="${SLURM_ARRAY_TASK_ID}"

# To avoid being contaminated by base env.
export SLURM_EXPORT_ENV=ALL

# Load specific version of matlab if requested.
module load MATLAB/${SP_MLVER}
matlab -nodisplay -nosplash <slurm_parallel/matlabWrapper.m
