#!/bin/bash

source ~/.nekrs_next_241127_1121_profile

export NEKRS_GPU_MPI=1
export UCX_RNDV_THRESH=1024
export UCX_RC_TM_ENABLE=y

# MPI is set up to use slurm to check allocation
# so this will run 1 rank per node allocated
mpirun -N 1 nekrs --setup hvt 1 2>&1 | tee log.run
