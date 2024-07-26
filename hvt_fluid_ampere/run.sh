#!/bin/bash

source ~/.nekrs_23-0_profile

cores=16

ulimit -s unlimited

which nrsmpi | tee log.nrsversion

nrsmpi hvt $cores 2>&1 | tee log.hvt
nrsvis hvt
