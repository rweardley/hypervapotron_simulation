#!/bin/bash

source $HOME/rupert/.nekrs_next_241115_profile
nrsmpi hvt 8 2>&1 | tee log.run

