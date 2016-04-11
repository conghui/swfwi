#!/bin/bash

# it is critical important to set this environment variable
# in order to keep everything in the corrent folder
export DATAPATH=`pwd`/

# only works in GNU or INTEL compiler
export OMP_NUM_THREADS=8

# possible option:
# fm,     essfwi,     enfwi,    noise
# fm-sw,  essfwi-sw,  enfwi-sw, noise-sw
#
scons task=fm
