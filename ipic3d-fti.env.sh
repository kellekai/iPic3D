#!/bin/bash
module purge
export LD_LIBRARY_PATH=""
export LIBRARY_PATH=""
export PATH=""
module load gcc/7.2.0 MPICH3 H5hut HDF5
export LD_LIBRARY_PATH=/home/kellekai/WORK/DEEP-EST/D6.1/apps/FTI-0.9.6/RELEASE/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/kellekai/WORK/DEEP-EST/D6.1/apps/FTI-0.9.6/RELEASE/lib:$LIBRARY_PATH

