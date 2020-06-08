#!/bin/bash
# If running on lxplus or HTCondor, uncomment to setup placet etc.:
#source /cvmfs/clicbp.cern.ch/x86_64-slc6-gcc62-opt/setup.sh

mkdir -p temp_run
cd temp_run
cp ../doublebeam_tuning_step1.tcl .
cp ../run_guinea_parallel* .

placet doublebeam_tuning_step1.tcl machine $1
