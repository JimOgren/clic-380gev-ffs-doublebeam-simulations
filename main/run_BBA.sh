#!/bin/bash
# If running on lxplus or HTCondor, uncomment to setup placet etc.:
source /cvmfs/clicbp.cern.ch/x86_64-slc6-gcc62-opt/setup.sh

mkdir -p temp_run
cd temp_run
cp ../singlebeam_BBA.tcl .
placet singlebeam_BBA.tcl machine $1
