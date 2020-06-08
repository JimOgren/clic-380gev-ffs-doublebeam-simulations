#!/bin/sh
mkdir -p temp_$3
cd temp_$3
cp ../run_guinea_parallel.tcl .
placet run_guinea_parallel.tcl thick_beam $1 pairs $2 beam_case $3
cd ../
rm -rf temp_$3
