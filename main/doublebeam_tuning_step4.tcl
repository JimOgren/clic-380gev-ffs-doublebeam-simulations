#############################################################################
#
# Tuning studies for double beam
#
# Step 4: incoherent pairs fast tuning knobs, machines below threshold
#
# Written by Jim Ogren, CERN, 2018-2020
#
#############################################################################
set t_1 [clock seconds]
global guinea_exec
set run_on_lxplus 0

if { $run_on_lxplus } {
   set script_dir /afs/cern.ch/work/j/jogren/TuningStudies/CLIC_FFS_380GeV/Clean_scripts/clic-380gev-ffs-doublebeam-simulations
   set guinea_exec /afs/cern.ch/eng/sl/clic-code/lx64slc5/guinea-pig/bin/guinea-old
} else {
   set script_dir /home/jim/GIT/clic-380gev-ffs-doublebeam-simulations
   set guinea_exec /home/jim/bin/guinea
}

set res_dir $script_dir/main/Results
set save_dir $res_dir/Results_step4

# Check if folders exist
if {![file exist $script_dir]} {
   puts "script_dir path does not exist!"
   exit
}

if {![file exist $res_dir]} {
   exec bash -c "mkdir -p $res_dir"
}
if {![file exist $save_dir]} {
   exec bash -c "mkdir -p $save_dir"
}

set e_initial 190
set e0 $e_initial

# command-line options
# sr          : synchrotron radiation on/off
# num_slice   : number of slices in multiparticle beam
# num_part    : number of particles per slice in multiparticle beam
# sigma       : misalignment in um
# sigmaM      : misalignment for multipoles in um
# sigmak      : relative magnetic strength errors
# sigmaroll   : roll misalignment in urad
# bpmres      : bpm resolution in um
# deltae      : energy difference for dfs
# machine     : machine seed
# loadmachine : load machine status from file on/off
# dfsmeasure  : to measure nominal dispersion and create model file
# wdisp       : gain for dispersion target steering (dts)
# iterdts     : number of dts iterations
# gaindts     : gain for dts
# beam_case   : which beam to load: 8, 20, 30 (vertical emittance)
# subset      : number of DOF used in random walk

array set args {
    sr 1
    num_slice 10
    num_part 2000
    sigma 10.0
    sigmaM 20.0
    sigmak 1e-4
    sigmaroll 100.0
    bpmres 0.020
    deltae 0.001
    machine 1
    loadmachine 0
    measure_response 0
    wdisp 0.71
    iterdts 30
    gaindts 0.5
    beam_case 20
	 subset 6
}

array set args $argv
set num_slice $args(num_slice)
set num_part $args(num_part)
set sigma $args(sigma)
set sigmaM $args(sigmaM)
set sigmak $args(sigmak)
set sigmaroll $args(sigmaroll)
set bpmres $args(bpmres)
set deltae $args(deltae)
set sr $args(sr)
set machine $args(machine)
set loadmachine $args(loadmachine)
set measure_response $args(measure_response)
set wdisp $args(wdisp)
set iterdts $args(iterdts)
set gaindts $args(gaindts)
set beam_case $args(beam_case)
set subset $args(subset)

#############################################################################
# load header files:
source $script_dir/scripts/make_beam.tcl
source $script_dir/scripts/wake_calc.tcl
source $script_dir/scripts/octave_functions.tcl
source $script_dir/scripts/beam_parameters.tcl

# load lattice and create beams
set beamline_list [list electron positron]
source $script_dir/scripts/loadlattice_double.tcl
source $script_dir/scripts/calc_lumi_jim.tcl

# Make two beams: 1e5 particles for checking lumi and 2e4 particles for faster tuning
source $script_dir/scripts/createbeams_double.tcl

FirstOrder 1

# Go through lattice file and find indices for different elements:
source $script_dir/scripts/element_indices.tcl

puts "Tuning with beamstrahlung - STEP 4 - Tuning with fast knobs"
puts "Machine = $machine"

source $script_dir/scripts/random_walk_methods.tcl
source $script_dir/scripts/setup_knobs_fast.tcl

#############################################################################
# Function for computing luminosity using the 1e5 particle beam
Octave {
	global IPind Sind Oind LL Lmeas
	IPind = IP;
	Sind = MIsext;
	Oind = MIoct;
	LL = [];
	Lmeas = 0;

	function [ ] = quick_check()
		global IPind LL Lmeas
		[E1,B1] = placet_test_no_correction("electron", "beam0s", "None", 1, 0, IPind);
	   [E2,B2] = placet_test_no_correction("positron", "beam0s", "None", 1, 0, IPind);
		[L, Lpeak, beam_size,beam_strahl] = get_lumi(B1, B2)
		e1 = beam_strahl(1);
		e2 = beam_strahl(2);
		sxp1 = beam_strahl(3);
		syp1 = beam_strahl(4);
		sxp2 = beam_strahl(5);
		syp2 = beam_strahl(6);
		printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
		printf("[e1, e2] \t= [%g, %g]\n",e1, e2);
		printf("[sxp1, syp1] \t= [%g, %g]\n",sxp1,syp1);
		printf("[sxp2, syp2] \t= [%g, %g]\n",sxp2,syp2);
		printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
		printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
		LL = [LL; L, Lpeak, beam_size, beam_strahl, Lmeas]
	end
}


# load machines after step 3
#############################################################################
load_beamline_status "electron" $res_dir/Results_step3/machine_status_electron_$machine.dat
load_beamline_status "positron" $res_dir/Results_step3/machine_status_positron_$machine.dat

SetBeamValues n $n_s n_slice $n_slice_s n_total $n_total_s
Octave { quick_check() }
SetBeamValues n $n_t n_slice $n_slice_t n_total $n_total_t

# Current luminosity
Octave {
	Tcl_Eval(["set Lumi ", num2str(LL(end,1))]);
}

# Sextupole tuning with knobs
#############################################################################
if { $Lumi < 8e33 } {
	puts "Tuning = True"
	sextupole_knobs_fast 15.0 $machine
}


SetBeamValues n $n_s n_slice $n_slice_s n_total $n_total_s
Octave { quick_check() }
SetBeamValues n $n_t n_slice $n_slice_t n_total $n_total_t

# Save data
#############################################################################
Octave {
	printf("\nSaving data...\n");
	fname = ["$save_dir", "/tuning_data_",num2str($machine),".dat"]
	eval(["save ", fname," LL "]);

	Lumis = load("LuminosityData.dat")	;
	fname = ["$save_dir", "/LuminosityData_" ,num2str($machine),".dat"]
	eval(["save ", fname," Lumis"]);

	STATUS = placet_element_get_attributes("electron");
	eval(["save -text ", "$save_dir", "/machine_status_electron_$machine.dat STATUS"]);
	STATUS = placet_element_get_attributes("positron");
	eval(["save -text ", "$save_dir", "/machine_status_positron_$machine.dat STATUS"]);
}
