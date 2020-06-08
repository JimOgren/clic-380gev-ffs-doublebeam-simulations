#############################################################################
#
# Tuning studies for double beam
#
# Step 0: beam-based alignment
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
set save_dir $res_dir/Results_BBA

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
#############################################################################
# command-line options
# sr          : synchrotron radiation on/off
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

array set args {
	sr 1
	sigma 10.0
	sigmaM 20.0
	sigmak 1e-4
	sigmaroll 100.0
	bpmres 0.020
	deltae 0.001
	machine 1
	loadmachine 0
	measure_response 1
	wdisp 0.71
	iterdts 30
	gaindts 0.5
	beam_case 20
}

array set args $argv
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

#############################################################################
# load header files:
source $script_dir/scripts/make_beam.tcl
source $script_dir/scripts/wake_calc.tcl
source $script_dir/scripts/octave_functions.tcl
source $script_dir/scripts/beam_parameters.tcl
source $script_dir/scripts/checkDispersionAndAlignment.tcl

# load lattice, create beamline "electron" and "positron"
set beamline_list [list electron positron]
source $script_dir/scripts/loadlattice_double.tcl

# functions for computing luminosity
source $script_dir/scripts/calc_lumi_jim.tcl

# make two beams: 1e5 particles for checking lumi and 2e4 particles for faster tuning
source $script_dir/scripts/createbeams_double.tcl

# go through lattice file and find indices for different elements:
source $script_dir/scripts/element_indices.tcl

# procedures for beam-based alignment
source $script_dir/scripts/beam_based_alignment_methods.tcl

# measure target dispersion
if { $measure_response } {
   source $script_dir/scripts/measure_nominal_dispersion_single_beam.tcl
}
# load target dispersion
source $script_dir/scripts/load_machine_model_single_beam.tcl

FirstOrder 1

puts "Step 0 - misalign machines and perform BBA"
puts "Machine = $machine\n"

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
		[L, Lpeak, beam_size,beam_strahl] = get_lumi(B1, B2);
		e1 = beam_strahl(1);
		e2 = beam_strahl(2);
		sxp1 = beam_strahl(3);
		syp1 = beam_strahl(4);
		sxp2 = beam_strahl(5);
		syp2 = beam_strahl(6);
		printf("\n[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
		printf("[e1, e2] \t= [%g, %g]\n",e1, e2);
		printf("[sxp1, syp1] \t= [%g, %g]\n",sxp1,syp1);
		printf("[sxp2, syp2] \t= [%g, %g]\n",sxp2,syp2);
		printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
		printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
		LL = [LL; L, Lpeak, beam_size, beam_strahl, Lmeas]
	end
}

# Procedure for misaligning a machine
proc my_survey { beamline_name } {
   global machine sigma sigmaM sigmaroll sigmak bpmres
   puts "Adding random static imperfections to $beamline_name beamline"
   Octave {
      randn("seed", $machine * 234);
      placet_element_set_attribute("$beamline_name", DI, "strength_x", 0.0);
      placet_element_set_attribute("$beamline_name", DI, "strength_y", 0.0);
      placet_element_set_attribute("$beamline_name", BI, "resolution", $bpmres);
      placet_element_set_attribute("$beamline_name", BI, "x", randn(size(BI)) * $sigma);
      placet_element_set_attribute("$beamline_name", BI, "y", randn(size(BI)) * $sigma);
      placet_element_set_attribute("$beamline_name", QI, "x", randn(size(QI)) * $sigma);
      placet_element_set_attribute("$beamline_name", QI, "y", randn(size(QI)) * $sigma);
      placet_element_set_attribute("$beamline_name", MI, "x", randn(size(MI)) * $sigmaM);
      placet_element_set_attribute("$beamline_name", MI, "y", randn(size(MI)) * $sigmaM);
      placet_element_set_attribute("$beamline_name", DI, "roll",randn(size(DI)) * $sigmaroll);
      placet_element_set_attribute("$beamline_name", QI, "roll",randn(size(QI)) * $sigmaroll);
      placet_element_set_attribute("$beamline_name", BI, "roll",randn(size(BI)) * $sigmaroll);
      placet_element_set_attribute("$beamline_name", MI, "roll",randn(size(MI)) * $sigmaroll);
      placet_element_set_attribute("$beamline_name", SI, "roll",randn(size(SI)) * $sigmaroll);

      QIK=placet_element_get_attribute("$beamline_name", QI, "strength");
      QIKE=QIK.*(1 + randn(size(QIK)) * $sigmak);
      MIK=placet_element_get_attribute("$beamline_name", MI, "strength");
      MIKE=MIK.*(1 + randn(size(MIK)) * $sigmak);
      ANG = placet_element_get_attribute("$beamline_name", SI, "angle");
		ANGE = ANG.*(1 + randn(size(ANG)) * $sigmak );


      placet_element_set_attribute("$beamline_name", QI, "strength", QIKE);
      placet_element_set_attribute("$beamline_name", MI, "strength", MIKE);
      placet_element_set_attribute("$beamline_name", SI, "angle", ANGE);
   }
}

# Step 0 Check luminosity of perfect machine
#############################################################################
SetBeamValues n $n_s n_slice $n_slice_s n_total $n_total_s
Octave { quick_check() }


# Step 1 misalign the machines
#############################################################################
my_survey "electron"
# to get different seed for positron beamline:
set machine [expr $machine +1000]
my_survey "positron"
set machine [expr $machine -1000]

# check luminosity of misaligned machine
Octave { quick_check() }
checkDispersionAndAlignment

# Step 2 perform beam-based alignment
#############################################################################
run_BBA_single_beam "electron"
run_BBA_single_beam "positron"

SetBeamValues n $n_s n_slice $n_slice_s n_total $n_total_s
Octave { quick_check() }
checkDispersionAndAlignment

# Save results and machine status
#############################################################################
Octave {
	printf("\nSaving data...\n");
	fname = ["$save_dir", "/tuning_data_",num2str($machine),".dat"]
	eval(["save ", fname," LL"]);

	STATUS = placet_element_get_attributes("electron");
	eval(["save -text ", "$save_dir", "/machine_status_electron_$machine.dat STATUS"]);
	STATUS = placet_element_get_attributes("positron");
	eval(["save -text ", "$save_dir", "/machine_status_positron_$machine.dat STATUS"]);
}
