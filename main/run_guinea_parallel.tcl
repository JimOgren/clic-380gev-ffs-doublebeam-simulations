# input offset in nm
ParallelThreads -num 1
array set args {
    thick_beam 0
    pairs 0
    beam_case 0
}

array set args $argv
set thick_beam $args(thick_beam)
set pairs $args(pairs)
set beam_case $args(beam_case)

puts " "
puts "Computing luminosity"
puts "thick_beam $thick_beam"
puts "pairs $pairs"
puts "beam_case $beam_case"


global guinea_exec
set run_on_lxplus 0
if { $run_on_lxplus } {
   set script_dir /afs/cern.ch/work/j/jogren/TuningStudies/CLIC_FFS_380GeV/Clean_scripts/clic-380gev-ffs-doublebeam-simulations
   set guinea_exec /afs/cern.ch/eng/sl/clic-code/lx64slc5/guinea-pig/bin/guinea-old
} else {
   set script_dir /home/jim/GIT/clic-380gev-ffs-doublebeam-simulations
   set guinea_exec /home/jim/bin/guinea
}

set e0 190
set match(charge) 5.2e9
set match(sigma_z) 70.0

if { $thick_beam == 1 } {
   set n_slice 50
   set n 2000
   set n_total [expr $n_slice*$n]
} else {
   set n_slice 20
   set n 1000
   set n_total [expr $n_slice*$n]
}

source $script_dir/scripts/calc_lumi_jim.tcl

Octave {
	str1 = "B1 = load(""../beams_temp/electron_${beam_case}.dat"");";
	str2 = "B2 = load(""../beams_temp/positron_${beam_case}.dat"");";
   eval(str1);
   eval(str2);

   if $pairs == 1
      [L, Lpeak, beam_size, beam_strahl, beam_defl, pairs] = get_lumi_full(B1, B2);
      save ../GP_output_temp/GP_output_${beam_case}.dat L Lpeak beam_size beam_strahl beam_defl pairs
   else
      [L, Lpeak, beam_size, beam_strahl, beam_defl] = get_lumi(B1, B2);
      save ../GP_output_temp/GP_output_${beam_case}.dat L Lpeak beam_size beam_strahl beam_defl
   end
}
