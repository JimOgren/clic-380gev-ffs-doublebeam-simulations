# Make two beams: 1e5 particles for checking lumi and 2e4 particles for faster tuning
####################
# Procedure for chaning the global variables of current beam
proc SetBeamValues { var1 v1 var2 v2 var3 v3 } {
	upvar $var1 myvar1
	upvar $var2 myvar2
	upvar $var3 myvar3
	set myvar1 $v1
	set myvar2 $v2
	set myvar3 $v3	
}

# create the two beams
puts "Creating beams"
set e0 $e_initial
set e2 [expr $e0 * (1.0 + $deltae)]
set e1 [expr $e0 * (1.0 - $deltae)]

set n_slice 1
set n 1
set n_total [expr $n_slice*$n]

make_beam_many beam0 $n_slice $n
exec echo "$e0 0 0 0 0 0" > particles.in
BeamRead -file particles.in -beam beam0
make_beam_many beam1 $n_slice $n
exec echo "$e1 0 0 0 0 0" > particles.in
BeamRead -file particles.in -beam beam1
make_beam_many beam2 $n_slice $n
exec echo "$e2 0 0 0 0 0" > particles.in
BeamRead -file particles.in -beam beam2

set n_slice_s 50
set n_s 2000
set n_total_s [expr $n_slice_s*$n_s]
SetBeamValues n $n_s n_slice $n_slice_s n_total $n_total_s
puts "Creating multiparticle beam with $n_total particles"
make_beam_many beam0s $n_slice_s $n_s

set n_slice_t 20
set n_t 1000
set n_total_t [expr $n_slice_t*$n_t]
SetBeamValues n $n_t n_slice $n_slice_t n_total $n_total_t
puts "Creating multiparticle beam with $n_total particles"
make_beam_many beam0t $n_slice_t $n_t
#make_beam_many beam1t $n_slice_t $n_t
#make_beam_many beam2t $n_slice_t $n_t

# Load beams from integrated simulation
Octave {
   disp("\nLoading beam from integrated simulation")
   disp('Loading beam with ey = 20 nm (end of main linac)')
   load $script_dir/beams/beam_start_of_FFS_ey20nm_1e5.dat
   BeamIn = beam;
   placet_set_beam("beam0s", BeamIn);
   
   load $script_dir/beams/beam_start_of_FFS_ey20nm_2e4.dat
   BeamIn = beam;
   placet_set_beam("beam0t", BeamIn);
   
#   BeamIn1 = BeamIn;
#   BeamIn2 = BeamIn;
#   BeamIn1(:,1) = BeamIn1(:,1)*(1.0+$deltae);
#   BeamIn2(:,1) = BeamIn2(:,1)*(1.0-$deltae);   
#   placet_set_beam("beam1t", BeamIn1);
#   placet_set_beam("beam2t", BeamIn2);
   disp(' ')
}
