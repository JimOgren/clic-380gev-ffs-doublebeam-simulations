source $script_dir/scripts/clic_basic_single.tcl
source $script_dir/scripts/clic_beam.tcl

set scale 1.0

set synrad $sr
set quad_synrad $sr
set mult_synrad $sr
set sbend_synrad $sr

puts "Load lattice"

SetReferenceEnergy $e0

foreach name $beamline_list {
   BeamlineNew
   source $script_dir/scripts/lattices/ffs_clic380gev_l6m_bx8_disp60.tcl

   Bpm -name "IP"
   Drift -name "D0" -length 6.016752

   Bpm
   BeamlineSet -name $name

   # Slice the sbends for more accurate sr in tracking
   Octave {
      SI = placet_get_number_list("$name", "sbend");
      placet_element_set_attribute("$name", SI, "thin_lens", int32(100));
   }
}
