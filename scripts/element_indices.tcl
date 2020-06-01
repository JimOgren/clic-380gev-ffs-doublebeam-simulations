Octave {
	IP = placet_get_name_number_list("electron", "IP");
	BI = placet_get_number_list("electron", "bpm");

	DI = placet_get_number_list("electron", "dipole");
	QI = placet_get_number_list("electron", "quadrupole");

	MI = placet_get_number_list("electron", "multipole");
	MS = placet_element_get_attribute("electron", MI, "strength");

	MI = MI(MS != 0);
	MS = placet_element_get_attribute("electron", MI, "strength");
	ML = placet_element_get_attribute("electron", MI, "length");
	MT = placet_element_get_attribute("electron", MI, "type");

	T = placet_element_get_attribute("electron", MI, "type");

	MIsext = MI(T<4);
	MSsext = placet_element_get_attribute("electron", MIsext, "strength");

	MIoct = MI(T==4);
	MSoct = placet_element_get_attribute("electron", MIoct, "strength");
  
	SI = placet_get_number_list("electron", "sbend");
	AS = placet_element_get_attribute("electron", SI, "angle");
}
