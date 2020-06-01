proc checkDispersionAndAlignment {} {
	global deltae
	Octave {
		if (~exist("Evec"))
			# initiate vectors
			Evec = [];
			misQvec = [];
			misMvec = [];
			misSvec = [];
			SbendRollvec = [];
			wrtbQvec = [];
			wrtbMvec = [];
			multStreng = [];
			angleSbends = [];
		end

		disp(' ')
		disp('==========================================================================')
		disp(' ')
		disp('Dispersion target WITHOUT multipoles - ELECTRON')
		MStemp = placet_element_get_attribute("electron", MI, "strength");
		placet_element_set_attribute("electron", MI, "strength", 0.0);
		eta_meas_arr = CheckDispersion("electron", "beam1", "beam2", $deltae, BI);
		eta_meas = [eta_meas_arr(:,1); eta_meas_arr(:,2)];
		[eta0_arr, eta_meas_arr, abs(eta0_arr - eta_meas_arr)]
		stdE = [sqrt(sum((eta_meas_arr(:,1) - eta0_arr(:,1)).^2)), ...
				sqrt(sum((eta_meas_arr(:,2) - eta0_arr(:,2)).^2))]
		placet_element_set_attribute("electron", MI, "strength", MStemp);
	
		disp('Dispersion target WITHOUT multipoles - POSITRON')
		MStemp = placet_element_get_attribute("positron", MI, "strength");
		placet_element_set_attribute("positron", MI, "strength", 0.0);
		eta_meas_arr = CheckDispersion("positron", "beam1", "beam2", $deltae, BI);
		eta_meas = [eta_meas_arr(:,1); eta_meas_arr(:,2)];
		[eta0_arr, eta_meas_arr, abs(eta0_arr - eta_meas_arr)]
		stdEP = [sqrt(sum((eta_meas_arr(:,1) - eta0_arr(:,1)).^2)), ...
				sqrt(sum((eta_meas_arr(:,2) - eta0_arr(:,2)).^2))]
		placet_element_set_attribute("positron", MI, "strength", MStemp);

		Evec = [Evec; stdE, stdEP];

		BPM_RES = placet_element_get_attribute("electron", BI(1), "resolution")
		BPM_RES = placet_element_get_attribute("positron", BI(1), "resolution")

		disp(' ')
		disp('Element misalignments')
		disp('Quadrupoles')
		xc = placet_element_get_attribute("electron", QI, "x");
		yc = placet_element_get_attribute("electron", QI, "y");
		rc = placet_element_get_attribute("electron", QI, "roll");
		xc2 = placet_element_get_attribute("positron", QI, "x");
		yc2 = placet_element_get_attribute("positron", QI, "y");
		rc2 = placet_element_get_attribute("positron", QI, "roll");
		disp('x_pos - y_pos - roll - x_pos - y_pos - roll')
		[xc, yc, rc, xc2, yc2, rc2]
		stdQ = [sqrt(sum(xc.^2)), sqrt(sum(yc.^2)), sqrt(sum(rc.^2)), sqrt(sum(xc2.^2)), sqrt(sum(yc2.^2)), sqrt(sum(rc2.^2))]
		misQvec = [misQvec; stdQ];

		disp(' ')
		disp('Multipoles')
		xc = placet_element_get_attribute("electron", MI, "x");
		yc = placet_element_get_attribute("electron", MI, "y");
		rc = placet_element_get_attribute("electron", MI, "roll");
		xc2 = placet_element_get_attribute("positron", MI, "x");
		yc2 = placet_element_get_attribute("positron", MI, "y");
		rc2 = placet_element_get_attribute("positron", MI, "roll");
		disp('x_pos - y_pos - roll - x_pos - y_pos - roll')
		[xc, yc, rc, xc2, yc2, rc2]
		stdM = [sqrt(sum(xc.^2)), sqrt(sum(yc.^2)), sqrt(sum(rc.^2)), sqrt(sum(xc2.^2)), sqrt(sum(yc2.^2)), sqrt(sum(rc2.^2))]
		misMvec = [misMvec; stdM];
	
		disp(' ')
		disp('Sbends')
		rc = placet_element_get_attribute("electron", SI, "roll");
		rc2 = placet_element_get_attribute("positron", SI, "roll");
		[rc, rc2]
		stdM = [sqrt(sum(rc.^2)), sqrt(sum(rc2.^2))]
		SbendRollvec = [SbendRollvec; stdM];

		disp(' ')
		disp('Misalignment w r t the beam')
		disp('Quadrupoles')
		[stdx, stdy] = checkAlignment("electron", "beam0", QI);
		stdQwrt = [stdx, stdy]
		[stdx, stdy] = checkAlignment("positron", "beam0", QI);
		stdQwrtP = [stdx, stdy]
		wrtbQvec = [wrtbQvec; stdQwrt, stdQwrtP];

		disp('Multipoles')
		[stdx, stdy] = checkAlignment("electron", "beam0", MI);
		stdMwrt = [stdx, stdy]
		[stdx, stdy] = checkAlignment("positron", "beam0", MI);
		stdMwrtP = [stdx, stdy]
		wrtbMvec = [wrtbMvec; stdMwrt, stdMwrtP];

	
		disp('Multipole strengths')
		MStemp = placet_element_get_attribute("electron", MI, "strength");
		[MStemp, MS, MStemp-MS]
		stdMS = sqrt(sum((MStemp-MS).^2))
		MStemp = placet_element_get_attribute("positron", MI, "strength");
		[MStemp, MS, MStemp-MS]
		stdMSP = sqrt(sum((MStemp-MS).^2))
		multStreng = [multStreng; stdMS, stdMSP];

	  	disp('SBend angles')
	  	AStemp = placet_element_get_attribute("electron", SI, "angle");
	  	stdAS = sqrt(sum(((AStemp-AS)./AS).^2));
	  	AStemp = placet_element_get_attribute("positron", SI, "angle");
	  	stdASP = sqrt(sum(((AStemp-AS)./AS).^2));
	  	
	  	angleSbends = [angleSbends; stdAS, stdASP];


		Evec
		misQvec
		misMvec
		wrtbQvec
		wrtbMvec
		multStreng
		angleSbends
	
		disp('==========================================================================')
		disp(' ')
	}
}


proc checkDispersionAndAlignment_2e4 {} {
	global deltae
	Octave {
		if (~exist("Evec"))
			# initiate vectors
			Evec = [];
			misQvec = [];
			misMvec = [];
			misSvec = [];
			SbendRollvec = [];
			wrtbQvec = [];
			wrtbMvec = [];
			multStreng = [];
			angleSbends = [];
		end

		disp(' ')
		disp('==========================================================================')
		disp(' ')
		disp('Dispersion target WITHOUT multipoles - ELECTRON')
		MStemp = placet_element_get_attribute("electron", MI, "strength");
		placet_element_set_attribute("electron", MI, "strength", 0.0);
		eta_meas_arr = CheckDispersion("electron", "beam1t", "beam2t", $deltae, BI);
		eta_meas = [eta_meas_arr(:,1); eta_meas_arr(:,2)];
		[eta0_arr, eta_meas_arr, abs(eta0_arr - eta_meas_arr)]
		stdE = [sqrt(sum((eta_meas_arr(:,1) - eta0_arr(:,1)).^2)), ...
				sqrt(sum((eta_meas_arr(:,2) - eta0_arr(:,2)).^2))]
		placet_element_set_attribute("electron", MI, "strength", MStemp);
	
		disp('Dispersion target WITHOUT multipoles - POSITRON')
		MStemp = placet_element_get_attribute("positron", MI, "strength");
		placet_element_set_attribute("positron", MI, "strength", 0.0);
		eta_meas_arr = CheckDispersion("positron", "beam1t", "beam2t", $deltae, BI);
		eta_meas = [eta_meas_arr(:,1); eta_meas_arr(:,2)];
		[eta0_arr, eta_meas_arr, abs(eta0_arr - eta_meas_arr)]
		stdEP = [sqrt(sum((eta_meas_arr(:,1) - eta0_arr(:,1)).^2)), ...
				sqrt(sum((eta_meas_arr(:,2) - eta0_arr(:,2)).^2))]
		placet_element_set_attribute("positron", MI, "strength", MStemp);

		Evec = [Evec; stdE, stdEP];

		BPM_RES = placet_element_get_attribute("electron", BI(1), "resolution")
		BPM_RES = placet_element_get_attribute("positron", BI(1), "resolution")

		disp(' ')
		disp('Element misalignments')
		disp('Quadrupoles')
		xc = placet_element_get_attribute("electron", QI, "x");
		yc = placet_element_get_attribute("electron", QI, "y");
		rc = placet_element_get_attribute("electron", QI, "roll");
		xc2 = placet_element_get_attribute("positron", QI, "x");
		yc2 = placet_element_get_attribute("positron", QI, "y");
		rc2 = placet_element_get_attribute("positron", QI, "roll");
		disp('x_pos - y_pos - roll - x_pos - y_pos - roll')
		[xc, yc, rc, xc2, yc2, rc2]
		stdQ = [sqrt(sum(xc.^2)), sqrt(sum(yc.^2)), sqrt(sum(rc.^2)), sqrt(sum(xc2.^2)), sqrt(sum(yc2.^2)), sqrt(sum(rc2.^2))]
		misQvec = [misQvec; stdQ];

		disp(' ')
		disp('Multipoles')
		xc = placet_element_get_attribute("electron", MI, "x");
		yc = placet_element_get_attribute("electron", MI, "y");
		rc = placet_element_get_attribute("electron", MI, "roll");
		xc2 = placet_element_get_attribute("positron", MI, "x");
		yc2 = placet_element_get_attribute("positron", MI, "y");
		rc2 = placet_element_get_attribute("positron", MI, "roll");
		disp('x_pos - y_pos - roll - x_pos - y_pos - roll')
		[xc, yc, rc, xc2, yc2, rc2]
		stdM = [sqrt(sum(xc.^2)), sqrt(sum(yc.^2)), sqrt(sum(rc.^2)), sqrt(sum(xc2.^2)), sqrt(sum(yc2.^2)), sqrt(sum(rc2.^2))]
		misMvec = [misMvec; stdM];
	
		disp(' ')
		disp('Sbends')
		rc = placet_element_get_attribute("electron", SI, "roll");
		rc2 = placet_element_get_attribute("positron", SI, "roll");
		[rc, rc2]
		stdM = [sqrt(sum(rc.^2)), sqrt(sum(rc2.^2))]
		SbendRollvec = [SbendRollvec; stdM];

		disp(' ')
		disp('Misalignment w r t the beam')
		disp('Quadrupoles')
		[stdx, stdy] = checkAlignment("electron", "beam0t", QI);
		stdQwrt = [stdx, stdy]
		[stdx, stdy] = checkAlignment("positron", "beam0t", QI);
		stdQwrtP = [stdx, stdy]
		wrtbQvec = [wrtbQvec; stdQwrt, stdQwrtP];

		disp('Multipoles')
		[stdx, stdy] = checkAlignment("electron", "beam0t", MI);
		stdMwrt = [stdx, stdy]
		[stdx, stdy] = checkAlignment("positron", "beam0t", MI);
		stdMwrtP = [stdx, stdy]
		wrtbMvec = [wrtbMvec; stdMwrt, stdMwrtP];

	
		disp('Multipole strengths')
		MStemp = placet_element_get_attribute("electron", MI, "strength");
		[MStemp, MS, MStemp-MS]
		stdMS = sqrt(sum((MStemp-MS).^2))
		MStemp = placet_element_get_attribute("positron", MI, "strength");
		[MStemp, MS, MStemp-MS]
		stdMSP = sqrt(sum((MStemp-MS).^2))
		multStreng = [multStreng; stdMS, stdMSP];

	  	disp('SBend angles')
	  	AStemp = placet_element_get_attribute("electron", SI, "angle");
	  	stdAS = sqrt(sum(((AStemp-AS)./AS).^2));
	  	AStemp = placet_element_get_attribute("positron", SI, "angle");
	  	stdASP = sqrt(sum(((AStemp-AS)./AS).^2));
	  	
	  	angleSbends = [angleSbends; stdAS, stdASP];


		Evec
		misQvec
		misMvec
		wrtbQvec
		wrtbMvec
		multStreng
		angleSbends
	
		disp('==========================================================================')
		disp(' ')
	}
}
