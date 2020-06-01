# Octave function definitions
#############################################################
	Octave {
	function RM = get_OrbitResponseMatrixMoveQuads(Beamline, Beam0, Quadrupole, BPM, LEV);
		Nbpms = length(BPM);
		Nquads = length(Quadrupole);
		RMxx = zeros(Nbpms, Nquads);
		RMxy = zeros(Nbpms, Nquads);
		RMyx = zeros(Nbpms, Nquads);
		RMyy = zeros(Nbpms, Nquads);

		for j = 1:Nquads
			x0 = placet_element_get_attribute(Beamline, Quadrupole(j), "x");
			placet_element_set_attribute(Beamline, Quadrupole(j), "x", x0 + LEV);
			placet_test_no_correction(Beamline, Beam0, "None");
			orb_temp = placet_get_bpm_readings(Beamline, BPM);
			RMxx(:,j) = orb_temp(:,1);
			RMyx(:,j) = orb_temp(:,2);

			placet_element_set_attribute(Beamline, Quadrupole(j), "x", x0 - LEV);
			placet_test_no_correction(Beamline, Beam0, "None");
			orb_temp = placet_get_bpm_readings(Beamline, BPM);
			RMxx(:,j) -= orb_temp(:,1);
			RMyx(:,j) -= orb_temp(:,2);
			placet_element_set_attribute(Beamline, Quadrupole(j), "x", x0);

			# Check causality! Only rows with bpms AFTER Quadrupole j should have nonzero elements
			Binds = 1:length(BPM);
			rows = Binds(BPM<Quadrupole(j));
			if(~isempty(rows))
				RMxx(rows,j) = zeros(size(RMxx(rows,j)));
				RMyx(rows,j) = zeros(size(RMyx(rows,j)));
			end
		end
		RMxx /= 2*LEV;
		RMyx /= 2*LEV;

		for j = 1:Nquads
			y0 = placet_element_get_attribute(Beamline, Quadrupole(j), "y");
			placet_element_set_attribute(Beamline, Quadrupole(j), "y", y0 + LEV);
			placet_test_no_correction(Beamline, Beam0, "None");
			orb_temp = placet_get_bpm_readings(Beamline, BPM);
			RMxy(:,j) = orb_temp(:,1);
			RMyy(:,j) = orb_temp(:,2);

			placet_element_set_attribute(Beamline, Quadrupole(j), "y", y0 - LEV);
			placet_test_no_correction(Beamline, Beam0, "None");
			orb_temp = placet_get_bpm_readings(Beamline, BPM);
			RMxy(:,j) -= orb_temp(:,1);
			RMyy(:,j) -= orb_temp(:,2);
			placet_element_set_attribute(Beamline, Quadrupole(j), "y", y0);

			# Check causality! Only rows with bpms AFTER Quadrupole j should have nonzero elements
			Binds = 1:length(BPM);
			rows = Binds(BPM<Quadrupole(j));
			if(~isempty(rows))
				RMxy(rows,j) = zeros(size(RMxy(rows,j)));
				RMyy(rows,j) = zeros(size(RMyy(rows,j)));
			end
		end
		RMxy /= 2*LEV;
		RMyy /= 2*LEV;

		RM = [RMxx, RMxy; RMyx, RMyy];
		#RM = [RMxx, zeros(size(RMxy)); zeros(size(RMyx)), RMyy];
	end
}

Octave {
	function RM = get_OrbitResponseMatrixRollQuads(Beamline, Beam0, Quadrupole, BPM, LEV);
		Nbpms = length(BPM);
		Nquads = length(Quadrupole);
		RMxr = zeros(Nbpms, Nquads);
		RMyr = zeros(Nbpms, Nquads);
		for j = 1:Nquads
			r0 = placet_element_get_attribute(Beamline, Quadrupole(j), "roll");
			placet_element_set_attribute(Beamline, Quadrupole(j), "roll", r0 + LEV);
			orb_temp = placet_get_bpm_readings(Beamline, BPM);
			RMxr(:,j) = orb_temp(:,1);
			RMyr(:,j) = orb_temp(:,2);

			placet_element_set_attribute(Beamline, Quadrupole(j), "roll", r0 - LEV);
			orb_temp = placet_get_bpm_readings(Beamline, BPM);
			RMxr(:,j) -= orb_temp(:,1);
			RMyr(:,j) -= orb_temp(:,2);
			placet_element_set_attribute(Beamline, Quadrupole(j), "roll", r0);

			# Check causality: Only rows with bpms AFTER Quadrupole j should have nonzero elements
			Binds = 1:length(BPM);
			rows = Binds(BPM<Quadrupole(j));
			if(~isempty(rows))
				RMxr(rows,j) = zeros(size(RMxr(rows,j)));
				RMyr(rows,j) = zeros(size(RMyr(rows,j)));
			end
		end
		RMxr /= 2*LEV;
		RMyr /= 2*LEV;
		RM = [RMxr; RMyr];
	end
}

Octave {
	function RM = get_DispersionResponseMatrixMoveQuads(Beamline, Beam1, Beam2, Quadrupole, BPM, LEV);
		Nbpms = length(BPM);
		Nquads = length(Quadrupole);
		RMxx = zeros(Nbpms, Nquads);
		RMxy = zeros(Nbpms, Nquads);
		RMyx = zeros(Nbpms, Nquads);
		RMyy = zeros(Nbpms, Nquads);
		for j = 1:Nquads
			x0 = placet_element_get_attribute(Beamline, Quadrupole(j), "x");
			placet_element_set_attribute(Beamline, Quadrupole(j), "x", x0 + LEV);
			eta_temp = MeasureDispersion(Beamline, Beam1, Beam2, $deltae, BPM);
			RMxx(:,j) = eta_temp(:,1);
			RMyx(:,j) = eta_temp(:,2);

			placet_element_set_attribute(Beamline, Quadrupole(j), "x", x0 - LEV);
			eta_temp = MeasureDispersion(Beamline, Beam1, Beam2, $deltae, BPM);
			RMxx(:,j) -= eta_temp(:,1);
			RMyx(:,j) -= eta_temp(:,2);
			placet_element_set_attribute(Beamline, Quadrupole(j), "x", x0);

			# Check causality! Only rows with bpms AFTER Quadrupole j should have nonzero elements
			Binds = 1:length(BPM);
			rows = Binds(BPM<Quadrupole(j));
			if(~isempty(rows))
				RMxx(rows,j) = zeros(size(RMxx(rows,j)));
				RMyx(rows,j) = zeros(size(RMyx(rows,j)));
			end
		end
		RMxx /= 2*LEV;
		RMyx /= 2*LEV;

		for j = 1:Nquads
			y0 = placet_element_get_attribute(Beamline, Quadrupole(j), "y");
			placet_element_set_attribute(Beamline, Quadrupole(j), "y", y0 + LEV);
			eta_temp = MeasureDispersion(Beamline, Beam1, Beam2, $deltae, BPM);
			RMxy(:,j) = eta_temp(:,1);
			RMyy(:,j) = eta_temp(:,2);

			placet_element_set_attribute(Beamline, Quadrupole(j), "y", y0 - LEV);
			eta_temp = MeasureDispersion(Beamline, Beam1, Beam2, $deltae, BPM);
			RMxy(:,j) -= eta_temp(:,1);
			RMyy(:,j) -= eta_temp(:,2);
			placet_element_set_attribute(Beamline, Quadrupole(j), "y", y0);

			# Check causality! Only rows with bpms AFTER Quadrupole j should have nonzero elements
			Binds = 1:length(BPM);
			rows = Binds(BPM<Quadrupole(j));
			if(~isempty(rows))
				RMxy(rows,j) = zeros(size(RMxy(rows,j)));
				RMyy(rows,j) = zeros(size(RMyy(rows,j)));
			end
		end
		RMxy /= 2*LEV;
		RMyy /= 2*LEV;

		RM = [RMxx, RMxy; RMyx, RMyy];
		#RM = [RMxx, zeros(size(RMxy)); zeros(size(RMyx)), RMyy];
	end
}

Octave {
	function RM = get_DispersionResponseMatrixRollQuads(Beamline, Beam1, Beam2, Quadrupole, BPM, LEV);
		Nbpms = length(BPM);
		Nquads = length(Quadrupole);
		RMxr = zeros(Nbpms, Nquads);
		RMyr = zeros(Nbpms, Nquads);
		for j = 1:Nquads
			r0 = placet_element_get_attribute(Beamline, Quadrupole(j), "roll");
			placet_element_set_attribute(Beamline, Quadrupole(j), "roll", r0 + LEV);
			eta_temp = MeasureDispersion(Beamline, Beam1, Beam2, $deltae, BPM);
			RMxr(:,j) = eta_temp(:,1);
			RMyr(:,j) = eta_temp(:,2);

			placet_element_set_attribute(Beamline, Quadrupole(j), "roll", r0 - LEV);
			eta_temp = MeasureDispersion(Beamline, Beam1, Beam2, $deltae, BPM);
			RMxr(:,j) -= eta_temp(:,1);
			RMyr(:,j) -= eta_temp(:,2);
			placet_element_set_attribute(Beamline, Quadrupole(j), "roll", r0);

			# Check causality: Only rows with bpms AFTER Quadrupole j should have nonzero elements
			Binds = 1:length(BPM);
			rows = Binds(BPM<Quadrupole(j));
			if(~isempty(rows))
				RMxr(rows,j) = zeros(size(RMxr(rows,j)));
				RMyr(rows,j) = zeros(size(RMyr(rows,j)));
			end
		end
		RMxr /= 2*LEV;
		RMyr /= 2*LEV;
		RM = [RMxr; RMyr];
	end
}


# Procedure for trajectory and dispersion target steering
#############################################################
proc run_BBA_single_beam { beamline } {
	puts "\nPerforming beam-based alignment on $beamline beamline"
	global wdisp iterdts gaindts deltae
	# Measure response matrices on actual machine:
	Octave {
		disp('Measure the BBA response matrix');
		MStemp = placet_element_get_attribute("$beamline", MI, "strength");
		placet_element_set_attribute("$beamline", MI, "strength", complex(0.0,0.0));

		dp = 1;
		RMorb = get_OrbitResponseMatrixMoveQuads("$beamline", "beam0", QI, BI, dp);
		RMdisp = get_DispersionResponseMatrixMoveQuads("$beamline", "beam1", "beam2", QI, BI, dp);
		placet_element_set_attribute("$beamline", MI, "strength", complex(MStemp,0.0));

		disp(" ")
		disp("START OF DTS - MOVING QUADS")
		disp(" ")
		disp("Turn off multipoles")
		MStemp = placet_element_get_attribute("$beamline", MI, "strength");
		placet_element_set_attribute("$beamline", MI, "strength", complex(0.0,0.0));

		# target orbit
		orb1 = zeros(size(eta0));

		# target excitations
		k1 = zeros(2*length(QI), 1);

		w_orb = 1;
		w_disp = $wdisp
		iter_dts = $iterdts;
		gain = $gaindts;
		limi = 1e4;

		# Measure orbit and dispersion
		placet_test_no_correction("$beamline", "beam0", "None");
		orb_meas_temp = placet_get_bpm_readings("$beamline", BI);
		orb_meas = [orb_meas_temp(:,1); orb_meas_temp(:,2)];
		eta_meas_temp = MeasureDispersion("$beamline", "beam1", "beam2", $deltae, BI);
		eta_meas = [eta_meas_temp(:,1); eta_meas_temp(:,2)];

		y_orb = orb_meas - orb1;
		y_disp = eta_meas - eta0;
		y1 = [w_orb*y_orb; w_disp*y_disp];
		chi2 = [sqrt(sum(y1.^2)), sqrt(sum(y_orb.^2)), sqrt(sum(y_disp.^2))]
		for i=1:iter_dts
			RR1 = [w_orb*RMorb; w_disp*RMdisp];
			[u, s, v] = svd(RR1);
			maxs = max(max(s));
			sinv = (zeros(size(s)))';
			count = 0;
			for k = 1:min(size(s))
				if (s(k,k) < maxs/limi)
					sinv(k,k) = 0;
					count += 1;
				else
					sinv(k,k) = 1/s(k,k);
				end
			end
			
			delta_pos = -gain*(v*sinv*u')*y1;
			max(abs(delta_pos))

			# set new values
			delta_Qx = delta_pos(1:length(QI));
			delta_Qy = delta_pos(length(QI)+1:2*length(QI));
			[delta_Qx, delta_Qy];

			placet_element_vary_attribute("$beamline", QI, "x", delta_Qx);
			placet_element_vary_attribute("$beamline", QI, "y", delta_Qy);

			# Measure orbit and dispersion
			placet_test_no_correction("$beamline", "beam0", "None");
			orb_meas_temp = placet_get_bpm_readings("$beamline", BI);
			orb_meas = [orb_meas_temp(:,1); orb_meas_temp(:,2)];
			eta_meas_temp = MeasureDispersion("$beamline", "beam1", "beam2", $deltae, BI);
			eta_meas = [eta_meas_temp(:,1); eta_meas_temp(:,2)];

			y_orb = orb_meas - orb1;
			y_disp = eta_meas - eta0;
			y1 = [w_orb*y_orb; w_disp*y_disp];
			chi2 = [sqrt(sum(y1.^2)), sqrt(sum(y_orb.^2)), sqrt(sum(y_disp.^2))]
		end
		
		disp("End of DTS")
		disp("Turn on multipoles")
		placet_element_set_attribute("$beamline", MI, "strength", complex(MStemp,0.0));
	}
}



# Procedure for trajectory and dispersion target steering using 2e4 particles
#############################################################
proc run_BBA_single_beam_2e4 { beamline } {
	puts "\nPerforming beam-based alignment on $beamline beamline using 2e4 particles"
	global wdisp iterdts gaindts deltae
	# Measure response matrices on actual machine:
	Octave {
		disp('Measure the BBA response matrix');
		MStemp = placet_element_get_attribute("$beamline", MI, "strength");
		placet_element_set_attribute("$beamline", MI, "strength", complex(0.0,0.0));

		dp = 1;
		RMorb = get_OrbitResponseMatrixMoveQuads("$beamline", "beam0t", QI, BI, dp);
		RMdisp = get_DispersionResponseMatrixMoveQuads("$beamline", "beam1t", "beam2t", QI, BI, dp);
		placet_element_set_attribute("$beamline", MI, "strength", complex(MStemp,0.0));

		disp(" ")
		disp("START OF DTS - MOVING QUADS")
		disp(" ")
		disp("Turn off multipoles")
		MStemp = placet_element_get_attribute("$beamline", MI, "strength");
		placet_element_set_attribute("$beamline", MI, "strength", complex(0.0,0.0));

		# target orbit
		orb1 = zeros(size(eta0));

		# target excitations
		k1 = zeros(2*length(QI), 1);

		w_orb = 1;
		w_disp = $wdisp
		iter_dts = $iterdts;
		gain = $gaindts;
		limi = 1e4;

		# Measure orbit and dispersion
		placet_test_no_correction("$beamline", "beam0t", "None");
		orb_meas_temp = placet_get_bpm_readings("$beamline", BI);
		orb_meas = [orb_meas_temp(:,1); orb_meas_temp(:,2)];
		eta_meas_temp = MeasureDispersion("$beamline", "beam1t", "beam2t", $deltae, BI);
		eta_meas = [eta_meas_temp(:,1); eta_meas_temp(:,2)];

		y_orb = orb_meas - orb1;
		y_disp = eta_meas - eta0;
		y1 = [w_orb*y_orb; w_disp*y_disp];
		chi2 = [sqrt(sum(y1.^2)), sqrt(sum(y_orb.^2)), sqrt(sum(y_disp.^2))]
		for i=1:iter_dts
			RR1 = [w_orb*RMorb; w_disp*RMdisp];
			[u, s, v] = svd(RR1);
			maxs = max(max(s));
			sinv = (zeros(size(s)))';
			count = 0;
			for k = 1:min(size(s))
				if (s(k,k) < maxs/limi)
					sinv(k,k) = 0;
					count += 1;
				else
					sinv(k,k) = 1/s(k,k);
				end
			end
			
			delta_pos = -gain*(v*sinv*u')*y1;
			max(abs(delta_pos))

			# set new values
			delta_Qx = delta_pos(1:length(QI));
			delta_Qy = delta_pos(length(QI)+1:2*length(QI));
			[delta_Qx, delta_Qy];

			placet_element_vary_attribute("$beamline", QI, "x", delta_Qx);
			placet_element_vary_attribute("$beamline", QI, "y", delta_Qy);

			# Measure orbit and dispersion
			placet_test_no_correction("$beamline", "beam0t", "None");
			orb_meas_temp = placet_get_bpm_readings("$beamline", BI);
			orb_meas = [orb_meas_temp(:,1); orb_meas_temp(:,2)];
			eta_meas_temp = MeasureDispersion("$beamline", "beam1t", "beam2t", $deltae, BI);
			eta_meas = [eta_meas_temp(:,1); eta_meas_temp(:,2)];

			y_orb = orb_meas - orb1;
			y_disp = eta_meas - eta0;
			y1 = [w_orb*y_orb; w_disp*y_disp];
			chi2 = [sqrt(sum(y1.^2)), sqrt(sum(y_orb.^2)), sqrt(sum(y_disp.^2))]
		end
		
		disp("End of DTS")
		disp("Turn on multipoles")
		placet_element_set_attribute("$beamline", MI, "strength", complex(MStemp,0.0));
	}
}


# Procedure for 1-2-1 steering
#############################################################
proc run_1_2_1_single_beam { beamline } {
	puts "\nPerforming 1-2-1 steering on the $beamline beamline"
	Octave {
		placet_test_no_correction("$beamline", "beam0", "None");
		C1 = placet_get_bpm_readings("$beamline");
		T = bpm_target_$beamline;
		diff_start = C1-T
		Bstart = placet_get_beam("beam0");
		[E,B] = placet_test_no_correction("$beamline", "beam0", "None",1,0,BI(1));
		shots_count = 0;
		
		for i = 1:length(QI)
		   j = i+1;
		   for k = 1:2
		      placet_test_no_correction("$beamline", "beam0", "None",1,0,BI(j));
		      shots_count = shots_count + 1;
		      C1 = placet_get_bpm_readings("$beamline");
		      px = placet_element_get_attribute("$beamline", QI(i), "x");
		      py = placet_element_get_attribute("$beamline", QI(i), "y");

		      placet_element_set_attribute("$beamline", QI(i), "x",px+1);
		      placet_test_no_correction("$beamline", "beam0", "None",1,0,BI(j));
		      shots_count = shots_count + 1;		      
		      C2 = placet_get_bpm_readings("$beamline");
		      R = [C2(j,1)-C1(j,1), C2(j,2)-C1(j,2)];

		      placet_element_set_attribute("$beamline", QI(i), "x",px);
		      placet_element_set_attribute("$beamline", QI(i), "y",py+1);
		      placet_test_no_correction("$beamline", "beam0", "None",1,0,BI(j));
		      shots_count = shots_count + 1;		      
		      C2 = placet_get_bpm_readings("$beamline");
		      R = [R; C2(j,1)-C1(j,1), C2(j,2)-C1(j,2); 1*eye(2)];
		      placet_element_set_attribute("$beamline", QI(i), "y",py);

		      step = pinv(R)*[T(j,:)-C1(j,:) zeros(1,2)]';
		      placet_element_vary_attribute("$beamline", QI(i), "x", step(1));
		      placet_element_vary_attribute("$beamline", QI(i), "y", step(2));      
		   end   
		end

		placet_test_no_correction("$beamline", "beam0", "None");
		C = placet_get_bpm_readings("$beamline");
		diff_start
		diff_end = C-T
		shots_count
	}
}

