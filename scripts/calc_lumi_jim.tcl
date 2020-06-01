# Set up GUINEA PIG
#############################################################################
# GUINEA EXEC PATH SHOULD BE SET GLOBALLY
array set gp_param "
	energy $e0
	particles [expr $match(charge)*1e-10]
	sigmaz $match(sigma_z)
	cut_x 400.0
	cut_y 20.0
	n_x 128
	n_y 256
	n_t 1
	charge_sign -1.0
	ecm_min [expr 2.0*$e0*0.99]"

proc write_guinea_correct {offsetx offsety offdy cuty ny} {
	global n_slice gp_param n_total
	set f [open acc.dat w]
	set offsetnew [expr $offsety + $offdy]

	puts $f "\$ACCELERATOR:: default_clic"
	puts $f "\{energy=$gp_param(energy);particles=$gp_param(particles);"
	puts $f "beta_x=8.0;beta_y=0.1;emitt_x=0.95;"
	puts $f "emitt_y=0.03;sigma_z=$gp_param(sigmaz);"
	puts $f "offset_y.1=$offsety;offset_y.2=$offsety;"
	puts $f "offset_x.1=$offsetx;offset_x.2=$offsetx;\}"

	puts $f "\$ACCELERATOR:: default_clic_old"
	puts $f "\{energy=$gp_param(energy);particles=$gp_param(particles);"
	puts $f "beta_x=7.0;beta_y=0.068;emitt_x=0.66;"
	puts $f "emitt_y=0.02;sigma_z=$gp_param(sigmaz);espread=0.01;dist_z=0;f_rep=100.0;"
	puts $f "offset_y.1=$offsety;offset_y.2=$offsety;"
	puts $f "offset_x.1=$offsetx;offset_x.2=$offsetx;"
	puts $f "n_b=154;waist_y=0;\}"

	puts $f "\$PARAMETERS:: default_simple"
	puts $f "\{n_x=$gp_param(n_x);n_y=$gp_param(n_y);n_z=$n_slice;"
	puts $f "n_t=$gp_param(n_t);n_m=$n_total;cut_x=$gp_param(cut_x);"
	puts $f "cut_y=$gp_param(cut_y);cut_z=3.0*sigma_z.1;"
	puts $f "offset_y.1=$offsety;offset_y.2=$offsety;"
	puts $f "offset_x.1=$offsetx;offset_x.2=$offsetx;"
	puts $f "force_symmetric=0;electron_ratio=1.0;do_photons=1;do_isr=0;do_lumi=0;beam_size=1;"
	puts $f "ecm_min=$gp_param(ecm_min);do_muons=0;"
	puts $f "do_coherent=1;grids=0;rndm_load=1;do_espread=0.0;"
	puts $f "do_pairs=0;track_pairs=0;store_beam=0;do_compt=0;photon_ratio=0.2;load_beam=3;"
	puts $f "do_hadrons=0;store_hadrons=0;do_jets=0;store_jets=0;store_photons=0;"
	puts $f "hist_ee_bins=134;hist_ee_max=2.01*energy.1;charge_sign=$gp_param(charge_sign);"
	puts $f "do_prod=0; prod_e=0.0; prod_scal=0.0; do_cross=0;do_eloss=1;ext_field=0;\}"

	puts $f "\$PARAMETERS:: expanded_simulation"
	puts $f "\{n_x=$gp_param(n_x);n_y=$gp_param(n_y);n_z=$n_slice;"
	puts $f "n_t=$gp_param(n_t);n_m=$n_total;cut_x=$gp_param(cut_x);"
	puts $f "cut_y=$gp_param(cut_y);cut_z=3.0*sigma_z.1;"
	puts $f "offset_y.1=$offsety;offset_y.2=$offsety;"
	puts $f "offset_x.1=$offsetx;offset_x.2=$offsetx;"
	puts $f "force_symmetric=0;electron_ratio=1.0;do_photons=1;store_photons=1;beam_size=1;"
	puts $f "ecm_min=$gp_param(ecm_min);photon_ratio=1.0;do_muons=0;"
	puts $f "do_coherent=0;grids=0;rndm_load=1;do_espread=0;do_eloss=1;"
	puts $f "do_pairs=0;track_pairs=0;store_beam=1;do_compt=0;load_beam=3;"
	puts $f "do_hadrons=0;store_hadrons=0;do_jets=0;store_jets=0;"
	puts $f "hist_ee_bins=1010;hist_ee_max=2.01*energy.1;charge_sign=$gp_param(charge_sign);\}"

	puts $f "\$PARAMETERS:: full_simulation"
	puts $f "\{n_x=$gp_param(n_x);n_y=$gp_param(n_y);n_z=$n_slice;"
	puts $f "n_t=$gp_param(n_t);n_m=$n_total;cut_x=$gp_param(cut_x);"
	puts $f "cut_y=$gp_param(cut_y);cut_z=3.0*sigma_z.1;"
	puts $f "offset_y.1=$offsety;offset_y.2=$offsety;"
	puts $f "offset_x.1=$offsetx;offset_x.2=$offsetx;"
	puts $f "force_symmetric=0;electron_ratio=1.0;do_photons=1;store_photons=1;beam_size=1;"
	puts $f "ecm_min=$gp_param(ecm_min);photon_ratio=1.0;do_muons=0;"
	puts $f "do_coherent=1;grids=7;rndm_load=1;do_espread=0;do_eloss=1;"
	puts $f "do_pairs=1;track_pairs=1;store_beam=1;do_compt=0;load_beam=3;"
	puts $f "do_hadrons=0;store_hadrons=0;do_jets=0;store_jets=0;"
	puts $f "hist_ee_bins=1010;hist_ee_max=2.01*energy.1;charge_sign=$gp_param(charge_sign);\}"

	puts $f "\$PARAMETERS:: full_simulation_old"
	puts $f "\{n_x=$gp_param(n_x);n_y=$gp_param(n_y);n_z=$n_slice;"
	puts $f "n_t=$gp_param(n_t);n_m=$n_total;cut_x=$gp_param(cut_x);"
	puts $f "cut_y=$gp_param(cut_y);cut_z=3.0*sigma_z.1;"
	puts $f "offset_y.1=$offsety;offset_y.2=$offsety;"
	puts $f "offset_x.1=$offsetx;offset_x.2=$offsetx;"
	puts $f "force_symmetric=0;electron_ratio=1.0;do_photons=1;do_isr=1;beam_size=1;"
	puts $f "ecm_min=$gp_param(ecm_min);photon_ratio=1.0;do_muons=0;do_trident=0;"
	puts $f "do_coherent=0;grids=0;rndm_load=0;do_espread=0;"
	puts $f "do_pairs=0;track_pairs=0;store_beam=1;do_compt=0;photon_ratio=0.2;load_beam=3;"
	puts $f "do_hadrons=0;store_hadrons=0;do_jets=0;store_jets=0;store_photons=1;"
	puts $f "hist_ee_bins=1010;hist_ee_max=2.02*energy.1;charge_sign=$gp_param(charge_sign);\}"
	close $f
}

# Return the bunch-crossing luminosity in m^-2
proc get_lumi_ee {name} {
	set l [exec grep lumi_ee= $name]
	set i1 [expr [string last "=" $l]+1]
	set i2 [expr [string last ";" $l]-1]
	return [string range $l $i1 $i2]
}

# Return the bunch-crossing luminosity in m^-2 for the 1% Energy peak
proc get_lumi_high {name} {
	set l [exec grep lumi_ee_high= $name]
	set i1 [expr [string last "=" $l]+1]
	set i2 [expr [string last ";" $l]-1]
	return [string range $l $i1 $i2]
}

# Return the average energy loss per particle [GeV] due to beamstrahlung
proc get_beamstrahlung {name} {
	set l [exec grep de1= $name]
	set i1 [expr [string first "=" $l]+1]
	set i2 [expr [string first ";" $l]-1]
	set e1 [string range $l $i1 $i2]

	set l [exec grep de2= $name]
	set i1 [expr [string last "=" $l]+1]
	set i2 [expr [string last ";" $l]-1]
	set e2 [string range $l $i1 $i2]
	return "$e1 $e2"
}

# Return the angles [urad] of the spent beams
proc get_beam_deflection {name} {
	set l [exec grep bpm_vx.1= $name]
	set i1 [expr [string first "=" $l]+1]
	set i2 [expr [string first ";" $l]-1]
	set vx1 [string range $l $i1 $i2]

	set l [exec grep bpm_vy.1= $name]
	set i1 [expr [string first "=" $l]+1]
	set i2 [expr [string first ";" $l]-1]
	set vy1 [string range $l $i1 $i2]

	set l [exec grep bpm_vx.2= $name]
	set i1 [expr [string first "=" $l]+1]
	set i2 [expr [string first ";" $l]-1]
	set vx2 [string range $l $i1 $i2]

	set l [exec grep bpm_vy.2= $name]
	set i1 [expr [string first "=" $l]+1]
	set i2 [expr [string first ";" $l]-1]
	set vy2 [string range $l $i1 $i2]

	return "$vx1 $vy1 $vx2 $vy2"
}

# return Total luminosity and peak luminosity
proc get_results {name} {
	return "[get_lumi_ee $name] [get_lumi_high $name]"
}

# TCL procedures for running GUINEA-PIG
##################################################################################
# Run default version
proc run_guinea {} {
	global gp_param guinea_exec
	# setup
	write_guinea_correct 0 0 0 0 0

	if { ![ file exist $guinea_exec] } {
		set guinea_exec guinea
		puts "Check Guinea exec path!"
	}

	# Run GUINEA-PIG
	if {[catch {exec $guinea_exec default_clic default_simple result.out}]} {
		puts "GUINEA-PIG failed. Possibly due to wrong grid size. I will try increase vertical grid size"
		set gp_param(n_y) [expr 3*$gp_param(n_y)]
		puts "gp_param(n_y) = $gp_param(n_y) "
		puts "gp_param(cut_y) = $gp_param(cut_y) "
		write_guinea_correct 0 0 0 0 0
		# Try to run with increased vertical grid
		if {[catch {exec $guinea_exec default_clic default_simple result.out}]} {
			puts "GUINEA-PIG failed again. I'm exiting."
			return {-1 -1}
		} else {
			return [get_results result.out]
		}
	} else {
		return [get_results result.out]
	}
}

# run expanded version
proc run_guinea_expanded {} {
	global gp_param guinea_exec
	# setup
	write_guinea_correct 0 0 0 0 0

	if { ![ file exist $guinea_exec] } {
	set guinea_exec guinea
	puts "Check Guinea exec path!"
	}

#	puts "gp_param(n_x) = $gp_param(n_x) "
#	puts "gp_param(cut_x) = $gp_param(cut_x) "
#	puts "gp_param(n_y) = $gp_param(n_y) "
#	puts "gp_param(cut_y) = $gp_param(cut_y) "

	# Run GUINEA-PIG
	if {[catch {exec $guinea_exec default_clic expanded_simulation result.out}]} {
		puts "GUINEA-PIG failed. Possibly due to wrong grid size. I will try increase vertical grid size"
		set gp_param(n_y) [expr 3*$gp_param(n_y)]
		puts "gp_param(n_y) = $gp_param(n_y) "
		puts "gp_param(cut_y) = $gp_param(cut_y) "
		write_guinea_correct 0 0 0 0 0
		# Try to run with increased vertical grid
		if {[catch {exec $guinea_exec default_clic expanded_simulation result.out}]} {
			puts "GUINEA-PIG failed again. I'm exiting."
			return {-1 -1}
		} else {
			return [get_results result.out]
		}
	} else {
		return [get_results result.out]
	}
}

# run full version including tracking of incoherent pairs
proc run_guinea_full {} {
	global gp_param guinea_exec
	# setup
	write_guinea_correct 0 0 0 0 0

	if { ![ file exist $guinea_exec] } {
	set guinea_exec guinea
	puts "Check Guinea exec path!"
	}

	# Run GUINEA-PIG
	if {[catch {exec $guinea_exec default_clic full_simulation result.out}]} {
		puts "GUINEA-PIG failed. Possibly due to wrong grid size. I will try increase vertical grid size"
		set gp_param(n_y) [expr 3*$gp_param(n_y)]
		puts "gp_param(n_y) = $gp_param(n_y) "
		puts "gp_param(cut_y) = $gp_param(cut_y) "
		write_guinea_correct 0 0 0 0 0
		# Try to run with increased vertical grid
		if {[catch {exec $guinea_exec default_clic full_simulation result.out}]} {
			puts "GUINEA-PIG failed again. I'm exiting."
			return {-1 -1}
		} else {
			return [get_results result.out]
		}
	} else {
		return [get_results result.out]
	}
}
# Octave functions for beam distributions
Octave {
	function [Bout] = cut_distribution(B,k)
		# Cut distribution in all dimension at +/- k*sigma_rms
		logvec1 = (B(:,1) > mean(B(:,1))-k*std(B(:,1))) .* (B(:,1) < mean(B(:,1))+k*std(B(:,1)));
		logvec2 = (B(:,2) > mean(B(:,2))-k*std(B(:,2))) .* (B(:,2) < mean(B(:,2))+k*std(B(:,2)));
		logvec3 = (B(:,3) > mean(B(:,3))-k*std(B(:,3))) .* (B(:,3) < mean(B(:,3))+k*std(B(:,3)));
		logvec4 = (B(:,4) > mean(B(:,4))-k*std(B(:,4))) .* (B(:,4) < mean(B(:,4))+k*std(B(:,4)));
		logvec5 = (B(:,5) > mean(B(:,5))-k*std(B(:,5))) .* (B(:,5) < mean(B(:,5))+k*std(B(:,5)));
		logvec6 = (B(:,6) > mean(B(:,6))-k*std(B(:,6))) .* (B(:,6) < mean(B(:,6))+k*std(B(:,6)));
		logvec = logical(logvec1.*logvec2.*logvec3.*logvec4.*logvec5.*logvec6);
		Bout = B(logvec,:);
	end
}

# Return beam parameters, mean and size after cut of distribution
Octave {
	function [X, Y, sx, sy] = get_beam_parameters(B)
		Bcut = cut_distribution(B,3);
		X = mean(Bcut(:,2));
		Y = mean(Bcut(:,3));
		sx = std(Bcut(:,2));
		sy = std(Bcut(:,3));
	end
}

# Octave function for adjusting grid parameters and centralizing the beams
##################################################################################
Octave {
	function [cutx, cuty, nx, ny] = adjust_grid_parameters(no_feedback)
		# the default setting = feedback activated to achieve head-on collision
		if nargin < 1
			no_feedback = 0;
		end

		# Load beams
		B1in = load("electron.ini");
		B2in = load("positron.ini");

		# Mean and rms of distribution after cut:
		[X1, Y1, sx1, sy1] = get_beam_parameters(B1in);
		[X2, Y2, sx2, sy2] = get_beam_parameters(B2in);

		# centralize the two distributions (mean of two beams should align with 0)
		# if no_feedback = 0 centralize each distirbution (perfect head-on collision)
		if no_feedback
			B1in(:,2) = B1in(:,2)-0.5*(X1+X2)*ones(size(B1in(:,2)));
			B2in(:,2) = B2in(:,2)-0.5*(X1+X2)*ones(size(B2in(:,2)));
			B1in(:,3) = B1in(:,3)-0.5*(Y1+Y2)*ones(size(B1in(:,3)));
			B2in(:,3) = B2in(:,3)-0.5*(Y1+Y2)*ones(size(B2in(:,3)));
		else
			B1in(:,2) = B1in(:,2)-X1*ones(size(B1in(:,2)));
			B2in(:,2) = B2in(:,2)-X2*ones(size(B2in(:,2)));
			B1in(:,3) = B1in(:,3)-Y1*ones(size(B1in(:,3)));
			B2in(:,3) = B2in(:,3)-Y2*ones(size(B2in(:,3)));
		end
		[X1, Y1, sx1, sy1] = get_beam_parameters(B1in);
		[X2, Y2, sx2, sy2] = get_beam_parameters(B2in);

		# save beams to file
		save_beam("electron.ini", B1in);
		save_beam("positron.ini", B2in);

		# distance from 0-axis to either distribution [nm]:
		xoff = 1e3*abs(X1-X2)/2;
		yoff = 1e3*abs(Y2-Y1)/2;

		# half-width of the two distributions
		wx = xoff+4*1e3*max([sx1,sx2]);
		wy = yoff+8*1e3*max([sy1,sy2]);

		# scale cell_width by sigma from distribution
		cell_size_x = 1e3*min([sx1,sx2])/10;
		cell_size_y = 1e3*min([sy1,sy2])/10;

		# find sufficient number of cells needed
		nx = 2; ny = 2;
		cutx = cell_size_x*nx/2;
		cuty = cell_size_y*ny/2;
		while cutx < wx
			nx = 2*nx;
			cutx = cell_size_x*nx/2;
		end
		while cuty < wy
			ny = 2*ny;
			cuty = cell_size_y*ny/2;
		end

		if nx > 512
			nx = 512;
			cutx = cell_size_x*nx/2;
		end
		if ny > 1024
			ny = 1024;
			cuty = cell_size_y*ny/2;
		end

		# Update parameters:
		Tcl_Eval(["set gp_param(cut_x) ",num2str(cutx)])
		Tcl_Eval(["set gp_param(cut_y) ",num2str(cuty)])
		Tcl_Eval(["set gp_param(n_x) ",num2str(nx)])
		Tcl_Eval(["set gp_param(n_y) ",num2str(ny)])
	end
}

# Octave function for retrieving photon parameters
##################################################################################
Octave {
	function	[sxp1,sxp2, syp1,syp2] = get_photon_parameters()
		B = load("photon1.dat");
		logvec2 = (B(:,2) > mean(B(:,2))-3*std(B(:,2))) .* (B(:,2) < mean(B(:,2))+3*std(B(:,2)));
		logvec3 = (B(:,3) > mean(B(:,3))-3*std(B(:,3))) .* (B(:,3) < mean(B(:,3))+3*std(B(:,3)));
		logvec = logical(logvec2.*logvec3);
		Bp1 = B(logvec,:);

		B = load("photon2.dat");
		logvec2 = (B(:,2) > mean(B(:,2))-3*std(B(:,2))) .* (B(:,2) < mean(B(:,2))+3*std(B(:,2)));
		logvec3 = (B(:,3) > mean(B(:,3))-3*std(B(:,3))) .* (B(:,3) < mean(B(:,3))+3*std(B(:,3)));
		logvec = logical(logvec2.*logvec3);
		Bp2 = B(logvec,:);

		sxp1 = 1e6*std(Bp1(:,2));
		sxp2 = 1e6*std(Bp2(:,2));
		syp1 = 1e6*std(Bp1(:,3));
		syp2 = 1e6*std(Bp2(:,3));
	end
}

# Octave function for retrieving pairs in BeamCal
##################################################################################
Octave {
	function Pout = pairs_in_BeamCal(P)
		# BeamCal parameters
		L = 3.181; 						# m, distance to start of BeamCal
		r_inner = 32e-3; 				# m
		r_outer = 150e-3; 			# m
		a_min = atan(r_inner/L); 	# rad
		a_max = atan(r_outer/L); 	# rad

		m = 0.511e-3; 					# electron mass in GeV;

		E = P(:,1);
		vx = P(:,2);
		vy = P(:,3);
		vz = P(:,4);

		# compute transverse momentum and polar angle
		P_T = sqrt(vx.^2+vy.^2).*sqrt(E.^2 - m^2);
		angle = atan(sqrt(vx.^2+vy.^2)./abs(vz));
		Pmin = 0.5*4*0.3*L*tan(angle);

		# Cut particles that hit the BeamCal
		ind_hit = logical((angle>=a_min).*(angle<=a_max).*(P_T>=Pmin));
		Pout = P(ind_hit,:);
	end

	function	[Etot1, Etot2, Ne1, Np1, Ne2, Np2] = get_pairs()
		# Load particles
		P = load('pairs.dat');

		# sort positive and negative vz:
		P1 = P(P(:,4)>0,:); 			# In BeamCal 1
		P2 = P(P(:,4)<0,:); 			# In BeamCal 2

		# Pairs that hit BeamCal 1:
		PBC1 = pairs_in_BeamCal(P1);
		E1 = PBC1(:,1);
		Ne1 = sum(E1>0); 					# number of electrons
		Np1 = sum(E1<0); 					# number of positrons
		Etot1 = sum(abs(E1)); 			# total energy

		# Pairs that hit BeamCal 1:
		PBC2 = pairs_in_BeamCal(P2);
		E2 = PBC2(:,1);
		Ne2 = sum(E2>0); 					# number of electrons
		Np2 = sum(E2<0); 					# number of positrons
		Etot2 = sum(abs(E2)); 			# total energy
	end
}


# Octave functions for running GP
##################################################################################
Octave {
	knobfile = fopen("Luminosities.dat","w"); fclose(knobfile);

	function [L, Lpeak] = get_lumi_simple(B1, B2)
		if nargin<2
			IP = placet_get_name_number_list("electron", "IP");
			[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IP);
			[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IP);
		end
		save_beam("electron.ini", B1);
		save_beam("positron.ini", B2);

		# center beams and adjust grid:
		adjust_grid_parameters();

		# run GUINEA-PIG
		Tcl_Eval("set res [run_guinea]");
		Tcl_Eval("set lumi_total [lindex \$res 0]");
		Tcl_Eval("set lumi_peak [lindex \$res 1]");
		L = str2num(Tcl_GetVar("lumi_total"));
		Lpeak = str2num(Tcl_GetVar("lumi_peak"));
	end
}

Octave {
	function [L, Lpeak, beam_size, beam_strahl, beam_defl] = get_lumi(B1, B2, no_feedback)
		if nargin<2
			IP = placet_get_name_number_list("electron", "IP");
			[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IP);
			[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IP);
			no_feedback = 0; # default = feedback on
		end
		if nargin == 2
			no_feedback = 0; # default = feedback on
		end
		save_beam("electron.ini", B1);
		save_beam("positron.ini", B2);

		# center beams and adjust grid:
		adjust_grid_parameters(no_feedback);

		# Run GUINEA-PIG
		#Tcl_Eval("set res [run_guinea_expanded 0.0 0.0]");
		Tcl_Eval("set res [run_guinea_expanded]");
		Tcl_Eval("set lumi_total [lindex \$res 0]");
		Tcl_Eval("set lumi_peak [lindex \$res 1]");
		L = str2num(Tcl_GetVar("lumi_total"));
		Lpeak = str2num(Tcl_GetVar("lumi_peak"));

		if L > 0
			# extract beamstrahlung parameters
			Tcl_Eval("set beamstrahl [get_beamstrahlung result.out]");
			Tcl_Eval("set e1 [lindex \$beamstrahl 0]");
			Tcl_Eval("set e2 [lindex \$beamstrahl 1]");
			e1 = str2num(Tcl_GetVar("e1"));
			e2 = str2num(Tcl_GetVar("e2"));
			# get photon beam divergences
			#[sxp1,sxp2, syp1,syp2] = get_photon_parameters();
			sxp1 = 0; sxp2 = 0; syp1 = 0; syp2 = 0;

			# extract beam deflection angles
			Tcl_Eval("set beamdefl [get_beam_deflection result.out]");
			Tcl_Eval("set vx1 [lindex \$beamdefl 0]");
			Tcl_Eval("set vy1 [lindex \$beamdefl 1]");
			Tcl_Eval("set vx2 [lindex \$beamdefl 2]");
			Tcl_Eval("set vy2 [lindex \$beamdefl 3]");
			vx1 = str2num(Tcl_GetVar("vx1"));
			vy1 = str2num(Tcl_GetVar("vy1"));
			vx2 = str2num(Tcl_GetVar("vx2"));
			vy2 = str2num(Tcl_GetVar("vy2"));
		else
			e1 = 0;
			e2 = 0;
			sxp1 = 0;
			syp1 = 0;
			sxp2 = 0;
			syp2 = 0;
			vx1 = 0;
			vy1 = 0;
			vx2 = 0;
			vy2 = 0
		end

		# compute beam sizes
		B1cut = cut_distribution(B1,3);
		B2cut = cut_distribution(B2,3);
		sx = [std(B1cut(:,2)), std(B2cut(:,2))];
		sy = [std(B1cut(:,3)), std(B2cut(:,3))];

		# define output variables
		beam_size = [sx(1), sy(1), sx(2), sy(2)];
		beam_strahl = [e1, e2, sxp1, syp1, sxp2, syp2];
		beam_defl = [vx1, vy1, vx2, vy2];

		# write to file:
      knobfile = fopen("Luminosities.dat","a");
		fprintf(knobfile,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", L, Lpeak, sx(1), sy(1), sx(2), sy(2), e1, e2, sxp1, syp1, sxp2, syp2, vx1, vy1, vx2, vy2,-1,-1,-1,-1,-1,-1);
		fclose(knobfile);
	end
}



# FULL simulation
Octave {
	function [L, Lpeak, beam_size, beam_strahl, beam_defl, pairs] = get_lumi_full(B1, B2, no_feedback)
		if nargin<2
			IP = placet_get_name_number_list("electron", "IP");
			[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IP);
			[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IP);
			no_feedback = 0; # default = feedback on
		end
		if nargin == 2
			no_feedback = 0; # default = feedback on
		end
		save_beam("electron.ini", B1);
		save_beam("positron.ini", B2);

		# center beams and adjust grid:
		adjust_grid_parameters(no_feedback);

		# Run GUINEA-PIG
		Tcl_Eval("set res [run_guinea_full]");
		Tcl_Eval("set lumi_total [lindex \$res 0]");
		Tcl_Eval("set lumi_peak [lindex \$res 1]");
		L = str2num(Tcl_GetVar("lumi_total"));
		Lpeak = str2num(Tcl_GetVar("lumi_peak"));
		if L > 0
			# extract beamstrahlung parameters
			Tcl_Eval("set beamstrahl [get_beamstrahlung result.out]");
			Tcl_Eval("set e1 [lindex \$beamstrahl 0]");
			Tcl_Eval("set e2 [lindex \$beamstrahl 1]");
			e1 = str2num(Tcl_GetVar("e1"));
			e2 = str2num(Tcl_GetVar("e2"));
			# get photon beam divergences
			#[sxp1,sxp2, syp1,syp2] = get_photon_parameters();
			sxp1 = 0; sxp2 = 0; syp1 = 0; syp2 = 0;

			# extract beam deflection angles
			Tcl_Eval("set beamdefl [get_beam_deflection result.out]");
			Tcl_Eval("set vx1 [lindex \$beamdefl 0]");
			Tcl_Eval("set vy1 [lindex \$beamdefl 1]");
			Tcl_Eval("set vx2 [lindex \$beamdefl 2]");
			Tcl_Eval("set vy2 [lindex \$beamdefl 3]");
			vx1 = str2num(Tcl_GetVar("vx1"));
			vy1 = str2num(Tcl_GetVar("vy1"));
			vx2 = str2num(Tcl_GetVar("vx2"));
			vy2 = str2num(Tcl_GetVar("vy2"));

			# get pairs in BeamCal
			[E1, E2, Ne1, Np1, Ne2, Np2] = get_pairs();
		else
			e1 = 0;
			e2 = 0;
			sxp1 = 0;
			syp1 = 0;
			sxp2 = 0;
			syp2 = 0;
			vx1 = 0;
			vy1 = 0;
			vx2 = 0;
			vy2 = 0;
			E1 = 0;
			E2 = 0;
			Ne1 = 0;
			Np1 = 0;
			Ne2 = 0;
			Np2 = 0;
		end

		# compute beam sizes
		B1cut = cut_distribution(B1,3);
		B2cut = cut_distribution(B2,3);
		sx = [std(B1cut(:,2)), std(B2cut(:,2))];
		sy = [std(B1cut(:,3)), std(B2cut(:,3))];

		# define output variables
		beam_size = [sx(1), sy(1), sx(2), sy(2)];
		beam_strahl = [e1, e2, sxp1, syp1, sxp2, syp2];
		beam_defl = [vx1, vy1, vx2, vy2];
		pairs =	[E1, E2, Ne1, Np1, Ne2, Np2];

		# write to file:
      knobfile = fopen("Luminosities.dat","a");
		fprintf(knobfile,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", L, Lpeak, sx(1), sy(1), sx(2), sy(2), e1, e2, sxp1, syp1, sxp2, syp2, vx1, vy1, vx2, vy2, E1, E2, Ne1, Np1, Ne2, Np2);
		fclose(knobfile);
	end
}
