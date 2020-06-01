################################################################
# Octave function definitions
################################################################
# measure dispersion (imperfect BPMs)
Octave {
	function Eta = MeasureDispersion(Beamline, Beam1, Beam2, delta, BI)
		placet_test_no_correction(Beamline, Beam2, "None"); Eta  = placet_get_bpm_readings(Beamline, BI);
		placet_test_no_correction(Beamline, Beam1, "None"); Eta -= placet_get_bpm_readings(Beamline, BI);
		Eta /= 2 * delta;
	end
}

# check dispersion (perfect BPMs)
Octave {
	function Eta = CheckDispersion(Beamline, Beam1, Beam2, delta, BI)
		placet_element_set_attribute(Beamline, BI, "resolution", 0.0);
		BIx = placet_element_get_attribute(Beamline, BI, "x");
		BIy = placet_element_get_attribute(Beamline, BI, "y");
		BIr = placet_element_get_attribute(Beamline, BI, "roll");
		placet_element_set_attribute(Beamline, BI, "x", 0.0);
		placet_element_set_attribute(Beamline, BI, "y", 0.0);
		placet_element_set_attribute(Beamline, BI, "roll", 0.0);

		try
			placet_test_no_correction(Beamline, Beam2, "None"); Eta  = placet_get_bpm_readings(Beamline, BI);
			placet_test_no_correction(Beamline, Beam1, "None"); Eta -= placet_get_bpm_readings(Beamline, BI);
			Eta /= 2 * delta;
		catch
			disp("Beam is probably lost");
			Eta = 1e4*ones(length(BI),2);
		end

		placet_element_set_attribute(Beamline, BI, "x", BIx);
		placet_element_set_attribute(Beamline, BI, "y", BIy);
		placet_element_set_attribute(Beamline, BI, "roll", BIr);
		placet_element_set_attribute(Beamline, BI, "resolution", $bpmres);
	end
}

Octave {
	function Eta = CheckDispersion2(Beamline, Beam1, Beam2, delta, Inds)
		Eta = zeros(length(Inds),2);
		for j = 1:length(Inds)
			[E,B] = placet_test_no_correction(Beamline, Beam2, "None", 1, 0, Inds(j));
			Eta(j,1) = mean(B(:,2));
			Eta(j,2) = mean(B(:,3));
			[E,B] = placet_test_no_correction(Beamline, Beam1, "None", 1, 0, Inds(j));
			Eta(j,1) = Eta(j,1) - mean(B(:,2));
			Eta(j,2) = Eta(j,2) - mean(B(:,3));
		end
		Eta /= 2 * delta;
	end
}

Octave {
	function [] = displayElementStatus(beamline,indices)
		fprintf("\nIndex\tName\tx\t\ty\t\troll\t\tstrengt\n");
		switch nargin
		case 1
			indices = 1:placet_get_name_number_list(beamline, "IP");
		end

		for j = 1:length(indices)
			fprintf("%i",j)
			fprintf("\t%s", placet_element_get_attribute(beamline, indices(j), "name"))
			fprintf("\t%e", placet_element_get_attribute(beamline, indices(j), "x"))
			fprintf("\t%e", placet_element_get_attribute(beamline, indices(j), "y"))
			fprintf("\t%e", placet_element_get_attribute(beamline, indices(j), "roll"))
			typ = placet_element_get_attribute(beamline, indices(j),"type_name");
			switch typ
			case "quadrupole"
				fprintf("\t%e", placet_element_get_attribute(beamline, indices(j), "strength"));
			case "multipole"
				fprintf("\t%e", placet_element_get_attribute(beamline, indices(j), "strength"));
			end
			fprintf("\n");
		end
	end
}

Octave {
	function chi2 = getTotalOffset(beamline,indices,prop)
		chi2 = 0;
		for j = 1:length(indices)
			temp = placet_element_get_attribute(beamline, indices(j), prop);
			chi2 = chi2 + temp^2;
		end
	end
}

Octave {
	function [stdx, stdy] = checkAlignment(beamline, beam, indices)
		xlem = zeros(size(indices));
		ylem = zeros(size(indices));
		xbeam = zeros(size(indices));
		ybeam = zeros(size(indices));
		xdiff = zeros(size(indices));
		ydiff = zeros(size(indices));
		try
			for j = 1:length(indices)
				xelem(j) = placet_element_get_attribute(beamline, indices(j), "x");
				yelem(j) = placet_element_get_attribute(beamline, indices(j), "y");
				[E,B] = placet_test_no_correction(beamline, beam, "None", 1, 0, indices(j));
				xbeam(j) = mean(B(:,2));
				ybeam(j) = mean(B(:,3));
				diffx(j) = xelem(j) - xbeam(j);
				diffy(j) = yelem(j) - ybeam(j);
			end
			disp('diff_x - diff_y - beam_x - beam_y')
			[diffx', diffy', xbeam', ybeam']
			stdx = sqrt(sum(diffx.^2));
			stdy = sqrt(sum(diffy.^2));
		catch
			disp("Beam is probably lost");
			stdx = 1e4;
			stdy = 1e4;
		end
	end
}

Octave {
	function [diffx, diffy] = getAlignment(beamline, beam, indices)
		xlem = zeros(size(indices));
		ylem = zeros(size(indices));
		xbeam = zeros(size(indices));
		ybeam = zeros(size(indices));
		xdiff = zeros(size(indices));
		ydiff = zeros(size(indices));
		try
			for j = 1:length(indices)
				xelem(j) = placet_element_get_attribute(beamline, indices(j), "x");
				yelem(j) = placet_element_get_attribute(beamline, indices(j), "y");
				[E,B] = placet_test_no_correction(beamline, beam, "None", 1, 0, indices(j));
				xbeam(j) = mean(B(:,2));
				ybeam(j) = mean(B(:,3));
				diffx(j) = xelem(j) - xbeam(j);
				diffy(j) = yelem(j) - ybeam(j);
			end
		catch
			disp("Beam is probably lost");
		end
	end
}


################################################################
# TCL procedure to save and load machine status
proc save_beamline_status { beamline_name file_name } {
	Octave {
		disp(" ");
		disp("Saving machine status to file $file_name")
		STATUS = placet_element_get_attributes("$beamline_name");
		save -text $file_name STATUS;
	}
}

proc load_beamline_status { beamline_name path_name } {
	Octave {
		if exist("$path_name", "file")
			clear STATUS;
			disp("reading machine status: $path_name")
			eval(['load ', "$path_name"]);
			placet_element_set_attributes("$beamline_name", STATUS);
			placet_beamline_refresh("$beamline_name")
		else
			disp("ERROR: machine status file not found");
			exit(1);
		end
	}
}

################################################################
# TCL procedure to save tuning data
proc save_tuning_data { file_name t_1 } {
	set t_2 [clock seconds]
	Octave {
		disp(" ");
		disp("Saving simulation data to file $file_name");
		disp(" ");
		disp("Simulation INFO:");
		disp(['Number of Luminosity measurements = ', num2str(Lmeas)]);
		t_sim = ($t_2 - $t_1)/60;
		disp(['Total simulation time = ', num2str(t_sim), ' minutes']);
		save $file_name Lumivec LumiPeakvec Lmeasvec sxvec syvec Evec EwMvec misQvec misMvec wrtbQvec wrtbMvec Lmeas t_sim
	}
}

################################################################
# Octave function to make Gaussian fit to core of beam
Octave {
	function [K,mu,sig] = fit_gauss(y,r)
		switch nargin
		case 1
			r = 1;
		end

		# make histogram
		[NN,XX]=hist(y,100);

		# Gaussian fit
		Y = NN';
		X = XX';

		# Fit only to core of histogram
		Y=Y(r:end+1-r);
		X=X(r:end+1-r);

		A = [X.^2, X, ones(size(X))];
		b = log(Y);
		xfit = A\b;

		# extract parameters
		sig = sqrt(-1/(2*xfit(1)));
		mu = xfit(2)*sig^2;
		K = exp(xfit(3) - xfit(2)^2/(2*xfit(1)));
	end
}

Octave {
	function [sx_core, sy_core, Kx, Ky, mux, muy] = get_core_beam_size(B)
		x = B(:,2);
		y = B(:,3);
		sx = std(x);
		sy = std(y);

		# remove mean and cut at 4 rms
		x = x-mean(x);
		x = x(x<4*sx);
		x = x(x>-4*sx);
		y = y-mean(y);
		y = y(y<4*sy);
		y = y(y>-4*sy);

		# fit only to core
		[Kx,mux,sx_core] = fit_gauss(x, 30);
		[Ky,muy,sy_core] = fit_gauss(y, 40);
	end
}

################################################################
# Cut beam distribition at k times original rms
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

################################################################
# Perform waist scan
Octave {
	function [dx, dy, sx, sy, l] = waist_scan(B)
		l = linspace(-600e-6, 600e-6, 1001);
		sx = zeros(size(l));
		sy = zeros(size(l));

		# Calculate beam matrix
		CM = cov(B);
		R = eye(6);

		# propagate beam matrix
		data = [];
		for i = 1:length(l)
			R(2,5) = l(i); R(3,6) = l(i);
			CM2 = R*CM*transpose(R);
			sx(i) = sqrt(CM2(2,2));
			sy(i) = sqrt(CM2(3,3));
		end

		[mx,ix] = min(sx);
		[my,iy] = min(sy);
		dx = l(ix); # distance to minimum waist_x
		dy = l(iy); # distance to minimum waist_y
	end
}
