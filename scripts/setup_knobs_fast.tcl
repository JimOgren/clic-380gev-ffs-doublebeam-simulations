# Load Knobs-file
##########################################################
Octave {
	load '$script_dir/main/Knobs_sr_${sr}.dat'
	global Knobs_sext Knobs_oct TEMP;
	Knobs_sext = Knobs_sextupole;
	Knobs_oct = Knobs_octupole;
	TEMP = [];

	saveLuminosityData = 1;
}

# Sextupole 1st order knobs
##########################################################
foreach axis { x y } {
	foreach knob { 1 2 3 4 5 6 } {
		Octave {
			function vary_electron_sextupole_knob${knob}${axis}(x)
				global Knobs_sext;
				knobs = x * Knobs_sext.V${axis}(:,$knob);
				placet_element_vary_attribute("electron", Knobs_sext.I${axis}, "$axis", knobs);
			end

			function vary_positron_sextupole_knob${knob}${axis}(x)
				global Knobs_sext;
				knobs = x * Knobs_sext.V${axis}(:,$knob);
				placet_element_vary_attribute("positron", Knobs_sext.I${axis}, "$axis", knobs);
			end
		}
	}
}


# Procedure: Sextupole 1st order knobs fast
##########################################################
proc sextupole_knobs_fast { range machine } {
	puts " "
	puts "Starting tuning with pairs. FAST knobs. GP Parallel exectuion."
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder
   foreach knob { 1 2 3 } {

		set axis x;
		puts "knobs = $knob, axis = $axis"
		puts "Scan electron beamline"

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range,$range,29);
			for jk=1:length(xvec)
				vary_electron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_electron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range,$range,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_electron_sextupole_knob${knob}${axis}(xmax);

			# save scan result:
			if saveLuminosityData
				fname = [$save_dir, "/Res_electron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
				eval(["save ", fname," Res kfit ymax xmax"]);
			end
		}

		##########################################################
		puts "Scan positron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range,$range,29);
			for jk=1:length(xvec)
				vary_positron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_positron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range,$range,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_positron_sextupole_knob${knob}${axis}(xmax);

			# save scan result:
			if saveLuminosityData
				fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
				eval(["save ", fname," Res kfit ymax xmax"]);
			end

		}

		##########################################################
		# Scan vertical direction
		set axis y;
		puts "knobs = $knob, axis = $axis"
		puts "Scan electron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range/5,$range/5,29);
			for jk=1:length(xvec)
				vary_electron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_electron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range/5,$range/5,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end

			vary_electron_sextupole_knob${knob}${axis}(xmax);

			if saveLuminosityData
				fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
				eval(["save ", fname," Res kfit ymax xmax"]);
			end
		}

		##########################################################
		puts "Scan positron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range/5,$range/5,29);
			for jk=1:length(xvec)
				vary_positron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_positron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range/5,$range/5,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_positron_sextupole_knob${knob}${axis}(xmax);

			if saveLuminosityData
				fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
				eval(["save ", fname," Res kfit ymax xmax"]);
			end
		}

		Octave {
			Lmeas = Lmeas + 4*29;
		}
	}
}



# Procedure: Sextupole 1st order knobs fast - 2nd version: knobs 4-6
##########################################################
proc sextupole_knobs_fast_2 { range machine } {
	puts " "
	puts "Starting tuning with pairs. FAST knobs. GP Parallel exectuion."
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder
   foreach knob { 4 5 6 } {

		set axis x;
		puts "knobs = $knob, axis = $axis"
		puts "Scan electron beamline"

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range,$range,29);
			for jk=1:length(xvec)
				vary_electron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_electron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range,$range,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_electron_sextupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_electron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		##########################################################
		puts "Scan positron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range,$range,29);
			for jk=1:length(xvec)
				vary_positron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_positron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range,$range,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_positron_sextupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		##########################################################
		# Scan vertical direction
		set axis y;
		puts "knobs = $knob, axis = $axis"
		puts "Scan electron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range/5,$range/5,29);
			for jk=1:length(xvec)
				vary_electron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_electron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range/5,$range/5,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_electron_sextupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_electron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		##########################################################
		puts "Scan positron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range/5,$range/5,29);
			for jk=1:length(xvec)
				vary_positron_sextupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_positron_sextupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range/5,$range/5,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_positron_sextupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		Octave {
			Lmeas = Lmeas + 4*29;
		}
	}
}





##########################################################
##########################################################
##########################################################

# Octupole 1st order knobs
##########################################################
foreach axis { x y } {
	foreach knob { 1 2 } {
		Octave {
			function vary_electron_octupole_knob${knob}${axis}(x)
				global Knobs_oct;
				knobs = x * Knobs_oct.V${axis}(:,$knob);
				placet_element_vary_attribute("electron", Knobs_oct.I${axis}, "$axis", knobs);
			end

			function vary_positron_octupole_knob${knob}${axis}(x)
				global Knobs_oct;
				knobs = x * Knobs_oct.V${axis}(:,$knob);
				placet_element_vary_attribute("positron", Knobs_oct.I${axis}, "$axis", knobs);
			end
		}
	}
}

# Procedure: Octupole 1st order knobs fast
##########################################################
proc octupole_knobs_fast { range machine } {
	puts " "
	puts "Starting tuning with pairs. Octupole FAST knobs. GP Parallel exectuion."
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder
   foreach knob { 1 2 } {

		set axis x;
		puts "knobs = $knob, axis = $axis"
		puts "Scan electron beamline"

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range,$range,29);
			for jk=1:length(xvec)
				vary_electron_octupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_electron_octupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range,$range,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_electron_octupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_electron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		##########################################################
		puts "Scan positron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range,$range,29);
			for jk=1:length(xvec)
				vary_positron_octupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_positron_octupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range,$range,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_positron_octupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		##########################################################
		# Scan vertical direction
		set axis y;
		puts "knobs = $knob, axis = $axis"
		puts "Scan electron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range/5,$range/5,29);
			for jk=1:length(xvec)
				vary_electron_octupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_electron_octupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range/5,$range/5,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end

			vary_electron_octupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_electron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		##########################################################
		puts "Scan positron beamline"
		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"
		Octave {
			xvec = linspace(-$range/5,$range/5,29);
			for jk=1:length(xvec)
				vary_positron_octupole_knob${knob}${axis}(+xvec(jk));

				# track:
				[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
				[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

				# save beams to temp-folder
				str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
				str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
				eval(str1);
				eval(str2);

				vary_positron_octupole_knob${knob}${axis}(-xvec(jk));
			end
			pause(5);
		}

		# execute GUINEA-PIG in parallel:
		exec bash -c "for i in `seq 1 1 10`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 11 1 20`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"
		exec bash -c "for i in `seq 21 1 29`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

		Octave {
			pause(5);
			xvec = linspace(-$range/5,$range/5,29);
			Res = [];
			for j = 1:length(xvec)
				eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
				Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
			end
			Res;
			Y = Res(:,18) + Res(:,19);
			[Res(:,1:2), Y]

			# make parabolic fit:
			xfit = linspace(xvec(1),xvec(end),200);
			A = [ones(size(xvec')), xvec', (xvec').^2];
			kfit = pinv(A)*Y;
			yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		   [ymax, ii] = max(yfit);
		   xmax = xfit(ii);
			disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      if 0
			if kfit(3) < 0
		      yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
		      [ymax, ii] = max(yfit);
		      xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);
	      else
	      	xmax = 0;
	      	ymax = Y(15);
	      	kfit = [mean(Y), 0, 0];
	      	disp(["No maximum from fit. Return x = ", num2str(xmax)]);
      	end
      	end
			vary_positron_octupole_knob${knob}${axis}(xmax);

			# save scan result:
			fname = [$save_dir, "/Res_positron_machine_",num2str($machine),"_knob${knob}${axis}.dat"];
			#eval(["save ", fname," Res kfit ymax xmax"]);
		}

		Octave {
			Lmeas = Lmeas + 4*29;
		}
	}
}
