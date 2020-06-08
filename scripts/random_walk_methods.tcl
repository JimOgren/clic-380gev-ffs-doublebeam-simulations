#############################################################################
Octave {
	global IPind Sind Oind
	IPind = IP;
	Sind = MIsext;
	Oind = MIoct;
	Lmeas = 0; % initialize luminosity measurement counter

	# create file for storing all luminosity measurements
	lumifile = fopen("LuminosityData.dat","w"); fclose(lumifile);

	function [] = write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
		sx1 = beam_size(1);		sy1 = beam_size(2);
		sx2 = beam_size(3);		sy2 = beam_size(4);
		e1 = beam_strahl(1);		e2 = beam_strahl(1);
		sxp1 = beam_strahl(3); 	syp1 = beam_strahl(4);
		sxp2 = beam_strahl(5); 	syp2 = beam_strahl(6);
		vx1 = beam_defl(1);		vy1 = beam_defl(2);
		vx2 = beam_defl(3);		vy2 = beam_defl(4);
		if nargin == 5
			E1 = -1; E2 = -1;
			Ne1 = -1; Np1 = -1;
			Ne2 = -1; Np2 = -1;
		else
			E1 = pairs(1); E2 = pairs(2);
			Ne1 = pairs(3); Np1 = pairs(4);
			Ne2 = pairs(5); Np2 = pairs(6);
		end

		# write to file:
      lumifile = fopen("LuminosityData.dat","a");
		fprintf(lumifile,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", L, Lpeak, sx1, sy1, sx2, sy2, e1, e2, sxp1, syp1, sxp2, syp2, vx1, vy1, vx2, vy2, E1, E2, Ne1, Np1, Ne2, Np2)
		fclose(lumifile);
	end
}

# RANDOM WALK METHODS
#############################################################################
#############################################################################
# Sextupoles only
Octave {
	global sel mov TEMP
	sel = 1;
	mov = 1;
	TEMP = [];

	# Move sextupoles
	###################################
	function [] = move_electron_sextupole(x)
		global sel mov Sind
		MIsext = Sind;
		ind = [MIsext, MIsext]';
		dir = [repmat("x", 6,1);repmat("y", 6,1)];
		s_temp = [];
		for i=1:length(sel)
			s_temp = [s_temp; placet_element_get_attribute("electron", ind(sel(i)), dir(sel(i)))];
		end

		for i = 1:length(sel)
			placet_element_set_attribute("electron",  ind(sel(i)), dir(sel(i)), s_temp(i) + x*mov(i));
		end
	end

	function [] = move_positron_sextupole(x)
		global sel mov Sind
		MIsext = Sind;
		ind = [MIsext, MIsext]';
		dir = [repmat("x", 6,1);repmat("y", 6,1)];
		s_temp = [];
		for i=1:length(sel)
			s_temp = [s_temp; placet_element_get_attribute("positron", ind(sel(i)), dir(sel(i)))];
		end

		for i = 1:length(sel)
			placet_element_set_attribute("positron",  ind(sel(i)), dir(sel(i)), s_temp(i) + x*mov(i));
		end
	end
}

# Procedure: random walk beamstrahlung maximizer
#############################################################################
#############################################################################
proc random_walk_beamstrahlung_power { subset } {
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder

	set NumOfPoints 7
	set MaxIterations 1500
	set count 1
   set flag 1
   set BSmax 0
   set BSlim 18.0
   set range 15.0

	Octave {
		printf("\n=================================================================\n")
		printf("=================================================================\n\n")
      printf("Starting Sextupole Random Walk - Beamstrahlung TOTAL POWER\n\n");
      printf("=================================================================\n")
      # check starting condition:
		[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
		[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);
		[L, Lpeak, beam_size, beam_strahl, beam_defl] = get_lumi(B1, B2);
		e1 = beam_strahl(1);
		e2 = beam_strahl(2);
   	Y_in = (e1+e2);
		printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
		printf("[e1, e2] \t= [%g, %g]\n",e1, e2);
		printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
		printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
		write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl)
   }

	while {$flag} {
		Octave {
			# Set random direction:
			sel = [];
			for j = 1:$subset
				temp = randi(12);
				if sum(temp==sel) < 1
					sel = [sel; temp];
				end
			end
			mov = randn(size(sel));
		}

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"

		if {rand() < 0.5} {
			puts "tune electron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_electron_sextupole(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"];
					eval(str1);
					eval(str2);

					move_electron_sextupole(-xvec(jk));
				end
				pause(5);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 0 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl];
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl)
				end
				Y = Res(:,8) + Res(:,9);
				[Res(:,1:2), Y, Res(:, 8:9), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_electron_sextupole(xmax);
				Tcl_Eval(["set BSmax ", num2str(ymax)]);
			}

		# tune positron beamline
		} else {
			puts "tune positron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_positron_sextupole(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
					eval(str1);
					eval(str2);

					move_positron_sextupole(-xvec(jk));
				end
				pause(10);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 0 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl];
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl)
				end
				Res;
				Y = Res(:,8) + Res(:,9);
				[Res(:,1:2), Y, Res(:, 8:9), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_positron_sextupole(xmax);
				Tcl_Eval(["set BSmax ", num2str(ymax)]);
			}
		}
		set count [expr $count+$NumOfPoints]
		puts "count = $count. BSmax = $BSmax"

		Octave {
			if ($count > $MaxIterations || $BSmax > $BSlim )
				Tcl_Eval("set flag 0");
			end
		}
	}

	Octave {
		printf("=================================================================\n\n")
		printf("Exiting beamstrahlung random walk...\n");
		printf("Total count = %i; BSmax = %g; \n", $count, $BSmax);
		printf("=================================================================\n")
		printf("=================================================================\n")
	}

	# remove temp folders:
	exec bash -c "rm -rf beams_temp"
	exec bash -c "rm -rf GP_output_temp"
	return $count
}



# Procedure: random walk pairs maximizer
#############################################################################
#############################################################################
proc random_walk_pairs { subset } {
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder

	set NumOfPoints 7
	set MaxIterations 2000
	set count 1
   set flag 1
   set Pmax 0
   set Plim 2500.0
   set range 10.0

	Octave {
		printf("\n=================================================================\n")
		printf("=================================================================\n\n")
      printf("Starting Sextupole Random Walk - PAIRS\n\n");
      printf("=================================================================\n")
      # check starting condition:
		[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
		[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);
		[L, Lpeak, beam_size,beam_strahl, beam_defl, pairs] = get_lumi_full(B1, B2);
		p1 = pairs(1);
		p2 = pairs(2);
   	Y_in = (p1+p2);
		printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
		printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
		printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
		printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
		write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
   }

	while {$flag} {
		Octave {
			# Set random direction:
			sel = [];
			for j = 1:$subset
				temp = randi(12);
				if sum(temp==sel) < 1
					sel = [sel; temp];
				end
			end
			mov = randn(size(sel));
		}

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"

		if {rand() < 0.5} {
			puts " "
			puts "tune electron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_electron_sextupole(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"];
					eval(str1);
					eval(str2);

					move_electron_sextupole(-xvec(jk));
				end
				pause(5);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
					p1 = pairs(1);
					p2 = pairs(2);
					printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
					printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
					printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
					printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
					printf("Y = %g at x = %g\n\n", (p1+p2), xvec(j));
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
				end
				Y = Res(:,18) + Res(:,19);
				[Res(:,1:2), Y, Res(:, 18:19), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_electron_sextupole(xmax);
				Tcl_Eval(["set Pmax ", num2str(ymax)]);
			}

		# tune positron beamline
		} else {
			puts " "
			puts "tune positron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_positron_sextupole(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
					eval(str1);
					eval(str2);

					move_positron_sextupole(-xvec(jk));
				end
				pause(10);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
					p1 = pairs(1);
					p2 = pairs(2);
					printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
					printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
					printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
					printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
					printf("Y = %g at x = %g\n\n", (p1+p2), xvec(j));
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
				end
				Res;
				Y = Res(:,18) + Res(:,19);
				[Res(:,1:2), Y, Res(:, 18:19), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_positron_sextupole(xmax);
				Tcl_Eval(["set Pmax ", num2str(ymax)]);
			}
		}
		set count [expr $count+$NumOfPoints]
		puts "count = $count. Pmax = $Pmax"

		Octave {
			if ($count > $MaxIterations || $Pmax > $Plim )
				Tcl_Eval("set flag 0");
			end
		}
	}

	Octave {
		printf("=================================================================\n\n")
		printf("Exiting pairs random walk...\n");
		printf("Total count = %i; Pmax = %g; \n", $count, $Pmax);
		printf("=================================================================\n")
		printf("=================================================================\n")
	}

	# remove temp folders:
	exec bash -c "rm -rf beams_temp"
	exec bash -c "rm -rf GP_output_temp"

	return $count

}



# Procedure: random walk pairs maximizer, LONG VERSION
#############################################################################
#############################################################################
proc random_walk_pairs_long { subset } {
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder

	set NumOfPoints 7
	set MaxIterations 200
	set count 1
   set flag 1
   set Pmax 0
   set Plim 6200.0
   set range 5.0

	Octave {
		printf("\n=================================================================\n")
		printf("=================================================================\n\n")
      printf("Starting Sextupole Random Walk - PAIRS LONG VERSION\n\n");
      printf("=================================================================\n")
      # check starting condition:
		[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
		[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);
		[L, Lpeak, beam_size,beam_strahl, beam_defl, pairs] = get_lumi_full(B1, B2);
		p1 = pairs(1);
		p2 = pairs(2);
   	Y_in = (p1+p2);
		printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
		printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
		printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
		printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
   }

	while {$flag} {
		Octave {
			# Set random direction:
			sel = [];
			for j = 1:$subset
				temp = randi(12);
				if sum(temp==sel) < 1
					sel = [sel; temp];
				end
			end
			mov = randn(size(sel));
		}

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"

		if {rand() < 0.5} {
			puts " "
			puts "tune electron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_electron_sextupole(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"];
					eval(str1);
					eval(str2);

					move_electron_sextupole(-xvec(jk));
				end
				pause(5);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
					p1 = pairs(1);
					p2 = pairs(2);
					printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
					printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
					printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
					printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
					printf("Y = %g at x = %g\n\n", (p1+p2), xvec(j));
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
				end
				Y = Res(:,18) + Res(:,19);
				[Res(:,1:2), Y, Res(:, 18:19), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_electron_sextupole(xmax);
				Tcl_Eval(["set Pmax ", num2str(ymax)]);
			}

		# tune positron beamline
		} else {
			puts " "
			puts "tune positron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_positron_sextupole(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
					eval(str1);
					eval(str2);

					move_positron_sextupole(-xvec(jk));
				end
				pause(10);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
					p1 = pairs(1);
					p2 = pairs(2);
					printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
					printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
					printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
					printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
					printf("Y = %g at x = %g\n\n", (p1+p2), xvec(j));
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
				end
				Res;
				Y = Res(:,18) + Res(:,19);
				[Res(:,1:2), Y, Res(:, 18:19), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_positron_sextupole(xmax);
				Tcl_Eval(["set Pmax ", num2str(ymax)]);
			}
		}
		set count [expr $count+$NumOfPoints]
		puts "count = $count. Pmax = $Pmax"

		Octave {
			if ($count > $MaxIterations || $Pmax > $Plim )
				Tcl_Eval("set flag 0");
			end
		}
	}

	Octave {
		printf("=================================================================\n\n")
		printf("Exiting pairs random walk...\n");
		printf("Total count = %i; Pmax = %g; \n", $count, $Pmax);
		printf("=================================================================\n")
		printf("=================================================================\n")
	}

	# remove temp folders:
	exec bash -c "rm -rf beams_temp"
	exec bash -c "rm -rf GP_output_temp"

	return $count
}





# RANDOM WALK METHODS - SEXTUPOLES AND QUADRUPOLES TOGETHER
#############################################################################
#############################################################################
Octave {
	# Move sextupoles and quadrupoles
	###################################
	function [] = move_electron_sq(x)
		global sel mov
      MIsext = [48 59 70 79 210 220];
      QI = [2 6 10 14 18 22 37 41 45 51 55 63 67 76 82 97 112 127 214 224];
      ind = [MIsext, MIsext, QI, QI]';
      dir = [repmat("x", 6,1);repmat("y", 6,1);repmat("x", 20,1);repmat("y", 20,1)];
		s_temp = [];
		for i=1:length(sel)
			s_temp = [s_temp; placet_element_get_attribute("electron", ind(sel(i)), dir(sel(i)))];
		end

		for i = 1:length(sel)
			placet_element_set_attribute("electron",  ind(sel(i)), dir(sel(i)), s_temp(i) + x*mov(i));
		end
	end

	function [] = move_positron_sq(x)
		global sel mov Sind
      MIsext = [48 59 70 79 210 220];
      QI = [2 6 10 14 18 22 37 41 45 51 55 63 67 76 82 97 112 127 214 224];
      ind = [MIsext, MIsext, QI, QI]';
      dir = [repmat("x", 6,1);repmat("y", 6,1);repmat("x", 20,1);repmat("y", 20,1)];
		s_temp = [];
		for i=1:length(sel)
			s_temp = [s_temp; placet_element_get_attribute("positron", ind(sel(i)), dir(sel(i)))];
		end

		for i = 1:length(sel)
			placet_element_set_attribute("positron",  ind(sel(i)), dir(sel(i)), s_temp(i) + x*mov(i));
		end
	end
}


# Procedure: random walk pairs maximizer, sextupoles and quadrupoles
#############################################################################
#############################################################################
proc random_walk_pairs_sq { subset } {
	exec bash -c "rm -rf temp_*" # remove temp folders
	exec bash -c "mkdir -p GP_output_temp/" # create temp folder
	exec bash -c "mkdir -p beams_temp/" # create temp folder

	set NumOfPoints 7
	set MaxIterations 1000
	set count 1
   set flag 1
   set Pmax 0
   set Plim 6200.0
   set range 0.5

	Octave {
		printf("\n=================================================================\n")
		printf("=================================================================\n\n")
      printf("Starting Sextupole and Quadrupole Random Walk - PAIRS LONG VERSION\n\n");
      printf("=================================================================\n")
      # check starting condition:
		[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
		[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);
		[L, Lpeak, beam_size,beam_strahl, beam_defl, pairs] = get_lumi_full(B1, B2);
		p1 = pairs(1);
		p2 = pairs(2);
   	Y_in = (p1+p2);
		printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
		printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
		printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
		printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
   }

	while {$flag} {
		Octave {
			# Set random direction:
			sel = [];
			for j = 1:$subset
				temp = randi(52);
				if sum(temp==sel) < 1
					sel = [sel; temp];
				end
			end
			mov = randn(size(sel));
		}

		# clean temp folders:
		exec bash -c "rm -f beams_temp/*"
		exec bash -c "rm -f GP_output_temp/*"

		if {rand() < 0.5} {
			puts " "
			puts "tune electron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_electron_sq(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"];
					eval(str1);
					eval(str2);

					move_electron_sq(-xvec(jk));
				end
				pause(5);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
					p1 = pairs(1);
					p2 = pairs(2);
					printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
					printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
					printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
					printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
					printf("Y = %g at x = %g\n\n", (p1+p2), xvec(j));
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
				end
				Y = Res(:,18) + Res(:,19);
				[Res(:,1:2), Y, Res(:, 18:19), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_electron_sq(xmax);
				Tcl_Eval(["set Pmax ", num2str(ymax)]);
			}

		# tune positron beamline
		} else {
			puts " "
			puts "tune positron beamline"
			Octave {
				xvec = linspace(-$range, $range, $NumOfPoints);
				for jk=1:length(xvec)
					move_positron_sq(+xvec(jk));

					# track:
					[E1,B1] = placet_test_no_correction("electron", "beam0t", "None", 1, 0, IPind);
					[E2,B2] = placet_test_no_correction("positron", "beam0t", "None", 1, 0, IPind);

					# save beams to temp-folder
					str1 = ["save_beam(""beams_temp/electron_",num2str(jk),".dat"", B1);"];
					str2 = ["save_beam(""beams_temp/positron_",num2str(jk),".dat"", B2);"]	;
					eval(str1);
					eval(str2);

					move_positron_sq(-xvec(jk));
				end
				pause(10);
			}

			# execute GUINEA-PIG in parallel:
			exec bash -c "for i in `seq 1 1 $NumOfPoints`  ; do ./run_guinea_parallel.sh 0 1 \$i & done"

			Octave {
				pause(5);
				Res = [];
				for j = 1:length(xvec)
					eval(["load GP_output_temp/GP_output_",num2str(j),".dat"])
					Res = [Res; xvec(j), L, Lpeak, beam_size, beam_strahl, beam_defl, pairs];
					p1 = pairs(1);
					p2 = pairs(2);
					printf("[L, Lpeak] \t= [%g, %g]\n",L,Lpeak);
					printf("[p1, p2] \t= [%g, %g]\n",p1, p2);
					printf("[sx1, sy1] \t= [%g, %g]\n",beam_size(1),beam_size(2));
					printf("[sx2, sy2] \t= [%g, %g]\n",beam_size(3),beam_size(4));
					printf("Y = %g at x = %g\n\n", (p1+p2), xvec(j));
					write_to_lumifile(L, Lpeak, beam_size, beam_strahl, beam_defl, pairs)
				end
				Res;
				Y = Res(:,18) + Res(:,19);
				[Res(:,1:2), Y, Res(:, 18:19), Res(:, 4:7)]

				# make parabolic fit:
				xfit = linspace(xvec(1),xvec(end),200);
				A = [ones(size(xvec')), xvec', (xvec').^2];
				kfit = pinv(A)*Y;
				yfit = kfit(1) + kfit(2)*xfit + kfit(3)*xfit.^2;
				[ymax, ii] = max(yfit);
				xmax = xfit(ii);
				disp(["Max value from fit = ", num2str(ymax), " at x = ", num2str(xmax)]);

				move_positron_sq(xmax);
				Tcl_Eval(["set Pmax ", num2str(ymax)]);
			}
		}
		set count [expr $count+$NumOfPoints]
		puts "count = $count. Pmax = $Pmax"

		Octave {
			if ($count > $MaxIterations || $Pmax > $Plim )
				Tcl_Eval("set flag 0");
			end
		}
	}

	Octave {
		printf("=================================================================\n\n")
		printf("Exiting pairs random walk...\n");
		printf("Total count = %i; Pmax = %g; \n", $count, $Pmax);
		printf("=================================================================\n")
		printf("=================================================================\n")
	}

	# remove temp folders:
	exec bash -c "rm -rf beams_temp"
	exec bash -c "rm -rf GP_output_temp"

	return $count

}
