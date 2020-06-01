# emittances unit is 10-7 m
array set match {
    alpha_x 7.421558072e-07
    alpha_y -0.0001911453068
    beta_x  32.49994211
    beta_y  8.999448416
}

set match(emitt_x)  9.0 
set match(emitt_y)  0.2
set match(charge) 5.2e9
set charge $match(charge)
set match(sigma_z) 70.0
set match(phase) 0.0
set match(e_spread) -1.0

if { $beam_case == 8 } {
   set match(emitt_x)  8.5 
   set match(emitt_y)  0.08
}

if { $beam_case == 20 } {
   set match(emitt_x)  9.0 
   set match(emitt_y)  0.2
}

if { $beam_case == 30 } {
   set match(emitt_x)  9.5 
   set match(emitt_y)  0.3
}
