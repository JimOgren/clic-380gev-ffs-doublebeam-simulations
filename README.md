# CLIC 380GeV FFS DoubleBeam Simulations

### CLIC FFS Simulations
These scripts simulates tuning procedures for the final-focus system (FFS) of the Compact Linear Collider (CLIC) particle accelerator. It is based on single-particle simulations for beam-based alignment (BBA) followed by sextupole alignment and knobs tuning based on a 1e5 multi-particle simulation. The luminosity is evaluated in a full beam-beam simulation. 

In this study, the two sides of the final-focus system were treated independently with independent imperfections. As guiding signals beamstrahlung and incoherent pairs instead of luminosity.

### Running the full simulation
The full simulation is defined in separate steps. To run the whole thing:
1) run /main/run_singlebeam_BBA.sh
2) run /main/run_

### Requirements
You need to have the tracking code PLACET and the beam-beam code GUINEA-PIG installed.

### Publication:

### Publication:
The results were presented in the paper: [J. Ogren et al., _Tuning the Compact Linear Collider 380 GeV final-focus system using realistic beam-beam signals_, Phys. Rev. Accel. Beams 23, 051002, May 2020.](https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.23.051002)

Written by Jim Ogren, 2018-2020, CERN, Geneva, Switzerland. 
