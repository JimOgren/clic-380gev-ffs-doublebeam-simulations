#
# Select the main linac phases. Here, only one is chosen.
#

set ph1 -7.0
set ph2 4.0
set ph3 30.0

#
# Define wavelength
#

set lambda [expr 0.025]

#
# Define gradient
#

set gradient [expr 0.1]

#
# Define longrange wakefield
# uncomment next command if single bunch only
#

set cav_modes {}
set N_mode 1
for {set i 0} {$i<$N_mode} {incr i} {
    lappend cav_modes [expr 0.15/10.25]
    lappend cav_modes 1.0
    lappend cav_modes 3
}

WakeSet waketmp $cav_modes
set ampl [waketmp 0.15]

if {[llength [info vars scale]]==0} {set scale 1.0}

set cav_modes {}
set N_mode 1
for {set i 0} {$i<$N_mode} {incr i} {
    set line "1.0 1.0 1.0"
    lappend cav_modes [expr 0.15/10.25]
    lappend cav_modes [expr $scale*6.6e3/$ampl]
    lappend cav_modes 3
}

#
# use this list to create fields
#

WakeSet wakelong $cav_modes

#
# define structure
#

InjectorCavityDefine -lambda $lambda \
    -wakelong wakelong \


#MainPhases -phases "$ph1 $ph2 $ph3"
