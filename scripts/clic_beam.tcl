set ch {1.0}

#
# Choice of Structure
#

#array set structure {a 3e-3 g 6.5883e-3 l 8.33333e-3 delta 0.29 delta_g -0.915e-3}
array set structure {a 2.75e-3 g 7e-3 l 8.33333e-3 delta 0.145 delta_g 0.333e-3}

#
# Loading resistive wall wakefield
#

#set fwake [open res.wake]
#for {set i 0} {$i<12} {incr i} {
#    gets $fwake line
#}
#set resw ""
#gets $fwake line
#while {![eof $fwake]} {
#    lappend resw "[lindex $line 0] [expr 1e-15*[lindex $line 6]]"
#    gets $fwake line
#}
#close $fwake
#SplineCreate resist_wake $resw

#
# transverse wakefield
# s is given in micro metres
# return value is in V/pCm^2
#

proc w_transv {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set tmp [expr $g/$l]
    set alpha [expr 1.0-0.4648*sqrt($tmp)-(1.0-2.0*0.4648)*$tmp]
    set s0 [expr $g/8.0*pow($a/($l*$alpha),2)]
    return [expr 4.0*377.0*3e8*$s0*1e-12/(acos(-1.0)*pow($a,4))*(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]
}

proc w_transv {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set s0 [expr 0.169*pow($a,1.79)*pow($g,0.38)*pow($l,-1.17)]
    return [expr 4.0*377.0*3e8*$s0*1e-12/(acos(-1.0)*pow($a,4))*(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]
}

proc w_transv {s} {
    global structure
    set l $structure(l)
    set a $structure(a)

    set g $structure(g)
    set s0 [expr 0.169*pow($a,1.79)*pow($g,0.38)*pow($l,-1.17)]
    set tmp [expr 4.0*377.0*3e8*$s0*1e-12/\
		 (acos(-1.0)*pow($a,4))\
		 *(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]

    set g [expr $structure(g)+0.5*$structure(delta_g)]
    set s0 [expr 0.169*pow($a*(1.0+0.5*$structure(delta)),1.79)*pow($g,0.38)*pow($l,-1.17)]
    set tmp [expr $tmp+4.0*377.0*3e8*$s0*1e-12/\
		 (acos(-1.0)*pow($a*(1.0+0.5*$structure(delta)),4))\
		 *(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]

    set g [expr $structure(g)-0.5*$structure(delta_g)]
    set s0 [expr 0.169*pow($a*(1.0-0.5*$structure(delta)),1.79)*pow($g,0.38)*pow($l,-1.17)]
    set tmp [expr $tmp+4.0*377.0*3e8*$s0*1e-12/\
		 (acos(-1.0)*pow($a*(1.0-0.5*$structure(delta)),4))\
		 *(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]

    set g [expr $structure(g)-$structure(delta_g)]
    set s0 [expr 0.169*pow($a*(1.0-$structure(delta)),1.79)*pow($g,0.38)*pow($l,-1.17)]
    set tmp [expr $tmp+0.5*4.0*377.0*3e8*$s0*1e-12/\
		 (acos(-1.0)*pow($a*(1.0-$structure(delta)),4))\
		 *(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]

    set g [expr $structure(g)+$structure(delta_g)]
    set s0 [expr 0.169*pow($a*(1.0+$structure(delta)),1.79)*pow($g,0.38)*pow($l,-1.17)]
    set tmp [expr $tmp+0.5*4.0*377.0*3e8*$s0*1e-12/\
		 (acos(-1.0)*pow($a*(1.0+$structure(delta)),4))\
		 *(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]

    set tmp [expr 0.25*$tmp]

#    set tmp [expr $tmp+0.333*1e3*[resist_wake [expr $s*1e-6]]]
##    set tmp [expr $tmp+10*0.333*1e3*[resist_wake [expr $s*1e-6]]]

    return $tmp
}

#puts [w_transv 10.0]

#exit

#
# longitudinal wakefield
# s is given in micro metres
# return value is in V/pCm
#

proc w_long_x {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set s0 [expr 0.169*pow($a,1.79)*pow($g,0.38)*pow($l,-1.17)]
#    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
puts $s0
    return [expr 377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]
}

proc w_long_x {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set tmp [expr $g/$l]
    set alpha [expr 1.0-0.4648*sqrt($tmp)-(1.0-2.0*0.4648)*$tmp]
    set s0 [expr $g/8.0*pow($a/($l*$alpha),2)]
    return [expr 1.05*377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]
}

proc w_long {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
    return [expr 377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]
}

proc w_long {s} {
    global structure
    set l $structure(l)

    set a $structure(a)
    set g $structure(g)
    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
    set tmp [expr 377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]

    set a [expr $structure(a)*(1.0-0.5*$structure(delta))]
    set g [expr $structure(g)-0.5*$structure(delta_g)]
    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
    set tmp [expr $tmp+377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]

    set a [expr $structure(a)*(1.0+0.5*$structure(delta))]
    set g [expr $structure(g)+0.5*$structure(delta_g)]
    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
    set tmp [expr $tmp+377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]

    set a [expr $structure(a)*(1.0-$structure(delta))]
    set g [expr $structure(g)-$structure(delta_g)]
    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
    set tmp [expr $tmp+0.5*377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]

    set a [expr $structure(a)*(1.0+$structure(delta))]
    set g [expr $structure(g)+$structure(delta_g)]
    set s0 [expr 0.41*pow($a,1.8)*pow($g,1.6)*pow($l,-2.4)]
    set tmp [expr $tmp+0.5*377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]

    return [expr $tmp/4.0]
}
