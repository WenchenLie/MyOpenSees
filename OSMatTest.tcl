proc MatTest {path N MatType para} {

	wipe
	model basic -ndm 2 -ndf 3
	node 1 0 0 0
	node 2 1 0 0
	uniaxialMaterial $MatType 1 {*}$para
	element truss 1 1 2 1 1
	fix 1 1 1 1
	fix 2 0 1 1
	timeSeries Path 1 -dt 1.0 -values $path
	pattern Plain 3 1 {
		sp 2 1 1.0
	}
	constraints Lagrange
	numberer RCM
	system BandGeneral
	test EnergyIncr 0.0001 80  
	algorithm Newton
	integrator LoadControl 0.0
	analysis Static
	set time [expr {[llength $path] - 1}]
	set dt [expr {double($time) / $N}]
	integrator LoadControl $dt
	set u {}
	set F {}
	for {set i 0} {$i < $N} {incr i} {
    	analyze 1
    	lappend u [expr [nodeDisp 2 1]]
    	lappend F [expr [eleForce 1 4]]
	}
	return [list $u $F]
}

# ---------------- input parameters ---------------
set path [list 0 5 -10 20 -20 30 -30 40 -40 0];  # 目标位移幅值
set N 600
set MatType "Steel02"
set para [list 100 100 0.02 18.5 0.925 0.15]

# ------------------ write txt --------------------
set result [MatTest $path $N $MatType $para]
set u [lindex $result 0]
set F [lindex $result 1]
set fileName1 "u.txt"
set fileName2 "F.txt"
set file1 [open $fileName1 "w"]
set file2 [open $fileName2 "w"]
set listStr1 [join $u "\n"]
set listStr2 [join $F "\n"]
puts $file1 $listStr1
puts $file2 $listStr2
close $file1
close $file2

# ------------- run python script to plot -----------
exec python Plot.py



