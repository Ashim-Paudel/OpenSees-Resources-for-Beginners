#Written: fmk 04/10

model testUniaxial

foreach {damage matTag fileOut} {nodamage 4 a damage 5 b} {

    set out [open $fileOut w]
    set plot [open plot$fileOut.m w]

    puts $plot "load $fileOut";
    puts $plot "hold off";
    puts $plot "plot($fileOut\(:,1\), $fileOut\(:,2\));";
    puts $plot "hold on; title(\"damage=$damage\"); grid on;";

    # material parameters
    set E 50
    set Fy 10
    set gap 0.2
    set eta 0.01

    uniaxialMaterial ElasticPPGap $matTag $E $Fy $gap $eta $damage

    uniaxialTest $matTag

    set strain 0.0
    set count 1
    foreach {strainIncr numIncr sym} {0.02 100 "x" -0.02 50 "o" 0.02 150 "r"} {
	for {set i 0} {$i < $numIncr} {incr i 1} {
	    set strain [expr $strain+$strainIncr]
	    strainUniaxialTest $strain
	    set stress [stressUniaxialTest]
	    set tangent [tangUniaxialTest]
	    puts $out "$strain $stress $tangent"
	}
	puts $plot "plot($fileOut\($count:[expr $count+$numIncr],1\), $fileOut\($count:[expr $count+$numIncr],2\), '$sym');"
	set count [expr $count+$numIncr-1]
    }

    close $out
    close $plot
}
