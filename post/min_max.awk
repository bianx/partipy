#!/bin/awk -f
BEGIN {min=1000; max=-1000}
{if ($3<min) {min=$3} if ($3>max) {max=$3} }
END {print "min=",min, "max=",max}
