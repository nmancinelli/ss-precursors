#!/bin/bash
lat=$1
lon=$2
#
echo $lat $lon
./GDM52_dispersion << ! > /dev/null
77
2
$lat $lon
99
!
#
