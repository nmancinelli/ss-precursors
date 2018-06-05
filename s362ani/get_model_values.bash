#!/bin/bash
#
lat=$1
lon=$2
depth=$3
cd PROGRAMS
bin/test_subshsv.exe << ! > ../tmp
$lat
$lon
$depth
!
cd ..
head -n 2 tmp | tail -n 1 > tmp2
