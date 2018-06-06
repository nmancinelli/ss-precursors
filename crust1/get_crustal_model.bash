#!/bin/bash
#
lat=$1
lon=$2

./getCN1point << ! > tmp
$lat $lon
q
!
head -n 13 tmp | tail -n 9 > tmp2
cat tmp


