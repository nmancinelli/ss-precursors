#!/bin/bash
#
lat=$1
lon=$2
#
./getpoint_map.1 << !
$lon $lat
!
#cat point.1
