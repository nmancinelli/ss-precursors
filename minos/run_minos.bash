#!/bin/bash
#
# Run a command
#
file=$1

./minos_bran2 << EOF
$file.txt
out.txt
$file.tmp
1.d-8 100
3
1 3000 1. 200. 0 0
EOF
cat out.txt | grep " s " | awk {'print $4,  $5,  $7'} > $file.out

#Clean up
rm out.txt
rm $file.tmp

exit
