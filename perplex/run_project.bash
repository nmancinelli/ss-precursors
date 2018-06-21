#!/bin/bash
#
project=$1
dat_file=$2
soln_file=$3
#
mkdir $project
cd $project
ln -s ../src .
ln -s ../*.dat .
ln -s ../$soln_file .
rm *.tab

src/vertex << ! > vertex.out
$dat_file
!
#
src/werami << ! > werami.out
$dat_file
2
2
n
13
n
14
n
0
y
10  150000
800 2200
100 100
0
!

python ../make_funs.py

cd ..
