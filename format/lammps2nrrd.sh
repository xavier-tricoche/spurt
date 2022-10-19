#!/bin/sh

input=$1
echo "processing LAMMPS file" $input

base=$2
echo "exporting NRRD files with base name" $base

natoms=$3
echo "there are" $natoms "atoms per time step"

nsteps=$4
echo "there are" $nsteps "steps total"

for ((i=0 ; i<nsteps; i++)); do
	skip=$(echo "scale=6; $i*(9+$natoms)+9" | bc);
	echo "skipping" $skip "lines"
	unu make -i $input -t float -s 5 $natoms -ls $skip -e ascii -o ${base}_t=`printf "%04d" $i`.nrrd
	echo ${base}_t=`printf "%04d" $i`.nrrd "exported"
done