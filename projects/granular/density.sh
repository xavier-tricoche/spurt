#!/bin/sh

path=/scratch/xmt/data/NJIT/tapping/gamma_2.75-tpour=0.5-5taps/

for ((i=$1 ; i<$2 ; i+= $3)); do
t=`echo "scale=5; $i/100000" | bc -l`
input=$path/steps/steps_`printf "%06d" $i`.nrrd
output=$path/analysis/new_density_`printf "%06d" $i`.nrrd
h=`/home/xmt/code/xavier/build/bin/bump_height -t $t -r 0.5 -g 2.75`
/home/xmt/code/xavier/build/bin/nrrd2txt -i $input -o tmp.txt
/scratch/xmt/bin/voro++ -px -pz 0.02 0 0.24 $h 0.7 0 0.24 tmp.txt
n=`cat tmp.txt.vol | wc -l`
echo "processing" $input
echo "at time =" $t ", height =" $h
echo $n "points in output of Voronoi computation"
/scratch/xmt/bin/unu make -i tmp.txt.vol -e ascii -t float -s 5 $n -o tmp.nrrd
/home/xmt/code/xavier/build/bin/density2D -i tmp.nrrd -o tmp2.nrrd -r 260 700 -t $t -rel 0.5 -g 2.75 -ymax 0.7
/scratch/xmt/bin/unu flip -i tmp2.nrrd -a 1 -o $output
echo $output "exported"
done

