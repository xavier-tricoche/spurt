# exit when any command fails
set -e

if [ $# -eq 0 ]; then
    echo "Need 2 arguments: input nrrd file and output filename"
    exit 1
fi

echo 'vprobe -i $1 -k vector -q jac -o J.nrrd'
vprobe -i $1 -k vector -q jac -o J.nrrd
echo 'unu crop -i J.nrrd -min 0 0 0 0 -max 2 M M M -o Jrow1.nrrd'
unu crop -i J.nrrd -min 0 0 0 0 -max 2 M M M -o Jrow1.nrrd
echo 'unu crop -i J.nrrd -min 3 0 0 0 -max 5 M M M -o Jrow2.nrrd'
unu crop -i J.nrrd -min 3 0 0 0 -max 5 M M M -o Jrow2.nrrd
echo 'unu crop -i J.nrrd -min 6 0 0 0 -max 8 M M M -o Jrow3.nrrd'
unu crop -i J.nrrd -min 6 0 0 0 -max 8 M M M -o Jrow3.nrrd

echo 'unu 2op x Jrow1.nrrd Jrow1.nrrd | unu project -a 0 -m sum -o C11.nrrd'
unu 2op x Jrow1.nrrd Jrow1.nrrd | unu project -a 0 -m sum -o C11.nrrd
echo 'unu 2op x Jrow1.nrrd Jrow2.nrrd | unu project -a 0 -m sum -o C12.nrrd'
unu 2op x Jrow1.nrrd Jrow2.nrrd | unu project -a 0 -m sum -o C12.nrrd
echo 'unu 2op x Jrow1.nrrd Jrow3.nrrd | unu project -a 0 -m sum -o C13.nrrd'
unu 2op x Jrow1.nrrd Jrow3.nrrd | unu project -a 0 -m sum -o C13.nrrd
echo 'unu 2op x Jrow2.nrrd Jrow2.nrrd | unu project -a 0 -m sum -o C22.nrrd'
unu 2op x Jrow2.nrrd Jrow2.nrrd | unu project -a 0 -m sum -o C22.nrrd
echo 'unu 2op x Jrow2.nrrd Jrow3.nrrd | unu project -a 0 -m sum -o C23.nrrd'
unu 2op x Jrow2.nrrd Jrow3.nrrd | unu project -a 0 -m sum -o C23.nrrd
echo 'unu 2op x Jrow3.nrrd Jrow3.nrrd | unu project -a 0 -m sum -o C33.nrrd'
unu 2op x Jrow3.nrrd Jrow3.nrrd | unu project -a 0 -m sum -o C33.nrrd

echo 'unu join -i C11.nrrd C12.nrrd C13.nrrd C22.nrrd C23.nrrd C33.nrrd -a 0 -incr | unu pad -min 0 0 0 0 -max M+1 M M M -b pad -v 1 | unu shuffle -a 0 -p 6 0 1 2 3 4 5 -o C7.nrrd'
unu join -i C11.nrrd C12.nrrd C13.nrrd C22.nrrd C23.nrrd C33.nrrd -a 0 -incr | unu pad -min 0 0 0 0 -max M+1 M M M -b pad -v 1 | unu shuffle -a 0 -p 6 0 1 2 3 4 5 -o C7.nrrd

echo 'vprobe -i C7.nrrd -k tensor -q eval0 -o lmax.nrrd'
vprobe -i C7.nrrd -k tensor -q eval0 -o lmax.nrrd

if [ -z $2 ]
then
    ftlename="${1}_ftle.nrrd"
else
    ftlename=${2}
fi

echo 'unu 1op log2 -i lmax.nrrd -o ' ${ftlename}
unu 1op log2 -i lmax.nrrd | unu 2op exists - 0 -o $ftlename

echo 'cleanup'
rm -f J.nrrd Jrow1.nrrd Jrow2.nrrd Jrow3.nrrd C11.nrrd C12.nrrd C13.nrrd C22.nrrd C23.nrrd C33.nrrd C7.nrrd lmax.nrrd
