# exit when any command fails
set -e

if [ $# -eq 0 ]; then
    echo "Need 2 arguments: input nrrd file and output filename"
    exit 1
fi


# first pad the 2d flow map to make it a (trivial) 3d flow map in 2D
echo 'unu pad -i $1 -min 0 0 0 -max M+1 M M -o __fmap32d.nrrd'
unu pad -i $1 -min 0 0 0 -max M+1 M M -o __fmap32d.nrrd

# next duplicate the flow map along a 3rd spatial dimension
echo 'unu join -i __fmap32d.nrrd __fmap32d.nrrd -a 3 -incr | unu axinfo -a 3 -sp 1 -c node -o __fmap3d.nrrd'
unu join -i __fmap32d.nrrd __fmap32d.nrrd -a 3 -incr | unu axinfo -a 3 -sp 1 -c node -o __fmap3d.nrrd

echo 'vprobe -i __fmap3d.nrrd -k vector -q jac -o __J.nrrd'
vprobe -i __fmap3d.nrrd -k vector -q jac -o __J.nrrd
echo 'unu crop -i __J.nrrd -min 0 0 0 0 -max 2 M M M -o __Jrow1.nrrd'
unu crop -i __J.nrrd -min 0 0 0 0 -max 2 M M M -o __Jrow1.nrrd
echo 'unu crop -i __J.nrrd -min 3 0 0 0 -max 5 M M M -o __Jrow2.nrrd'
unu crop -i __J.nrrd -min 3 0 0 0 -max 5 M M M -o __Jrow2.nrrd
echo 'unu crop -i __J.nrrd -min 6 0 0 0 -max 8 M M M -o __Jrow3.nrrd'
unu crop -i __J.nrrd -min 6 0 0 0 -max 8 M M M -o __Jrow3.nrrd

echo 'unu 2op x __Jrow1.nrrd __Jrow1.nrrd | unu project -a 0 -m sum -o __C11.nrrd'
unu 2op x __Jrow1.nrrd __Jrow1.nrrd | unu project -a 0 -m sum -o __C11.nrrd
echo 'unu 2op x __Jrow1.nrrd __Jrow2.nrrd | unu project -a 0 -m sum -o __C12.nrrd'
unu 2op x __Jrow1.nrrd __Jrow2.nrrd | unu project -a 0 -m sum -o __C12.nrrd
echo 'unu 2op x __Jrow1.nrrd __Jrow3.nrrd | unu project -a 0 -m sum -o __C13.nrrd'
unu 2op x __Jrow1.nrrd __Jrow3.nrrd | unu project -a 0 -m sum -o __C13.nrrd
echo 'unu 2op x __Jrow2.nrrd __Jrow2.nrrd | unu project -a 0 -m sum -o __C22.nrrd'
unu 2op x __Jrow2.nrrd __Jrow2.nrrd | unu project -a 0 -m sum -o __C22.nrrd
echo 'unu 2op x __Jrow2.nrrd __Jrow3.nrrd | unu project -a 0 -m sum -o __C23.nrrd'
unu 2op x __Jrow2.nrrd __Jrow3.nrrd | unu project -a 0 -m sum -o __C23.nrrd
echo 'unu 2op x __Jrow3.nrrd __Jrow3.nrrd | unu project -a 0 -m sum -o __C33.nrrd'
unu 2op x __Jrow3.nrrd __Jrow3.nrrd | unu project -a 0 -m sum -o __C33.nrrd

echo 'unu join -i __C11.nrrd __C12.nrrd __C13.nrrd __C22.nrrd __C23.nrrd __C33.nrrd -a 0 -incr | unu pad -min 0 0 0 0 -max M+1 M M M -b pad -v 1 | unu shuffle -a 0 -p 6 0 1 2 3 4 5 -o __C7.nrrd'
unu join -i __C11.nrrd __C12.nrrd __C13.nrrd __C22.nrrd __C23.nrrd __C33.nrrd -a 0 -incr | unu pad -min 0 0 0 0 -max M+1 M M M -b pad -v 1 | unu shuffle -a 0 -p 6 0 1 2 3 4 5 -o __C7.nrrd

echo 'vprobe -i __C7.nrrd -k tensor -q eval0 | unu slice -a 2 -p 0 -o __lmax.nrrd'
vprobe -i __C7.nrrd -k tensor -q eval0 | unu slice -a 2 -p 0 -o __lmax.nrrd

if [ -z $2 ]
then
    ftlename="${1}_ftle.nrrd"
else
    ftlename=${2}
fi

echo 'unu 1op log2 -i __lmax.nrrd | unu slice -a 3 -p 0 -o ' ${ftlename}
unu 1op log2 -i __lmax.nrrd | unu 2op exists - 0 -o $ftlename

echo 'cleanup'
rm -f __fmap32d.nrrd __fmap23d.nrrd __fmap3d.nrrd __J.nrrd __Jrow1.nrrd __Jrow2.nrrd __Jrow3.nrrd __C11.nrrd __C12.nrrd __C13.nrrd __C22.nrrd __C23.nrrd __C33.nrrd __C7.nrrd __lmax.nrrd
