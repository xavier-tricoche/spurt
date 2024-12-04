#!/bin/sh

# executable name
exec_dir=/home/xmt/code/xavier/build/bin
exec=${exec_dir}/single_source_field

# relevant directories
aniso_dir=${HOME}/aniso/
in_dir=${aniso_dir}/data
rbf_dir=${aniso_dir}/rbf
dense_dir=${aniso_dir}/dense

# base names
rbf_base=rbfls_weights
dense_base=dense_rbfls

for name in `ls ${in_dir}`; 
do
    base=${name/.dat/}
    input=${in_dir}/${name}
    rbf=${rbf_dir}/${rbf_base}_${base}
    dense=${dense_dir}/dense_${base}.nrrd
    tt=${dense_dir}/dense_${base}_tt.nrrd
    nabla=${dense_dir}/dense_${base}_nabla.nrrd
    ${exec} -i ${input} -w ${rbf} -k r3 -a high -v -o ${dense} \
    -r 400 400 -b 99.75 105 26 32
    echo ${dense} " successfully exported (1/3)"
    unu slice -i ${dense} -a 0 -p 0 -o ${tt}
    echo ${tt} " successfully exported (2/3)"
    unu crop -i ${dense} -min 1 0 0 -max M M M -o ${nabla}
    echo ${nabla} " successfully exported (3/3)"
    echo "\
    
    "
done
