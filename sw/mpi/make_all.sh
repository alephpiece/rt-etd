#!/bin/bash

set -e

versions=(origin)
cc="mpicc"
cflags="-O3 -msimd"
clibs="-lm"

for ver in ${versions[@]}
do
    echo "====> Version ${ver}"
    sourcefile="rt_etd_mpi_${ver}.c"
    exefile="etd_mpi_${ver}"
    ${cc} ${cflags} ${sourcefile} ${clibs} -o ${exefile}
done
