#!/usr/bin/env bash

set -e

versions=(origin)
cxx="mpicxx"
cxxflags="-O3 -fopenmp -fopt-info-vec"

for ver in ${versions[@]}
do
    echo "====> Version ${ver}"
    sourcefile="rt_etd_mpi_${ver}.cxx"
    exefile="etd_mpi_${ver}"
    ${cxx} ${cxxflags} ${sourcefile} -o ${exefile}
done
