#!/usr/bin/env bash

set -e
#source /apps/spack/share/spack/setup-env.sh
#spack load gcc

versions=(origin v1 v1a v2 v2a v3 v4)
cxx="g++"
cxxflags="-O3 -march=native -fopenmp -fopt-info-vec"

for ver in ${versions[@]}
do
    echo "====> Version ${ver}"
    sourcefile="rt_etd_${ver}.cxx"
    exefile="etd_${ver}"
    ${cxx} ${cxxflags} ${sourcefile} -o ${exefile}
done
