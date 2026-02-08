#!/usr/bin/env bash

set -e
#source /apps/spack/share/spack/setup-env.sh
#spack load gcc

# Arguments
versions=${1:-"origin"}
istart=${2:$((-2**16))}
iend=${3:-$((2**21))}

# Executable
exefile="etd_mpi_${ver}"

echo
echo "====> Version: ${ver}"
echo "====> exe    : ${exefile}"

for ((n=istart;n<=iend;n*=2))
do
    mpirun -np 32 ./${exefile} ${n}
done
