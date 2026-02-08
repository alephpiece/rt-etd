#!/usr/bin/env bash
#SBATCH -J singlecore
#SBATCH -p Balerion
#SBATCH -n 1
#SBATCH --exclusive

set -e
#source /apps/spack/share/spack/setup-env.sh
#spack load gcc

# Arguments
# versions: origin, v1, v2, etc. (special: all)
versions=${1:-"origin"}
istart=${2:-1}
iend=${3:-$((2**20))}

if [[ ${versions} == *"all"* ]]
then
    versions=(origin v1 v1a v2 v2a v3)
fi

echo -e "\nThe following versions will be run:"
echo -e "\t${versions[@]}"

for ver in ${versions[@]}
do
    # Executable
    exefile="etd_${ver}"

    echo
    echo "====> Version: ${ver}"
    echo "====> exe    : ${exefile}"

    for ((n=istart;n<=iend;n*=2))
    do
        ./${exefile} ${n}
    done
done
