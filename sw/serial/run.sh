#!/bin/bash

set -e

# Arguments
# versions: origin, v1, v2, etc. (special: all)
versions=${1:-"origin"}
istart=${2:-$((2**7))}
iend=${3:-$((2**16))}

if [[ ${versions} == *"all"* ]]
then
    versions=(origin v1 v2)
fi

echo -e "\nThe following versions will be run:"
echo -e "\t${versions[@]}"

logbasedir=log/$(date "+%Y%m%d")

for ver in ${versions[@]}
do
    # Executable
    exefile="etd_${ver}"

    # Log directory
    logdir=${logbasedir}/${exefile}
    mkdir -p ${logdir} &> /dev/null

    echo
    echo "====> Version: ${ver}"
    echo "====> exe    : ${exefile}"
    echo "====> log dir: ${logdir}"

    for ((n=istart;n<=iend;n*=2))
    do
        echo "====> N      : ${n}"
        bsub -q q_sw_share -n 1 -o "${logdir}/N${n}_$(date "+%H%M").txt" ./${exefile} ${n}
    done
done
