#!/bin/bash

set -e

# Arguments
versions=${1:-"origin"}
istart=${2:-$((2**6))}
iend=${3:-$((2**12))}

if [[ ${versions} == *"all"* ]]
then
    versions=(origin)
fi

echo -e "\nThe following versions will be run:"
echo -e "\t${versions[@]}"

logbasedir=log/$(date "+%Y%m%d")

for ver in ${versions[@]}
do

    # Executable
    exefile="etd_mpi_${ver}"

    # Log directory
    logdir=${logbasedir}/${exefile}
    mkdir -p ${logdir} &> /dev/null

    echo
    echo "====> Version: ${ver}"
    echo "====> exe    : ${exefile}"
    echo "====> log dir: ${logdir}"

    #n=$((2**23))    # strong scaling
    n=8192  # weak scaling

    for ((np=istart;np<=iend;np*=2))
    do
        #N=${n}  # strong scaling
        N=$((n*np)) # weak scaling
        nodes=$(((np+3)/4))

        echo
        echo "====> N               : ${N}"
        echo "====> n per processes : ${n}"
        echo "====> nodes           : ${nodes}"
        echo "====> processes       : ${np} (host)"
        bsub -q q_sw_share -n ${np} -share_size 1024 -o "${logdir}/N${N}_n${n}_np${np}_nodes${nodes}_$(date "+%H%M%S").txt" ./${exefile} ${N}
    done
done
