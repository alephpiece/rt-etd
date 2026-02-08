#!/bin/bash

set -e

versions=(origin v1 v2)
cc="sw5cc"
cflags="-host -O3 -msimd"
clibs="-lm"

for ver in ${versions[@]}
do
    echo "====> Version ${ver}"
    sourcefile="rt_etd_${ver}.c"
    exefile="etd_${ver}"
    ${cc} ${cflags} ${sourcefile} ${clibs} -o ${exefile}
done
