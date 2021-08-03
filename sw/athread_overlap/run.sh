#!/bin/bash

n=${1:-4096}
bsub -b -q q_sw_share -n 1 -cgsp 64 -share_size 256 -I ./etd_athread ${n}
