#!/bin/bash

 /usr/sw-mpp/bin/bsub -I -b -m 1 -p -q q_sw_err_slow  \
   -host_stack 1024 -share_size 7000 -n 1 -np 1 -cgsp 64 -o bsub.out -J "fd" \
   ./fd


