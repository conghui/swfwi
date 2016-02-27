#!/bin/bash

export OMP_NUM_THREADS=2

mpirun -np 2 ./mpifwi vin=smvel.rsf  shots=shots.rsf \
  grads=grads.rsf objs=objs.rsf illums=illums.rsf niter=10 vout=vsnaps.rsf
