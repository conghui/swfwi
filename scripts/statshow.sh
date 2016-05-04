#!/bin/bash

for s in absobjs.rsf norobjs.rsf l1norm.rsf l2norm.rsf; do
  echo -e "\n$s:"
  sfdisfil < $s col=7 | tail -n1
done
