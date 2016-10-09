#!/bin/bash

data=shots.rsf
mean=0
variance=1

nshots=`sfget parform=n n3 < $data`
echo "nshots: $nshots"

subfiles=""

for ((i=0; i < $nshots; i++)); do
  echo "processing shot: $i"
  sfwindow < $data f3=$i n3=1 | sfnoise seed=$i var=$variance mean=$mean > data-noise-$i.rsf
  subfiles="${subfiles} data-noise-$i.rsf"
done

sfcat axis=3 $subfiles > shots-noise.rsf --out=`pwd`/rsf/shots-noise.rsf

sfrm $subfiles
