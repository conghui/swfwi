#!/bin/bash
#set -x

if [[ $# -ne 1 ]]; then
  echo "please specify the version. (e.g. v1.0)"
  exit
fi

tag=$1
jobdirname=job
cd ..
rootdir=`pwd`

# move job to another dir
mv $jobdirname ../

cd ..                                 # go to top dir
targetdir=$rootdir-$tag
mv $rootdir $targetdir                # rename it
tar czf ${targetdir}.tgz ${targetdir##*/}/  # create a tarball
mv $targetdir $rootdir                # rename back
mv $jobdirname $rootdir               # move job back
