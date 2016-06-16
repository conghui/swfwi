#!/bin/bash
if [ ! -d pic ] ; then
	mkdir pic
fi
./produce.sh
./pdf2jpg.sh
./jpg2mp4.sh
