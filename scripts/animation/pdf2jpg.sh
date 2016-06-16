#!/bin/bash

echo "Converting pdf to jpg, waiting..."
index=0
for pdffile in pic/*.pdf; do
  jpgfile=${pdffile/.pdf/.jpg}
	index=$index+1
  convert $pdffile $jpgfile &
	if (($index % 16 == 0));
	then
		wait;
	fi
done
wait
