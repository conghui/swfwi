#!/bin/bash
index=1
for((i=0;$i<300;i+=1));
do
	j=`expr $i \+ 1`
	jump=`expr $i / 30`
	begin=`expr $jump \* 30`
	judge=`echo \($j \- $begin\) \% \($jump \+ 1\) | bc` 
	if [ $judge = 0 ] ; then
		echo "Producing ${j}th velocity"
		k=`printf "%03d" $index`
		index=`expr $index \+ 1`
		sfwindow f3=$i n3=1 < vsnaps.rsf > pic/$k.rsf
		sfgrey title="Velocity $j" color=j allpos=y pclip=100 bias=1500 gainpanel=1 scalebar=y barreverse=y barunit=m/s barlabel=Velocity < pic/$k.rsf > pic/$k.vpl
		vpconvert format=pdf pic/$k.vpl > /dev/null
	fi
done
