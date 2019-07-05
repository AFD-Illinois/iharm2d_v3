#!/bin/sh
#
# do a series of runs at factor-of-two increasing resolution
#
lab=A   # label for this series of runs
#
for res in 64 128 256 512 1024
do
	# stick the resolution in as N1 in the script
	scr="s/N1RES/$res/"
	sed $scr < decs.template > decs.h
	#
	make newrun
	#
	./harm &> harm.out
	#
	resdir="run$lab$res"
	mkdir $resdir
	mv harm.out images dumps $resdir
	cp *.c *.h $resdir
	#
done
