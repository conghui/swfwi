main: add.o wrap.o main.o
	swaCC -o main  main.o wrap.o add.o

add.o : add.c add.h
	swacc -c add.c

wrap.o : wrap.cpp wrap.h
	mpiCC -host -c wrap.cpp

main.o : main.cpp
	mpiCC -host -c main.cpp

run:
	/usr/sw-mpp/bin/bsub -I -b -m 1 -p -q q_sw_acc -host_stack 1024 -share_size 7000 -n 1 -cgsp 64 ./main

clean:
	rm *.o main -f *.a *.out

