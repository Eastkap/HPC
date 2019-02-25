CC=mpicc
IDIR=.
CFLAGS=-I$(IDIR) -lm
EXEC=mandel
nb=1

arg1=800 800 0.35 0.355 0.353 0.358 200
arg2=800 800 -0.736 -0.184 -0.735 -0.183 500
arg3=800 800 -0.736 -0.184 -0.735 -0.183 300
arg4=800 800 -1.48478 0.00006 -1.48440 0.00044 100
arg5=800 800 -1.5 -0.1 -1.3 0.1 10000
arg6=8000 8000 -1.5 -0.1 -1.3 0.1 10000



all0:
	mpicc -o mandel mandel.c

all1:
	mpicc -O1 -o  mandel mandel.c

all2:
	mpicc -O2 -o  mandel mandel.c

all3:
	mpicc -O3 -o  mandel mandel.c

run1:
	mpirun $(EXEC) $(arg1)

run2:
	mpirun  $(EXEC) $(arg2)

run3:
	mpirun  $(EXEC) $(arg3)

run4:
	mpirun  $(EXEC) $(arg4)

run5:
	mpirun  $(EXEC) $(arg5)

run6:
	echo $(arg6)
	mpirun $(EXEC) $(arg6)

runone:
	make run1
	make run2
	make run3
	make run4
	make run5
	make run6
