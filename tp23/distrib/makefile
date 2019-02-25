CC=mpicc
IDIR=.
CFLAGS=-I$(IDIR) -lm
EXEC=mandel
nb=50

DEPS=$(wildcard *.h)
SOURCE=$(wildcard *.c)
OBJECTS=$(SOURCE:%.c=%.o)

arg1=800 800 0.35 0.355 0.353 0.358 200
arg2=800 800 -0.736 -0.184 -0.735 -0.183 500
arg3=800 800 -0.736 -0.184 -0.735 -0.183 300
arg4=800 800 -1.48478 0.00006 -1.48440 0.00044 100
arg5=800 800 -1.5 -0.1 -1.3 0.1 10000
arg6=800 800 -1.5 -0.1 -1.3 0.1 100000

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $^ $(CFLAGS)

$(EXEC): $(SOURCE)
	$(CC) -O3 -g -o $@ $^ $(CFLAGS)


clean:
	rm $(EXEC)

run:
	make run1

run1:
	mpirun -n $(nb) -hostfile hostfile $(EXEC) $(arg1)

run2:
	mpirun -n $(nb) -hostfile hostfile $(EXEC) $(arg2)

run3:
	mpirun -n $(nb) -hostfile hostfile $(EXEC) $(arg3)

run4:
	mpirun -n $(nb) -hostfile hostfile $(EXEC) $(arg4)

run5:
	mpirun -n $(nb) -hostfile hostfile $(EXEC) $(arg5)

run6:
	echo $(arg6)
	mpirun -n $(nb) -hostfile hostfile $(EXEC) $(arg6)
