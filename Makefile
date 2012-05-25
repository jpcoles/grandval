CC=nvcc
CFLAGS=-Xcompiler "-Wall -O3 -fopenmp" --generate-code arch=compute_20,code=sm_21 --ptxas-options="-v"
LDFLAGS=-Xlinker "-lm -fopenmp -lpng"
SRC=$(wildcard *.cu)
OBJS=$(patsubst %.cu,%.o, $(SRC))

all: $(OBJS)
	echo $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o grandval

%.o: %.cu
	echo $(OBJS) $(SRC)
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
