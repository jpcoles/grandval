CC=nvcc
#CFLAGS=-Xcompiler "-Wall -O3 -fopenmp"  -v
LDFLAGS_PNG=$(shell libpng-config --ldflags)
CFLAGS_PNG=$(shell libpng-config --cflags)
CFLAGS=-g -G -Xcompiler "-Wall -O3 $(CFLAGS_PNG)" --generate-code arch=compute_12,code=sm_12 --ptxas-options="-v"
#CFLAGS=-Xcompiler "-Wall -O3 $(CFLAGS_PNG)" --generate-code arch=compute_20,code=sm_21 --ptxas-options="-v"
LDFLAGS=-Xlinker "-lm $(LDFLAGS_PNG)"
#LDFLAGS=-Xlinker "-lm -fopenmp"
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
