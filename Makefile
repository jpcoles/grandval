CC=nvcc
#CUDA_VERSION=30
CUDA_VERSION=12
LDFLAGS_PNG=$(shell libpng-config --ldflags)
CFLAGS_PNG=$(shell libpng-config --cflags)
CFLAGS=-m64 -g -G -Xcompiler "-Wall -O3 $(CFLAGS_PNG)" --generate-code arch=compute_$(CUDA_VERSION),code=sm_$(CUDA_VERSION) --ptxas-options="-v"
LDFLAGS=-m64 -Xlinker "-lm $(LDFLAGS_PNG)"
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
