CC=nvcc
HOSTNAME=$(shell hostname)
ifeq ($(HOSTNAME),mayall.physik.uzh.ch)
	CUDA_VERSION=30
	ARCH=64
else
	CUDA_VERSION=12
	ARCH=32
endif
LDFLAGS_PNG=$(shell libpng-config --ldflags)
CFLAGS_PNG=$(shell libpng-config --cflags)
CFLAGS=-m$(ARCH) -g -G -Xcompiler "-Wall -O3 $(CFLAGS_PNG)" --generate-code arch=compute_$(CUDA_VERSION),code=sm_$(CUDA_VERSION) --ptxas-options="-v"
LDFLAGS=-m$(ARCH) -Xlinker "-lm $(LDFLAGS_PNG)"
SRC=$(wildcard *.cu)
OBJS=$(patsubst %.cu,%.o, $(SRC))

all: $(OBJS)
	echo $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o grandval

%.o: %.cu
	echo $(OBJS) $(SRC)
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o
