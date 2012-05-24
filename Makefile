CC=nvcc
CFLAGS=-Xcompiler "-Wall -O3 -fopenmp"
LDFLAGS=-Xlinker "-lm -fopenmp"
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
