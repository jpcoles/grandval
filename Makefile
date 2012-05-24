CC=nvcc
CFLAGS=#-Wall -O3 -fopenmp
LDFLAGS=#-lm -fopenmp
SRC=$(wildcard *.c*)
#OBJS=$(SRC:.c*=.o)

all: $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o grandval

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

%.o: %.cu %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm $(OBJS)
