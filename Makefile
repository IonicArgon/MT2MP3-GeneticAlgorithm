CC=gcc
CFLAGS=-lm -Ofast
DFLAGS=-Wall -Wextra -W -pedantic -O0 -g -pg -lm
ELITISM=

DEPS = functions.h
OBJS = GA.o OF.o functions.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(ELITISM)

GA: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(ELITISM)

.PHONY: clean

clean:
	rm -f *.o GA