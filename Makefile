# Starflood

CC ?= gcc

CFLAGS = -fopenmp -ggdb -Og -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow

LDFLAGS = -lm

.PHONY: all
all: build/starflood

build/starflood: build/main.o build/rng.o build/simulation.o build/solver.o build/visualization.o
	$(CC) $(CFLAGS) build/*.o -o $@ $(LDFLAGS)

build/main.o:
	$(CC) $(CFLAGS) -c -o $@ src/main.c

build/rng.o:
	$(CC) $(CFLAGS) -c -o $@ src/rng.c

build/simulation.o:
	$(CC) $(CFLAGS) -c -o $@ src/simulation.c

build/solver.o:
	$(CC) $(CFLAGS) -c -o $@ src/solver.c

build/visualization.o:
	$(CC) $(CFLAGS) -c -o $@ src/visualization.c

.PHONY: clean
clean:
	rm -fv build/*.o build/starflood
