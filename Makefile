# Starflood

CC = gcc

CFLAGS = -fopenmp -ggdb -march=x86-64 -mtune=generic -Og -pedantic -std=c11 -Wall -Wextra -Wshadow

LDFLAGS = -lm

.PHONY: all
all: starflood

starflood: src/main.o src/simulation.o src/visualization.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ src/main.c src/simulation.c src/visualization.c

src/main.o:
	$(CC) -c $(CFLAGS) -o $@ src/main.c

src/simulation.o:
	$(CC) -c $(CFLAGS) -o $@ src/simulation.c

src/visualization.o:
	$(CC) -c $(CFLAGS) -o $@ src/visualization.c

.PHONY: clean
clean:
	rm -fv src/main.o src/simulation.o src/visualization.o starflood
