CC="gcc"

CFLAGS="-fopenmp -g -march=native -Og -pedantic -std=c99 -Wall -Wextra -Wshadow"

all: build/starflood
	$(CC) src/main.c -o

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)
