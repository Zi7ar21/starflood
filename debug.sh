g++ -O0 -ggdb3 -pedantic -std=c++11 -Wall -Wextra -Wshadow -c src/image.cpp -o build/image.o
g++ -O0 -ggdb3 -pedantic -std=c++11 -Wall -Wextra -Wshadow -c src/main.cpp -o build/main.o
g++ -O0 -ggdb3 -pedantic -std=c++11 -Wall -Wextra -Wshadow -o build/starflood.out build/image.o build/main.o
