CC=g++

all: build/mm.o build/main.o
	$(CC) build/*.o -o WIsH

build/mm.o: mm.cpp
	$(CC) -c mm.cpp -o build/mm.o -std=c++11

build/main.o: main.cpp
	$(CC) -c main.cpp -o build/main.o -std=c++11

clean:
	rm -f *.o WIsH
