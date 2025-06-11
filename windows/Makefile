CXX=-i686-w64-mingw32-g++
CFLAGS=-c -Wall -std=c++20 -g

all: read_binary

read_binary: main.o read_binary.o model.o
	$(CXX) main.o read_binary.o model.o -o read_binary -static

main.o: main.cpp read_binary.hpp
	$(CXX) $(CFLAGS) main.cpp -static

read_binary.o: read_binary.cpp read_binary.hpp
	$(CXX) $(CFLAGS) read_binary.cpp -static

model.o: model.cpp read_binary.hpp
	$(CXX) $(CFLAGS) model.cpp -static

clean:
	rm -rf *.o 
