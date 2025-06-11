CXX=g++
CFLAGS=-c -Wall -std=c++20 -g

all: read_binary

read_binary: main.o read_binary.o model.o
	$(CXX) main.o read_binary.o model.o -o read_binary

main.o: main.cpp read_binary.hpp
	$(CXX) $(CFLAGS) main.cpp

read_binary.o: read_binary.cpp read_binary.hpp
	$(CXX) $(CFLAGS) read_binary.cpp

model.o: model.cpp read_binary.hpp
	$(CXX) $(CFLAGS) model.cpp

clean:
	rm -rf *.o 
