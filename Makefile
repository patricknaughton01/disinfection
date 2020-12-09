CC=g++
CFLAGS=-Wall -I. -std=c++11 -g
LIBS=-L/usr/local/lib/ -lKrisLibrary -llog4cxx -laprutil-1 \
	-lapr-1

all: helper.o plane.o plane_finder.o merge_triangle_mesh.o main.o
	$(CC) $(CFLAGS) -o finder $^ $(LIBS)

build:
	g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` \
	merge_triangle_mesh.cpp plane_finder.cpp plane.cpp helper.cpp pybind.cpp \
	-o merge_triangle_mesh`python3-config --extension-suffix` $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $^

.PHONY:clean
clean:
	rm -f *.o finder
