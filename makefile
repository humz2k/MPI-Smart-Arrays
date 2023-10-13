main:
	mpicxx -O3 smart-map-test.cpp -o smart-map-test.o -fopenmp

N ?= 2

run:
	mpirun -n $(N) ./smart-map-test.o

clean:
	rm *.o