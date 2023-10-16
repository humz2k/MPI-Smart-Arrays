main:
	mpicxx -O3 test.cpp -o test.o -fopenmp

N ?= 2

run:
	mpirun -n $(N) ./test.o

clean:
	rm *.o