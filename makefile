main:
	mpicxx -O3 fft_test.cpp -o fft_test.o -fopenmp

N ?= 8

run:
	mpirun -n $(N) ./fft_test.o

clean:
	rm *.o