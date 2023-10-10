main:
	mpicxx smart-array.cpp -o smart-array.o

N ?= 8

run:
	mpirun -n $(N) ./smart-array.o

clean:
	rm *.o