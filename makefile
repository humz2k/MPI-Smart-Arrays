.PHONY: main
main:
	mpicxx -O3 -std=c++20 test.cpp -o test.o -fopenmp

clean:
	rm *.o