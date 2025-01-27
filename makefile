.PHONY: main
main:
	mpic++ -O3 -std=c++17 test.cpp -o test.o -fopenmp

clean:
	rm *.o