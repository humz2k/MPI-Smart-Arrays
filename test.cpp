#include "mpi-smart-arrays.hpp"

int main(){
    MPI_Init(NULL,NULL);

    smartarray::SmartArray<int> arr(100, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}