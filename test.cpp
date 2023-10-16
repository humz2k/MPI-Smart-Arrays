#include <mpi.h>
#include "mpi-map.hpp"

class Map_1{
    public:
        MPI_Comm comm;
        int rank;
        int size;
        int n;
    
    Map_1(MPI_Comm comm_, int n_) : comm(comm_), n(n_){
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&size);
    };
    ~Map_1(){};

    map_return_t map(int i){
        int n_per_rank = n/size;
        int rank_to_get = i / n_per_rank;
        return make_map_return((rank_to_get)%size,(i),i);
    }
};

int main(){
    MPI_Init(NULL,NULL);

    int n = 8;
    int world_rank;
    int world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    
    Map_1 map(MPI_COMM_WORLD,n);

    double start = MPI_Wtime();
    SmartMap<Map_1> my_smart_map(MPI_COMM_WORLD,map,n);
    double end = MPI_Wtime();
    printf("setup took %g for rank %d\n",end-start,world_rank);

    int* in = (int*)malloc(sizeof(int) * n);
    int* out = (int*)malloc(sizeof(int) * n);

    for (int i = 0; i < n; i++){
        in[i] = world_rank;
    }

    my_smart_map.forward(in,out);

    if(world_rank == 1){
        for (int i = 0; i < n; i++){
            printf("rank %d idx %d = %d\n",world_rank,i,out[i]);
        }
    }
    
    free(in);
    free(out);
    MPI_Finalize();
    return 0;
}