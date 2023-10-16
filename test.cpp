#include <mpi.h>
#include "gets.hpp"

class Map_1{
    public:
        MPI_Comm comm;
        int rank;
        int size;
    
    Map_1(MPI_Comm comm_) : comm(comm_){
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&size);
    };
    ~Map_1(){};

    map_return_t map(int i){
        int n_per_rank = 256*256*256/8;
        int rank_to_get = i / n_per_rank;
        return make_map_return((rank_to_get)%size,(i),i);
    }
};

int main(){
    MPI_Init(NULL,NULL);

    Map_1 map(MPI_COMM_WORLD);

    get_t* init = find_gets(256*256*256,map);

    if (!map.rank){
        get_t* cur = init;
        while (cur){
            printf("rank %d: %d items at idx %d mapped from rank %d idx %d with stride %d\n",map.rank,cur->n,cur->dest,cur->rank,cur->src,cur->stride);
            cur = cur->next;
        }
    }

    free_get(init);

    MPI_Finalize();
    return 0;
}