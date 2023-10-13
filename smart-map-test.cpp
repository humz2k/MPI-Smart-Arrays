#include "mpi-smart-map.hpp"

class map_1{
    private:
        MPI_Comm comm;
        int comm_rank;
        int comm_size;
        int n;
    
    public:
        inline map_1(MPI_Comm comm_, int n_) : comm(comm_), n(n_){
            MPI_Comm_rank(comm,&comm_rank);
            MPI_Comm_size(comm,&comm_size);
        }

        inline ~map_1(){

        }

        inline map_return_t get(int idx){
            int out_rank = (comm_rank + 1)%comm_size;
            int out_idx = (idx + (idx/5)*2)%n;
            return make_map_return(out_rank,out_idx);
        }

};

void test(){
    int n = 10;
    
    map_1 map(MPI_COMM_WORLD,n);
    
    int* arr = (int*)malloc(sizeof(int) * n);

    SmartMap<int,map_1> my_smart_map(MPI_COMM_WORLD,arr,n,map);

    free(arr);
}

int main(){
    MPI_Init(NULL,NULL);
    
    test();

    MPI_Finalize();
    return 0;
};