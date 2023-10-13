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
            int n_per = n / comm_size;

            int out_rank = idx / n_per;//(comm_rank + 1)%comm_size;
            int out_idx = idx % n_per;//(idx + (idx/5)*2)%n;
            return make_map_return(out_rank,out_idx);
        }

};

void test(){
    int n = 128;
    
    map_1 map(MPI_COMM_WORLD,n);
    
    int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    int* in = (int*)malloc(sizeof(int) * n);
    int* out = (int*)malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++){
        in[i] = i;
        if(!world_rank)printf("in[%d] = %d\n",i,in[i]);
    }

    SmartMap<map_1> my_smart_map(MPI_COMM_WORLD,n,map);

    my_smart_map.execute(in,out);

    for (int i = 0; i < n; i++){
        if(!world_rank)printf("out[%d] = %d\n",i,out[i]);
    }

    free(in);
    free(out);
}

int main(){
    MPI_Init(NULL,NULL);
    
    test();

    MPI_Finalize();
    return 0;
};