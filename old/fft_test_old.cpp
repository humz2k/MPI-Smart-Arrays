#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>
#include "hvector.hpp"
#include "smart-array.hpp"
#include <cassert>
#include <math.h>

//#define READ_BLOCKS
#define MAP_CACHE_SZ 4
int map_cache[MAP_CACHE_SZ];
int map_cache_ptr = 0;

inline int map_1(int i, hvec3<int>& ng, hvec3<int>& local_grid_size, hvec3<int>& dims, int nlocal){

    int x = i / (ng.y * ng.z);
    int y = (i/ng.z)%ng.y;
    int z = i%ng.z;

    int coord_x = x / local_grid_size.x;
    int coord_y = y / local_grid_size.y;
    int coord_z = z / local_grid_size.z;
    int rank = coord_x * dims.y * dims.z + coord_y * dims.z + coord_z;
    int local_x = x % local_grid_size.x;
    int local_y = y % local_grid_size.y;
    int local_z = z % local_grid_size.z;
    int local_idx = local_x * local_grid_size.y * local_grid_size.z + local_y * local_grid_size.z + local_z;
    int global_idx = rank * nlocal + local_idx;

    return global_idx;
}

void test(){
    hvec3<int> ng = make_hvec3(256,256,256);
    int world_rank;
    int world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    int _dims[3] = {0,0,0};
    MPI_Dims_create(world_size,3,_dims);
    hvec3<int> dims = make_hvec3(_dims[0],_dims[1],_dims[2]);
    hvec3<int> local_grid_size = (ng + (dims - 1)) / dims;

    hvec3<int> coords;
    coords[0] = world_rank / (dims[1] * dims[2]);
    coords[1] = (world_rank - coords[0]*dims[1]*dims[2]) / dims[2];
    coords[2] = (world_rank - coords[0]*dims[1]*dims[2]) - coords[1] * dims[2];

    //printf("rank %d coords = %d %d %d, local_grid = %d %d %d\n",world_rank,coords[0],coords[1],coords[2],local_grid_size[0],local_grid_size[1],local_grid_size[2]);
    if(!world_rank){
        printf("dims = [%d %d %d]\nlocal_grid_size = [%d %d %d]\nng = [%d %d %d]\n",dims.x,dims.y,dims.z,local_grid_size.x,local_grid_size.y,local_grid_size.z,dims.x,dims.y,dims.z);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);


    int nlocal = local_grid_size.x*local_grid_size.y*local_grid_size.z;

    assert(nlocal % world_size == 0);

    if(!world_rank){
        printf("Optimal memory: %d | Actual memory: %d | Diff percent: %g\n",ng.x*ng.y*ng.z,nlocal*world_size,(((double)(nlocal*world_size - ng.x*ng.y*ng.z))/((double)ng.x*ng.y*ng.z))*100);
    }

    int* raw = (int*)malloc(sizeof(int)*nlocal);

    SmartArray<int> arr(MPI_COMM_WORLD,raw,nlocal,local_grid_size.z,world_size*4);

    MPI_Barrier(MPI_COMM_WORLD);

    int idx = 0;
    for (int i = 0; i < local_grid_size[0]; i++){
        for (int j = 0; j < local_grid_size[1]; j++){
            for (int k = 0; k < local_grid_size[2]; k++){
                int x = i + local_grid_size[0] * coords[0];
                int y = j + local_grid_size[1] * coords[1];
                int z = k + local_grid_size[2] * coords[2];
                int global_idx = x * ng.y * ng.z + y * ng.z + z;
                raw[idx] = global_idx;
                //if(world_rank == 1){
                //    printf("%d = %d\n",idx,global_idx);
                //}
                idx++;
            }
        }
    }
    arr.lock_writes();
    MPI_Barrier(MPI_COMM_WORLD);

    int n_pencils = ng.x * ng.y;
    int pencils_per_rank = (n_pencils + (world_size - 1)) / world_size;
    int my_pencil_start = pencils_per_rank * ng.z * world_rank;
    int my_pencil_end = pencils_per_rank * ng.z + my_pencil_start;
    if (my_pencil_end > ng.x*ng.y*ng.z){
        printf("!!!\n");
        my_pencil_end = ng.x*ng.y*ng.z;
    }

    printf("rank %d start = %d end = %d\n",world_rank,my_pencil_start,my_pencil_end);

    int nffts = (my_pencil_end - my_pencil_start) / ng.z;
    assert(nffts * ng.z == my_pencil_end - my_pencil_start);
    int local_sz = nffts*ng.z;
    int* local_arr = (int*)malloc(sizeof(int)*local_sz);
    arr.start_timer();
    //map_1 map;

    #ifdef READ_BLOCKS
    for (int i = my_pencil_start; i < my_pencil_end; i += arr.get_cache_block_size()){
        //int x = i / (ng.y * ng.z);
        //int y = (i - x*ng.y*ng.z) / ng.z;
        //int z = (i - x*ng.y*ng.z) - y * ng.z;
        int out = map_1(i,ng,local_grid_size,dims,nlocal);
        arr.read_block(out,&local_arr[i-my_pencil_start]);
        //local_arr[i - my_pencil_start] = arr.read(out);
        //printf("rank %d getting %d %d %d\n",world_rank,x,y,z);
    }
    #else
    //#pragma omp parallel for
    for (int i = my_pencil_start; i < my_pencil_end; i++){
        //int x = i / (ng.y * ng.z);
        //int y = (i/ng.z)%ng.y;
        //int z = i%ng.z;
        
        int my_map = map_1(i,ng,local_grid_size,dims,nlocal);
        //#pragma omp critical
        int out = arr.read(my_map);

        local_arr[i - my_pencil_start] = out;
    }
    #endif

    arr.stop_timer();

    MPI_Barrier(MPI_COMM_WORLD);

    /*if(world_rank == 0){
        for(int i = 0; i < local_sz; i++){
            printf("%d = %d\n",i,local_arr[i]);
        }
    }*/

    //MPI_Barrier(MPI_COMM_WORLD);

    free(local_arr);

    free(raw);

    arr.print_stats();
}

int main(){
    MPI_Init(NULL,NULL);
    test();
    MPI_Finalize();
    return 0;
}