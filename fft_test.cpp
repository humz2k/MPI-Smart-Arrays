#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>
#include "hvector.hpp"
#include "mpi-map.hpp"
#include <cassert>
#include <math.h>
#include "fft_maps.hpp"

struct cell{
    int x;
    int y;
    int z;
};

struct double2{
    double x;
    double y;
};

inline cell make_cell(int x, int y, int z){
    cell out;
    out.x = x;
    out.y = y;
    out.z = z;
    return out;
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

    if(!world_rank){
        printf("dims = [%d %d %d]\nlocal_grid_size = [%d %d %d]\nng = [%d %d %d]\n",dims.x,dims.y,dims.z,local_grid_size.x,local_grid_size.y,local_grid_size.z,dims.x,dims.y,dims.z);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    int nlocal = local_grid_size.x*local_grid_size.y*local_grid_size.z;

    assert(nlocal % world_size == 0);

    if(!world_rank){
        printf("Optimal memory: %d | Actual memory: %d | Diff percent: %g\n",ng.x*ng.y*ng.z,nlocal*world_size,(((double)(nlocal*world_size - ng.x*ng.y*ng.z))/((double)ng.x*ng.y*ng.z))*100);
    }

    cell* buff1 = (cell*)malloc(sizeof(cell)*nlocal);
    cell* buff2 = (cell*)malloc(sizeof(cell)*nlocal);

    MPI_Barrier(MPI_COMM_WORLD);

    int idx = 0;
    for (int i = 0; i < local_grid_size[0]; i++){
        for (int j = 0; j < local_grid_size[1]; j++){
            for (int k = 0; k < local_grid_size[2]; k++){
                int x = i + local_grid_size[0] * coords[0];
                int y = j + local_grid_size[1] * coords[1];
                int z = k + local_grid_size[2] * coords[2];
                int global_idx = x * ng.y * ng.z + y * ng.z + z;
                buff1[idx] = make_cell(x,y,z);
                //buff1[idx].y = global_idx;//make_cell(x,y,z);
                idx++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    map_1 map1(MPI_COMM_WORLD,ng,local_grid_size,dims,coords,nlocal);
    SmartMap<map_1> my_smart_map_1(MPI_COMM_WORLD,map1,nlocal);

    map_2 map2(MPI_COMM_WORLD,ng,local_grid_size,dims,coords,nlocal);
    SmartMap<map_2> my_smart_map_2(MPI_COMM_WORLD,map2,nlocal);

    map_3 map3(MPI_COMM_WORLD,ng,local_grid_size,dims,coords,nlocal);
    SmartMap<map_3> my_smart_map_3(MPI_COMM_WORLD,map3,nlocal);

    int print_rank = 0;

    double start = MPI_Wtime();

    my_smart_map_1.forward(buff1,buff2);
    MPI_Barrier(MPI_COMM_WORLD);
    /*if(world_rank == print_rank){
        printf("\n\n");
        for (int i = 0; i < nlocal; i++){
            printf("rank %d: buff2[%d] = [%d,%d,%d]\n",world_rank,i,buff2[i].x,buff2[i].y,buff2[i].z);
        }
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);
    my_smart_map_2.forward(buff2,buff1);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    /*if(world_rank == print_rank){
        printf("\n\n");
        for (int i = 0; i < nlocal; i++){
            printf("rank %d: buff1[%d] = [%d,%d,%d]\n",world_rank,i,buff1[i].x,buff1[i].y,buff1[i].z);
        }
    }*/
    //MPI_Barrier(MPI_COMM_WORLD);
    my_smart_map_3.forward(buff1,buff2);

    
    //MPI_Barrier(MPI_COMM_WORLD);
    
    //if(world_rank == print_rank){
    //    printf("\n\n");
    //    for (int i = 0; i < 100; i++){
    //        printf("rank %d: buff2[%d] = [%d,%d,%d]\n",world_rank,i,buff2[i].x,buff2[i].y,buff2[i].z);
    //    }
    //}

    double end = MPI_Wtime();

    printf("Time = %g\n",end-start);

    my_smart_map_3.backward(buff2,buff1);

    /*if(world_rank == print_rank){
        printf("\n\n");
        for (int i = 0; i < nlocal; i++){
            printf("rank %d: buff1[%d] = [%d,%d,%d]\n",world_rank,i,buff1[i].x,buff1[i].y,buff1[i].z);
        }
    }*/

    my_smart_map_2.backward(buff1,buff2);

    /*if(world_rank == print_rank){
        printf("\n\n");
        for (int i = 0; i < nlocal; i++){
            printf("rank %d: buff2[%d] = [%d,%d,%d]\n",world_rank,i,buff2[i].x,buff2[i].y,buff2[i].z);
        }
    }*/

    my_smart_map_1.backward(buff2,buff1);

    if(world_rank == print_rank){
        printf("\n\n");
        for (int i = 0; i < 200; i++){
            printf("rank %d: buff1[%d] = [%d,%d,%d]\n",world_rank,i,buff1[i].x,buff1[i].y,buff1[i].z);
        }
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    free(buff1);
    free(buff2);
}

int main(){
    MPI_Init(NULL,NULL);
    test();
    MPI_Finalize();
    return 0;
}