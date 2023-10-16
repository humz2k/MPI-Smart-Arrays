#ifndef _FFT_MAP_
#define _FFT_MAP_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>
#include "hvector.hpp"
#include "mpi-map.hpp"
#include <cassert>
#include <math.h>

class map_1{
    private:
        MPI_Comm comm;
        hvec3<int> ng;
        hvec3<int> local_grid_size;
        hvec3<int> dims;
        hvec3<int> coords;
        int world_rank;
        int world_size;
        int nlocal;
        int total_pencils;
        int pencils_per_rank;
        int total_per_rank;
    
    public:
        inline map_1(MPI_Comm comm_, hvec3<int> ng_, hvec3<int> local_grid_size_, hvec3<int> dims_, hvec3<int> coords_, int nlocal_) : comm(comm_), ng(ng_), local_grid_size(local_grid_size_), dims(dims_), coords(coords_), nlocal(nlocal_){
            MPI_Comm_size(comm,&world_size);
            MPI_Comm_rank(comm,&world_rank); 
            total_pencils = ng.x * ng.y;
            pencils_per_rank = total_pencils / world_size;
            total_per_rank = pencils_per_rank * ng.z;

        }

        inline ~map_1(){

        }

        inline map_return_t map(int i){
            
            int this_pencil = i / ng.z;
            int this_pencil_start = world_rank * pencils_per_rank + this_pencil;
            int pencil_x = this_pencil_start / ng.y;
            int pencil_y = this_pencil_start % ng.y;
            int pencil_z = i%ng.z;

            int idx = pencil_x * ng.y * ng.z + pencil_y * ng.z + pencil_z;
            
            int x = idx / (ng.y * ng.z);
            int y = (idx/ng.z)%ng.y;
            int z = idx%ng.z;

            int coord_x = x / local_grid_size.x;
            int coord_y = y / local_grid_size.y;
            int coord_z = z / local_grid_size.z;
            int rank = coord_x * dims.y * dims.z + coord_y * dims.z + coord_z;
            int local_x = x % local_grid_size.x;
            int local_y = y % local_grid_size.y;
            int local_z = z % local_grid_size.z;
            int local_idx = local_x * local_grid_size.y * local_grid_size.z + local_y * local_grid_size.z + local_z;

            //if(!world_rank)
            //    printf("rank %d mapping %d to rank %d idx %d\n",world_rank,i,rank,local_idx);

            return make_map_return(rank,local_idx,i);
        }

};

class map_2{
    private:
        MPI_Comm comm;
        hvec3<int> ng;
        hvec3<int> local_grid_size;
        hvec3<int> dims;
        hvec3<int> coords;
        int world_rank;
        int world_size;
        int nlocal;
        int total_pencils;
        int pencils_per_rank;
        int total_per_rank;
    
    public:
        inline map_2(MPI_Comm comm_, hvec3<int> ng_, hvec3<int> local_grid_size_, hvec3<int> dims_, hvec3<int> coords_, int nlocal_) : comm(comm_), ng(ng_), local_grid_size(local_grid_size_), dims(dims_), coords(coords_), nlocal(nlocal_){
            MPI_Comm_size(comm,&world_size);
            MPI_Comm_rank(comm,&world_rank);
            total_pencils = ng.x * ng.z;
            pencils_per_rank = total_pencils / world_size;
            total_per_rank = pencils_per_rank * ng.y;

        }

        inline ~map_2(){

        }

        inline map_return_t map(int i){
            
            int this_pencil = i / ng.y;
            int this_pencil_start = world_rank * pencils_per_rank + this_pencil;
            int pencil_x = this_pencil_start / ng.z;
            int pencil_z = this_pencil_start % ng.z;
            int pencil_y = i%ng.y;

            int idx = pencil_x * ng.y * ng.z + pencil_y * ng.z + pencil_z;

            int pencil_id = pencil_x * ng.y + pencil_y;
            int rank = pencil_id / ((ng.x*ng.y)/world_size);
            int local_idx = idx % (((ng.x*ng.y)/world_size) * ng.z);

            return make_map_return(rank,local_idx,i);
            
        }

};

class map_3{
    private:
        MPI_Comm comm;
        hvec3<int> ng;
        hvec3<int> local_grid_size;
        hvec3<int> dims;
        hvec3<int> coords;
        int world_rank;
        int world_size;
        int nlocal;
        int total_pencils;
        int pencils_per_rank;
        int total_per_rank;
    
    public:
        inline map_3(MPI_Comm comm_, hvec3<int> ng_, hvec3<int> local_grid_size_, hvec3<int> dims_, hvec3<int> coords_, int nlocal_) : comm(comm_), ng(ng_), local_grid_size(local_grid_size_), dims(dims_), coords(coords_), nlocal(nlocal_){
            MPI_Comm_size(comm,&world_size);
            MPI_Comm_rank(comm,&world_rank);
            total_pencils = ng.z * ng.y;
            pencils_per_rank = total_pencils / world_size;
            total_per_rank = pencils_per_rank * ng.x;

        }

        inline ~map_3(){

        }

        inline map_return_t map(int i){
            
            int this_pencil = i / ng.x;
            int this_pencil_start = world_rank * pencils_per_rank + this_pencil;
            int pencil_z = this_pencil_start / ng.y;
            int pencil_y = this_pencil_start % ng.y;
            int pencil_x = i%ng.x;

            int idx = pencil_x * ng.y * ng.z + pencil_y * ng.z + pencil_z;

            int pencil_id = pencil_x * ng.z + pencil_z;
            int rank = pencil_id / ((ng.x*ng.z)/world_size);
            int local_idx = idx % (((ng.x*ng.z)/world_size) * ng.y);

            return make_map_return(rank,local_idx,i);
            
        }

};


#endif