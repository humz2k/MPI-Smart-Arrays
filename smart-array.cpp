#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>
#include "hvector.hpp"
#include <cassert>
#include <math.h>
#include <vector>
#include <algorithm>

#define SMARR_VERBOSE

typedef struct {
    int tag;
    int cache_block;
    int offset;
    int cache_idx;
    int rank;
    int rank_idx;
    int rank_block_start;
} cache_addr_t;

inline cache_addr_t make_cache_addr_t(int tag, int cache_block, int offset, int cache_idx, int rank, int rank_idx, int rank_block_start){
    cache_addr_t out;
    out.tag = tag;
    out.cache_block = cache_block;
    out.offset = offset;
    out.cache_idx = cache_idx;
    out.rank = rank;
    out.rank_idx = rank_idx;
    out.rank_block_start = rank_block_start;
    return out;
}

template<class T>
class SmartArray{
    private:
        T* raw;
        T* cache;
        bool* valid;
        int* tag;
        int comm_size;
        int comm_rank;
        MPI_Comm comm;
        MPI_Win win;
        int global_n;
        int local_n;
        int cache_block_size;
        int max_cache_blocks;
        int n_cache_blocks;
        int cache_n;
        bool active;
        int hits;
        int misses;
        int evictions;
        int rerecvs;
        size_t bytes_transferred;
        size_t bytes_used;
        std::vector<int> re_recv_tag;
    
    public:
        inline T* get_raw(){
            return raw;
        }

        inline int get_comm_rank(){
            return comm_rank;
        }

        inline int get_comm_size(){
            return comm_size;
        }

        inline int get_local_n(){
            return local_n;
        }

        inline int get_global_n(){
            return global_n;
        }

        inline void log(const char* format, ...)
        {
            if (active){
                va_list argptr;
                va_start(argptr, format);
                vfprintf(stdout, format, argptr);
                va_end(argptr);
            }
        }

        inline void lock_writes(){
            MPI_Win_fence(0,win);
        }

        inline SmartArray(MPI_Comm comm_, int global_n_, int cache_block_size_ = 5, int max_cache_blocks_ = 10) : comm(comm_), global_n(global_n_),
                                                                                                            cache_block_size(cache_block_size_),
                                                                                                            max_cache_blocks(max_cache_blocks_),
                                                                                                            hits(0),misses(0),evictions(0),rerecvs(0),
                                                                                                            bytes_transferred(0),bytes_used(0){
            MPI_Comm_size(comm,&comm_size);
            MPI_Comm_rank(comm,&comm_rank);
            local_n = (global_n + (comm_size-1)) / comm_size;
            active = !comm_rank;
            n_cache_blocks = (global_n + (cache_block_size - 1) / cache_block_size) < max_cache_blocks ? (global_n + (cache_block_size - 1) / cache_block_size) : max_cache_blocks;
            cache_n = n_cache_blocks * cache_block_size;
            assert(local_n % cache_block_size == 0);
            raw = (T*)malloc(sizeof(T) * local_n);
            cache = (T*)malloc(sizeof(T) * cache_n);
            valid = (bool*)malloc(sizeof(bool) * n_cache_blocks);
            tag = (int*)malloc(sizeof(int) * n_cache_blocks);

            MPI_Win_create(raw,sizeof(T)*local_n,1,MPI_INFO_NULL,comm,&win);

            for (int i = 0; i < n_cache_blocks; i++){
                valid[i] = false;
            }

            MPI_Win_fence(0,win);

            #ifdef SMARR_VERBOSE
            log("Initialed with:\n   comm_size = %d\n   global_n = %d\n   local_n = %d\n   cache_block_size = %d\n   max_cache_blocks = %d\n   n_cache_blocks = %d\n   cache_n = %d\n   cache \% = %g\%\n",comm_size,global_n,local_n,cache_block_size,max_cache_blocks,n_cache_blocks,cache_n,(((double)cache_n)/(double)local_n)*100);//,(((double)(cache_n * comm_size))/(double)global_n)*100);
            #endif

        }

        inline ~SmartArray(){
            #ifdef SMARR_VERBOSE
            log("Finalizing smart array!\n");
            #endif
            MPI_Win_fence(0,win);
            MPI_Win_free(&win);
            free(raw);
            free(cache);
            free(valid);
            free(tag);
        }

        inline cache_addr_t get_map(int idx){
            int rank = idx / local_n;
            int rank_idx = idx % local_n;
            int rank_block_start = (rank_idx / cache_block_size) * cache_block_size;

            int tag = idx / cache_block_size;
            int cache_block = (tag % n_cache_blocks);// + rank) % n_cache_blocks;
            int offset = idx % cache_block_size;
            int cache_idx = cache_block * cache_block_size + offset;
            
            return make_cache_addr_t(tag,cache_block,offset,cache_idx,rank,rank_idx,rank_block_start);
        }

        inline bool in_cache(cache_addr_t addr){
            if (!(valid[addr.cache_block]))return false;
            if (tag[addr.cache_block] != addr.tag)return false;
            return true;
        }

        inline bool is_local(cache_addr_t addr){
            return (addr.rank == comm_rank);
        }

        inline T read(int idx){
            cache_addr_t addr = get_map(idx);
            if (in_cache(addr)){
                hits++;
                #ifdef SMARR_VERBOSE
                //printf("rank %d has idx %d in cache (cache_idx = %d)!\n",comm_rank,idx,addr.cache_idx);
                #endif
                bytes_used += sizeof(T);
                return cache[addr.cache_idx];
            }
            if (is_local(addr)){
                hits++;
                bytes_used += sizeof(T);
                return raw[addr.rank_idx];
            }
            misses++;
            if (valid[addr.cache_block])evictions++;
            bytes_transferred += sizeof(T) * cache_block_size;
            bytes_used += sizeof(T);
            if(std::find(re_recv_tag.begin(), re_recv_tag.end(), addr.tag) != re_recv_tag.end()) {
                if(comm_rank == 0)
                    printf("re-recv tag %d\n",addr.tag);
                rerecvs++;
            } else {
                if(comm_rank == 0)
                    printf("recv tag %d\n",addr.tag);
                re_recv_tag.push_back(addr.tag);
            }

            #ifdef SMARR_VERBOSE
            //printf("rank %d does not have idx %d in cache!\n",comm_rank,idx);
            #endif
            MPI_Request req;
            MPI_Rget(&cache[addr.cache_block*cache_block_size],cache_block_size * sizeof(T),MPI_BYTE,addr.rank,addr.rank_block_start * sizeof(T),cache_block_size * sizeof(T),MPI_BYTE,win,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);
            valid[addr.cache_block] = true;
            tag[addr.cache_block] = addr.tag;
            return cache[addr.cache_idx];
        }

        /*inline void write(int idx, T data){
            cache_addr_t addr = get_map(idx);
            #ifdef SMARR_VERBOSE
            //printf("rank %d writing to idx %d!\n",comm_rank,idx);
            #endif
            MPI_Put(&data,sizeof(T),MPI_BYTE,addr.rank,addr.rank_idx*sizeof(T),sizeof(T),MPI_BYTE,win);
        }*/

        inline void print_stats(){
            MPI_Barrier(comm);
            for (int i = 0; i < comm_size; i++){
                if (comm_rank == i){
                    printf("rank %d:\n   hits = %d\n   misses = %d\n   evictions = %d\n   rerecvs = %d\n   bytes_transferred = %lu\n   bytes_used = %lu\n   bytes percent used = %g\n",comm_rank,hits,misses,evictions,rerecvs,bytes_transferred,bytes_used,(((double)bytes_used))/((double)bytes_transferred)*100);
                }
                MPI_Barrier(comm);
            }
            MPI_Barrier(comm);
        }
};

typedef struct {
    int rank;
    int idx;
    int global_idx;
} map_t;

inline map_t map_1(int x, int y, int z, hvec3<int> local_grid_size, hvec3<int> dims){
    int coord_x = x / local_grid_size[0];
    int coord_y = y / local_grid_size[1];
    int coord_z = z / local_grid_size[2];
    int rank = coord_x * dims[1] * dims[2] + coord_y * dims[2] + coord_z;
    int local_x = x % local_grid_size[0];
    int local_y = y % local_grid_size[1];
    int local_z = z % local_grid_size[2];
    int local_idx = local_x * local_grid_size[1] * local_grid_size[2] + local_y * local_grid_size[2] + local_z;
    map_t out;
    out.rank = rank;
    out.idx = local_idx;
    int global_idx = rank * local_grid_size[0] * local_grid_size[1] * local_grid_size[2] + local_idx;
    out.global_idx = global_idx;
    return out;
}

void test(){
    hvec3<int> ng = make_hvec3(8,8,8);
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

    if(!world_rank){
        printf("Optimal memory: %d | Actual memory: %d | Diff percent: %g\n",ng.x*ng.y*ng.z,nlocal*world_size,(((double)(nlocal*world_size - ng.x*ng.y*ng.z))/((double)ng.x*ng.y*ng.z))*100);
    }

    SmartArray<int> arr(MPI_COMM_WORLD,nlocal * world_size,local_grid_size.z,world_size);

    MPI_Barrier(MPI_COMM_WORLD);

    int* raw = arr.get_raw();

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

    for (int i = my_pencil_start; i < my_pencil_end; i++){
        int x = i / (ng.y * ng.z);
        int y = (i - x*ng.y*ng.z) / ng.z;
        int z = (i - x*ng.y*ng.z) - y * ng.z;
        map_t out = map_1(x,y,z,local_grid_size,dims);
        local_arr[i - my_pencil_start] = arr.read(out.global_idx);
        //printf("rank %d getting %d %d %d\n",world_rank,x,y,z);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank == 0){
        for(int i = 0; i < local_sz; i++){
            printf("%d = %d\n",i,local_arr[i]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(local_arr);

    //int n_pencils = ng.x*ng.y;
    //int n_per = (n_pencils + (world_size - 1)) / world_size;
    /*int n_x_per = (ng.x + (world_size - 1))/world_size;

    int local_pencil_size = n_per*ng.z;
    int* local_pencils = (int*)malloc(sizeof(int)*local_pencil_size);

    int x_start = n_x_per * world_rank;
    int x_end = x_start + n_x_per;
    idx = 0;
    for (int x = x_start; x < x_end; x++){
        for (int y = 0; y < ng.y; y++){
            for (int z = 0; z < ng.z; z++){
                map_t out = map_1(x,y,z,local_grid_size,dims);
                if(world_rank == 0){
                    printf("recv from rank %d idx %d, global %d\n",out.rank,out.idx,out.global_idx);
                }
                local_pencils[idx] = arr.read(out.global_idx);
                idx++;
                //printf("rank %d needs to get idx %d\n",world_rank,out.global_idx);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank == 1){
        for(int i = 0; i < local_pencil_size; i++){
            printf("%d = %d\n",i,raw[i]);
        }
    }*/
    
    //free(local_pencils);

    arr.print_stats();
}

int main(){
    MPI_Init(NULL,NULL);
    test();
    MPI_Finalize();
    return 0;
}