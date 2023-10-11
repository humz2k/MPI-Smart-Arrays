#ifndef _SMARR_ARRAY_
#define _SMARR_ARRAY_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>
#include "hvector.hpp"
#include <cassert>
#include <math.h>

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
        size_t bytes_transferred;
        size_t bytes_used;
        double start;
        double end;
        bool free_raw;
    
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
        
        inline int get_cache_block_size(){
            return cache_block_size;
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
                                                                                                            hits(0),misses(0),evictions(0),
                                                                                                            bytes_transferred(0),bytes_used(0),free_raw(true){
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

        inline SmartArray(MPI_Comm comm_, T* raw_, int local_n_, int cache_block_size_ = 5, int max_cache_blocks_ = 10) : comm(comm_), local_n(local_n_),
                                                                                                            cache_block_size(cache_block_size_),
                                                                                                            max_cache_blocks(max_cache_blocks_),
                                                                                                            hits(0),misses(0),evictions(0),
                                                                                                            bytes_transferred(0),bytes_used(0),free_raw(false),raw(raw_){
            MPI_Comm_size(comm,&comm_size);
            MPI_Comm_rank(comm,&comm_rank);
            //local_n = (global_n + (comm_size-1)) / comm_size;
            global_n = local_n * comm_size;
            active = !comm_rank;
            n_cache_blocks = (global_n + (cache_block_size - 1) / cache_block_size) < max_cache_blocks ? (global_n + (cache_block_size - 1) / cache_block_size) : max_cache_blocks;
            cache_n = n_cache_blocks * cache_block_size;
            assert(local_n % cache_block_size == 0);
            //raw = (T*)malloc(sizeof(T) * local_n);
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
            if(free_raw)
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
            int cache_block = (tag % n_cache_blocks + rank) % n_cache_blocks;
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
            if (is_local(addr)){
                hits++;
                return raw[addr.rank_idx];
            }
            if (in_cache(addr)){
                hits++;
                return cache[addr.cache_idx];
            }
            misses++;
            if (valid[addr.cache_block])evictions++;

            MPI_Request req;
            MPI_Rget(&cache[addr.cache_block*cache_block_size],cache_block_size * sizeof(T),MPI_BYTE,addr.rank,addr.rank_block_start * sizeof(T),cache_block_size * sizeof(T),MPI_BYTE,win,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);
            valid[addr.cache_block] = true;
            tag[addr.cache_block] = addr.tag;
            return cache[addr.cache_idx];
        }

        inline void read_block(int idx, T* out){
            cache_addr_t addr = get_map(idx);
            T* read_from;
            if (is_local(addr)){
                hits++;
                read_from = &raw[addr.rank_block_start];
            } 
            else {
                if (in_cache(addr)){
                    hits++;
                    read_from = &cache[addr.cache_block*cache_block_size];
                } else {
                    misses++;
                    if(valid[addr.cache_block])evictions++;
                    MPI_Request req;
                    MPI_Rget(&cache[addr.cache_block*cache_block_size],cache_block_size * sizeof(T),MPI_BYTE,addr.rank,addr.rank_block_start * sizeof(T),cache_block_size * sizeof(T),MPI_BYTE,win,&req);
                    MPI_Wait(&req,MPI_STATUS_IGNORE);
                    valid[addr.cache_block] = true;
                    tag[addr.cache_block] = addr.tag;
                    read_from = &cache[addr.cache_block*cache_block_size];
                }
            }
            for (int i = 0; i < cache_block_size; i++){
                out[i] = read_from[i];
            }
        }

        inline void start_timer(){
            hits = 0;
            misses = 0;
            start = MPI_Wtime();
        }

        inline void stop_timer(){
            end = MPI_Wtime();
        }

        inline void print_stats(){
            bytes_used = (hits + misses) * sizeof(T);
            bytes_transferred = misses * sizeof(T) * cache_block_size;
            MPI_Barrier(comm);
            for (int i = 0; i < comm_size; i++){
                if (comm_rank == i){
                    printf("rank %d:\n   time = %g\n   hits = %d\n   misses = %d\n   evictions = %d\n   bytes_transferred = %lu\n   bytes_used = %lu\n   bytes percent used = %g\n",comm_rank,end-start,hits,misses,evictions,bytes_transferred,bytes_used,(((double)bytes_used))/((double)bytes_transferred)*100);
                }
                MPI_Barrier(comm);
            }
            MPI_Barrier(comm);
        }
};

#endif