#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>

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
    
    public:
        inline T* get_raw(){
            return raw;
        }

        inline int get_comm_rank(){
            return comm_rank;
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

        inline SmartArray(MPI_Comm comm_, int global_n_, int cache_block_size_ = 5, int max_cache_blocks_ = 10) : comm(comm_), global_n(global_n_),
                                                                                                            cache_block_size(cache_block_size_),
                                                                                                            max_cache_blocks(max_cache_blocks_){
            MPI_Comm_size(comm,&comm_size);
            MPI_Comm_rank(comm,&comm_rank);
            local_n = (global_n + (comm_size-1)) / comm_size;
            active = !comm_rank;
            n_cache_blocks = (global_n + (cache_block_size - 1) / cache_block_size) < max_cache_blocks ? (global_n + (cache_block_size - 1) / cache_block_size) : max_cache_blocks;
            cache_n = n_cache_blocks * cache_block_size;
            
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
            log("Initialed with:\n   comm_size = %d\n   global_n = %d\n   local_n = %d\n   cache_block_size = %d\n   max_cache_blocks = %d\n   n_cache_blocks = %d\n   cache_n = %d\n",comm_size,global_n,local_n,cache_block_size,max_cache_blocks,n_cache_blocks,cache_n);
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
            int tag = idx / cache_block_size;
            int cache_block = tag % n_cache_blocks;
            int offset = idx % cache_block_size;
            int cache_idx = cache_block * cache_block_size + offset;
            int rank = idx / local_n;
            int rank_idx = idx % local_n;
            int rank_block_start = (rank_idx / cache_block_size) * cache_block_size;
            return make_cache_addr_t(tag,cache_block,offset,cache_idx,rank,rank_idx,rank_block_start);
        }

        inline bool in_cache(cache_addr_t addr){
            if (!(valid[addr.cache_block]))return false;
            if (tag[addr.cache_block] != addr.tag)return false;
            return true;
        }

        inline T read(int idx){
            cache_addr_t addr = get_map(idx);
            if (in_cache(addr)){
                #ifdef SMARR_VERBOSE
                printf("rank %d has idx %d in cache (cache_idx = %d)!\n",comm_rank,idx,addr.cache_idx);
                #endif
                return cache[addr.cache_idx];
            }
            #ifdef SMARR_VERBOSE
            printf("rank %d does not have idx %d in cache!\n",comm_rank,idx);
            #endif
            MPI_Request req;
            MPI_Rget(&cache[addr.cache_block*cache_block_size],cache_block_size * sizeof(T),MPI_BYTE,addr.rank,addr.rank_block_start * sizeof(T),cache_block_size * sizeof(T),MPI_BYTE,win,&req);
            MPI_Wait(&req,MPI_STATUS_IGNORE);
            valid[addr.cache_block] = true;
            tag[addr.cache_block] = addr.tag;
            return cache[addr.cache_idx];
        }
};


void test(){
    SmartArray<int> arr(MPI_COMM_WORLD,100);
    int* raw = arr.get_raw();
    int start = arr.get_comm_rank() * arr.get_local_n();
    for (int i = 0; i < arr.get_local_n(); i++){
        raw[i] = start + i;
    }

    if (!arr.get_comm_rank()){
        for (int i = 0; i < 100; i++){
            printf("%d: %d\n",i,arr.read(i));
        }
        /*for (int i = 0; i < 100; i++){
            cache_addr_t val = arr.get_map(i);
            printf("%d : tag = %d, block = %d, offset = %d, rank = %d, rank_idx = %d\n",i,val.tag,val.cache_block,val.offset, val.rank, val.rank_idx);
        }*/
    }
}

int main(){
    MPI_Init(NULL,NULL);
    test();
    MPI_Finalize();
    return 0;
}