#ifndef _SMAR_MAP_
#define _SMAR_MAP_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "hvector.hpp"
#include <vector>

struct map_return_t{
    int rank;
    int idx;
};

inline map_return_t make_map_return(int rank, int idx){
    map_return_t out;
    out.rank = rank;
    out.idx = idx;
    return out;
}

class ContigGet{
    private:
        int rank;
        int start;
        int n;
        bool init;
    
    public:
        inline ContigGet() : rank(-1), start(-1), n(-1), init(false){

        }

        inline ContigGet(int rank_, int start_, int n_ = 1) : rank(rank_), start(start_), n(n_), init(true){
            
        }

        inline ContigGet(map_return_t tmp, int n_ = 1) : rank(tmp.rank), start(tmp.idx), n(n_), init(true){
            printf("init contig get, rank %d start %d n %d\n",rank,start,n);
        }

        inline ~ContigGet(){}

        inline bool add(map_return_t in){
            if (!init)return false;
            if (in.rank != rank)return false;
            if (in.idx != (start + n))return false;
            n++;
            return true;
        }

        inline bool is_init(){
            return init;
        }

        inline bool is_secondary(map_return_t in){
            if (!init)return false;
            if (in.rank != rank)return false;
            int off = in.idx - start;
            return ((off >= 0) && (off < n));
        }
};

class SecondaryGet{
    private:
        ContigGet& origin;
        int offset;
    
    public:
        inline SecondaryGet(ContigGet& origin_, int offset_) : origin(origin_), offset(offset_){

        }

        inline ~SecondaryGet(){}
};

template<class T, class MapT>
class SmartMap{
    private:
        int comm_size;
        int comm_rank;
        MPI_Comm comm;
        std::vector<ContigGet> gets;
        std::vector<SecondaryGet> secondaryget;
        MapT map;
        T* raw;
        int n;

    public:
        inline void setup_map(){
            if(comm_rank == 0){
                printf("Setting up map!\n");
                ContigGet cur;
                for (int i = 0; i < n; i++){
                    map_return_t out = map.get(i);
                    if (cur.add(out)){
                        printf("found contiguous map for idx %d\n",i);
                        continue;
                    }
                    bool found_secondary = false;
                    for (ContigGet m : gets){
                        if (m.is_secondary(out)){
                            printf("found secondary map for idx %d\n",i);
                            found_secondary = true;
                            break;
                        }
                    }
                    if(found_secondary){
                        continue;
                    }
                    cur = ContigGet(out);
                    gets.push_back(cur);
                    printf("map %d from rank %d idx %d\n",i,out.rank,out.idx);
                }
            }
            
        }

        inline SmartMap(MPI_Comm comm_, T* raw_, int n_, MapT map_) : comm(comm_), raw(raw_), n(n_), map(map_){
            MPI_Comm_size(comm,&comm_size);
            MPI_Comm_rank(comm,&comm_rank);
            setup_map();
        }

        inline ~SmartMap(){}

};

#endif