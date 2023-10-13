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
        int get_start;
        int n;
        int dest_start;
        bool init;
    
    public:
        inline ContigGet() : init(false){

        }

        inline ContigGet(map_return_t tmp, int out_idx) : rank(tmp.rank), get_start(tmp.idx), dest_start(out_idx), n(1), init(true){
            //printf("init contig get, rank %d start %d n %d, local_start = %d\n",rank,get_start,n,dest_start);
        }

        inline ~ContigGet(){
            //printf("destructor called\n");
        }

        inline bool add(map_return_t in){
            if (!init)return false;
            if (in.rank != rank)return false;
            if (in.idx != (get_start + n))return false;
            n++;
            return true;
        }

        inline bool is_init(){
            return init;
        }

        inline bool is_secondary(map_return_t in){
            
            //printf("Check Secondary! %d %d, %d %d %d\n",in.rank,in.idx,rank,get_start,n);
            if (!init)return false;
            if (in.rank != rank)return false;
            int off = in.idx - get_start;
            return ((off >= 0) && (off < n));
        }

        inline int get_secondary_offset(int i){
            return i - get_start;
        }

        inline int get_rank(){
            return rank;
        }
};

class SecondaryGet{
    private:
        ContigGet& origin;
        int offset;
        int n;
        int dest_start;
    
    public:
        inline SecondaryGet(ContigGet& origin_, int offset_, int dest_) : origin(origin_), offset(offset_), n(1), dest_start(dest_){

        }

        inline ~SecondaryGet(){

        }

        inline bool add(map_return_t in){
            if (!origin.is_secondary(in))return false;
            if (in.idx != (n+offset))return false;
            n++;
            return true;
        }

};

struct send_recv_t{
    int rank;
    int start;
    int n;
    int dest;
};

template<class T, class MapT>
class SmartMap{
    private:
        int comm_size;
        int comm_rank;
        MPI_Comm comm;
        MapT map;
        T* raw;
        send_recv_t** recv_data;
        send_recv_t** send_data;
        int n;

    public:
        inline void setup_map(){

            std::vector<ContigGet> gets;
            std::vector<SecondaryGet> sgets;

            for (int i = 0; i < n; i++){
                map_return_t out = map.get(i);
                if (gets.size() != 0){
                    if (gets[gets.size()-1].add(out)){
                        continue;
                    }
                }
                bool found_secondary = false;

                for (int j = 0; j < sgets.size(); j++){
                    if (sgets[j].add(out)){
                        break;
                    }
                }

                if (found_secondary)continue;
                
                for (int j = 0; j < gets.size(); j++){
                    if (gets[j].is_secondary(out)){
                        SecondaryGet tmp(gets[j],gets[j].get_secondary_offset(out.idx),i);
                        sgets.push_back(tmp);
                        found_secondary = true;
                        break;
                    }
                }

                if (found_secondary)continue;

                gets.push_back(ContigGet(out,i));
            }

            int nrecvs[comm_size] = {0};
            for (int j = 0; j < gets.size(); j++){
                nrecvs[gets[j].get_rank()]++;
            }

            //int* recvdata[comm_size];

            for (int i = 0; i < comm_size; i++){
                if(nrecvs[i] == 0)continue;
                printf("rank %d expecting %d recvs from rank %d\n",comm_rank,nrecvs[i],i);
                recv_data[i] = (send_recv_t*)malloc(sizeof(send_recv_t)*nrecvs[i]);
            }
            
            int nsends[comm_size] = {0};

            MPI_Alltoall(nrecvs,1,MPI_INT,nsends,1,MPI_INT,comm);

            for (int i = 0; i < comm_size; i++){
                if(nsends[i] == 0)continue;
                printf("rank %d sending %d blocks to rank %d\n",comm_rank,nsends[i],i);
            }

            for (int i = 0; i < comm_size; i++){
                if(nrecvs[i] == 0)continue;
                //free(recvdata[i]);
            }
            
        }

        inline SmartMap(MPI_Comm comm_, T* raw_, int n_, MapT map_) : comm(comm_), raw(raw_), n(n_), map(map_){
            MPI_Comm_size(comm,&comm_size);
            MPI_Comm_rank(comm,&comm_rank);
            recv_data = (send_recv_t**)malloc(sizeof(send_recv_t*)*comm_size);
            send_data = (send_recv_t**)malloc(sizeof(send_recv_t*)*comm_size);
            for (int i = 0; i < comm_size; i++){
                recv_data[i] = NULL;
                send_data[i] = NULL;
            }
            setup_map();
        }

        inline ~SmartMap(){
            for (int i = 0; i < comm_size; i++){
                if(send_data[i])free(send_data[i]);
                if(recv_data[i])free(recv_data[i]);
            }
            free(send_data);
            free(recv_data);
        }

};

#endif