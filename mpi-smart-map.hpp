#ifndef _SMAR_MAP_
#define _SMAR_MAP_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "hvector.hpp"
#include <vector>
#include <string.h>

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
            //printf("add?\n");
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

        inline int send_start(){
            return get_start;
        }

        inline int recv_start(){
            return dest_start;
        }

        inline int send_n(){
            return n;
        }
};

class UnifiedGet{
    private:
        std::vector<ContigGet> my_gets;
        int rank;
    
    public:
        inline UnifiedGet(int rank_) : rank(rank_){

        }

        inline ~UnifiedGet(){};

        inline void add(ContigGet get){
            my_gets.push_back(get);
        }

        inline int get_n(){
            return my_gets.size();
        }

        inline void fill_sends(int* origin, int* n, int* dest){
            for (int i = 0; i < my_gets.size(); i++){
                origin[i] = my_gets[i].send_start();
                n[i] = my_gets[i].send_n();
                dest[i] = my_gets[i].recv_start();
            }
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
            //printf("ADD?\n");
            n++;
            return true;
        }

        inline int get_source(){
            return offset;
        }

        inline int get_count(){
            return n;
        }

        inline int get_idx(){
            return dest_start;
        }

};

struct send_recv_t{
    int rank;
    int start;
    int n;
    int dest;
};

template<class MapT>
class SmartMap{
    private:
        int comm_size;
        int comm_rank;
        MPI_Comm comm;
        MapT map;
        int total_recvs;
        int total_sends;
        int total_secondary;
        int total_local;
        int n;

        int* send_ranks;
        int* send_idxs;
        int* send_counts;
        int* send_tags;

        int* recv_ranks;
        int* recv_idxs;
        int* recv_counts;
        int* recv_tags;

        int* secondary_idxs;
        int* secondary_counts;
        int* secondary_sources;
        
        int* local_idxs;
        int* local_counts;
        int* local_sources;

        int* send_per_rank;
        int* recv_per_rank;

        int total_recv_size;

    public:
        inline void reduce(std::vector<ContigGet>& gets, std::vector<SecondaryGet>& sgets, std::vector<ContigGet>& local_gets){
            bool reset_local = true;
            bool reset_global = false;
            for (int i = 0; i < n; i++){
                map_return_t out = map.get(i);
                if (out.rank == comm_rank){
                    if ((local_gets.size() != 0) && (!reset_local)){
                        if (local_gets[local_gets.size()-1].add(out)){
                            //if(!comm_rank)printf("Adding contig! %d, %d\n",out.idx,i);
                            continue;
                        }
                    }
                    //if(!comm_rank)printf("Should be new contig? %d, %d\n",out.idx,i);
                    local_gets.push_back(ContigGet(out,i));
                    reset_local = false;
                    reset_global = true;
                    continue;
                }

                reset_local = true;

                if (!reset_global){
                    if (gets.size() != 0){
                        if (gets[gets.size()-1].add(out)){
                            continue;
                        }
                    }
                }

                reset_global = false;
                

                bool found_secondary = false;

                for (int j = 0; j < sgets.size(); j++){
                    if (sgets[j].add(out)){
                        found_secondary = true;
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
        }

        inline void unify(std::vector<ContigGet>& gets, std::vector<UnifiedGet>& unified_gets){
            unified_gets.reserve(comm_size);
            for (int i = 0; i < comm_size; i++){
                unified_gets.push_back(UnifiedGet(i));
            }
            for (auto m : gets){
                unified_gets[m.get_rank()].add(m);
            }
        }

        inline int fill_nrecvs(int nrecvs[], std::vector<UnifiedGet>& unified_gets){
            int total = 0;
            for (int i = 0; i < comm_size; i++){
                int n = unified_gets[i].get_n();
                nrecvs[i] = n;
                total += n;
            }
            return total;
        }

        inline void fill_origin(int origin[], int n[], int dests[], int nrecvs[], std::vector<UnifiedGet>& unified_gets){
            int* origin_ptr = origin;
            int* n_ptr = n;
            int* dest_ptr = dests;

            for (int i = 0; i < comm_size; i++){
                if (nrecvs[i] == 0)continue;
                unified_gets[i].fill_sends(origin_ptr,n_ptr,dest_ptr);
                origin_ptr += nrecvs[i];
                n_ptr += nrecvs[i];
                dest_ptr += nrecvs[i];
            }
        }

        inline int fill_nsends(int nrecvs[], int nsends[]){
            MPI_Alltoall(nrecvs,1,MPI_INT,nsends,1,MPI_INT,comm);
            int total = 0;
            for (int i = 0; i < comm_size; i++){
                total += nsends[i];
            }
            return total;
        }

        inline void send_origin(int r_origin[], int r_n[], int s_origin[], int s_n[], int nrecvs[], int nsends[]){

            int s_disp[comm_size];
            int r_disp[comm_size];
            s_disp[0] = 0;
            r_disp[0] = 0;
            for (int i = 1; i < comm_size; i++){
                s_disp[i] = s_disp[i-1] + nrecvs[i-1];
                r_disp[i] = r_disp[i-1] + nsends[i-1];
            }

            MPI_Alltoallv(r_origin,nrecvs,s_disp,MPI_INT,s_origin,nsends,r_disp,MPI_INT,comm);
            MPI_Alltoallv(r_n,nrecvs,s_disp,MPI_INT,s_n,nsends,r_disp,MPI_INT,comm);

        }

        inline void fill_send_info(int s_origin[], int s_n[], int nsends[]){
            send_counts = (int*)malloc(sizeof(int)*total_sends);
            send_idxs = (int*)malloc(sizeof(int)*total_sends);
            send_ranks = (int*)malloc(sizeof(int)*total_sends);
            send_tags = (int*)malloc(sizeof(int)*total_sends);

            int* s_origin_ptr = s_origin;
            int* s_n_ptr = s_n;
            int* send_counts_ptr = send_counts;
            int* send_idxs_ptr = send_idxs;
            int* send_ranks_ptr = send_ranks;
            int* send_tags_ptr = send_tags;

            for (int i = 0; i < comm_size; i++){
                if (nsends[i] == 0)continue;
                for (int j = 0; j < nsends[i]; j++){
                    *send_counts_ptr++ = *s_n_ptr++;
                    *send_idxs_ptr++ = *s_origin_ptr++;
                    *send_ranks_ptr++ = i;
                    *send_tags_ptr++ = j;
                }
            }
        }

        inline void fill_recv_info(int r_dest[], int r_n[], int nrecvs[]){
            recv_counts = (int*)malloc(sizeof(int)*total_recvs);
            recv_idxs = (int*)malloc(sizeof(int)*total_recvs);
            recv_ranks = (int*)malloc(sizeof(int)*total_recvs);
            recv_tags = (int*)malloc(sizeof(int)*total_recvs);

            int* r_dest_ptr = r_dest;
            int* r_n_ptr = r_n;
            int* recv_counts_ptr = recv_counts;
            int* recv_idxs_ptr = recv_idxs;
            int* recv_ranks_ptr = recv_ranks;
            int* recv_tags_ptr = recv_tags;

            for (int i = 0; i < comm_size; i++){
                if (nrecvs[i] == 0)continue;
                for (int j = 0; j < nrecvs[i]; j++){
                    *recv_counts_ptr++ = *r_n_ptr++;
                    *recv_idxs_ptr++ = *r_dest_ptr++;
                    *recv_ranks_ptr++ = i;
                    *recv_tags_ptr++ = j;
                }
            }
        }

        inline void fill_secondary_info(std::vector<SecondaryGet>& sgets){
            total_secondary = sgets.size();
            secondary_counts = (int*)malloc(sizeof(int) * total_secondary);
            secondary_idxs = (int*)malloc(sizeof(int) * total_secondary);
            secondary_sources = (int*)malloc(sizeof(int) * total_secondary);

            for (int i = 0; i < total_secondary; i++){
                secondary_counts[i] = sgets[i].get_count();
                secondary_idxs[i] = sgets[i].get_idx();
                secondary_sources[i] = sgets[i].get_source();
            }
        }

        inline void fill_local_info(std::vector<ContigGet>& local_gets){
            total_local = local_gets.size();
            local_counts = (int*)malloc(sizeof(int) * total_local);
            local_idxs = (int*)malloc(sizeof(int) * total_local);
            local_sources = (int*)malloc(sizeof(int) * total_local);

            for (int i = 0; i < total_local; i++){
                local_counts[i] = local_gets[i].send_n();
                local_idxs[i] = local_gets[i].recv_start();
                local_sources[i] = local_gets[i].send_start();
                //printf("local %d %d %d\n",local_counts[i],local_idxs[i],)
            }
        }

        inline void setup_map(){
            
            std::vector<ContigGet> gets;
            std::vector<ContigGet> local_gets;
            std::vector<SecondaryGet> sgets;
            std::vector<UnifiedGet> unified_gets;

            reduce(gets,sgets,local_gets);

            fill_local_info(local_gets);

            fill_secondary_info(sgets);
            
            unify(gets,unified_gets);

            int nrecvs[comm_size];
            int nsends[comm_size];

            total_recvs = fill_nrecvs(nrecvs,unified_gets);

            int r_origin[total_recvs];
            int r_n[total_recvs];
            int r_dests[total_recvs];

            fill_origin(r_origin,r_n,r_dests,nrecvs,unified_gets);

            fill_recv_info(r_dests,r_n,nrecvs);

            total_sends = fill_nsends(nrecvs,nsends);

            int s_origin[total_sends];
            int s_n[total_sends];

            send_origin(r_origin,r_n,s_origin,s_n,nrecvs,nsends);

            fill_send_info(s_origin,s_n,nsends);

            for (int i = 0; i < comm_size; i++){
                send_per_rank[i] = 0;
                recv_per_rank[i] = 0;
            }

            for (int i = 0; i < total_sends; i++){
                int count = send_counts[i];
                int rank = send_ranks[i];
                send_per_rank[rank] += count;
            }
            
            total_recv_size = 0;
            for (int i = 0; i < total_recvs; i++){
                int count = recv_counts[i];
                int rank = recv_ranks[i];
                recv_per_rank[rank] += count;
                total_recv_size += count;
            }

            for (int i = 0; i < comm_size; i++){
                printf("rank %d send %d to   rank %d\n",comm_rank,send_per_rank[i],i);
                printf("rank %d recv %d from rank %d\n",comm_rank,recv_per_rank[i],i);
            }

            /*for (int i = 0; i < total_sends; i++){
                printf("rank %d sending %d from %d to rank %d (tag = %d)\n",comm_rank,send_counts[i],send_idxs[i],send_ranks[i],send_tags[i]);
            }

            for (int i = 0; i < total_recvs; i++){
                printf("rank %d recieving %d to %d from rank %d (tag = %d)\n",comm_rank,recv_counts[i],recv_idxs[i],recv_ranks[i],recv_tags[i]);
            }

            for (int i = 0; i < total_secondary; i++){
                printf("rank %d mapping %d items from %d to %d\n",comm_rank,secondary_counts[i],secondary_sources[i],secondary_idxs[i]);
            }*/
            
        }

        inline SmartMap(MPI_Comm comm_, int n_, MapT map_) : comm(comm_), n(n_), map(map_),send_counts(NULL),send_idxs(NULL),send_ranks(NULL),send_tags(NULL),
                                                                        recv_counts(NULL),recv_idxs(NULL),recv_ranks(NULL),recv_tags(NULL),
                                                                        secondary_counts(NULL), secondary_idxs(NULL), secondary_sources(NULL),
                                                                        local_counts(NULL),local_idxs(NULL),local_sources(NULL),send_per_rank(NULL),recv_per_rank(NULL){
            MPI_Comm_size(comm,&comm_size);
            MPI_Comm_rank(comm,&comm_rank);
            send_per_rank = (int*)malloc(sizeof(int)*comm_size);
            recv_per_rank = (int*)malloc(sizeof(int)*comm_size);
            setup_map();
            MPI_Barrier(comm);
            for (int i = 0; i < comm_size; i++){
                if(comm_rank == i){
                    printf("rank %d: total_sends = %d, total_recvs = %d, total_secondary = %d, total_local = %d\n",comm_rank,total_sends,total_recvs,total_secondary,total_local);
                }
                MPI_Barrier(comm);
            }
            MPI_Barrier(comm);
        }

        inline ~SmartMap(){
            if(send_counts)free(send_counts);
            if(send_idxs)free(send_idxs);
            if(send_ranks)free(send_ranks);
            if(send_tags)free(send_tags);
            if(recv_counts)free(recv_counts);
            if(recv_idxs)free(recv_idxs);
            if(recv_ranks)free(recv_ranks);
            if(recv_tags)free(recv_tags);
            if(secondary_counts)free(secondary_counts);
            if(secondary_idxs)free(secondary_idxs);
            if(secondary_sources)free(secondary_sources);
            if(local_counts)free(local_counts);
            if(local_idxs)free(local_idxs);
            if(local_sources)free(local_sources);
            if(send_per_rank)free(send_per_rank);
            if(recv_per_rank)free(recv_per_rank);
        }

        template<class T>
        inline void forward(T* in, T* out){
            
            int n_buffered = 0;

            for (int i = 0; i < total_sends; i++){
                int n = send_counts[i];
                int in_idx = send_idxs[i];
                T* in_buff = &in[in_idx];
                int dest = send_ranks[i];
                T* out_buff = &out[n_buffered];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }
            int local_start = n_buffered;

            n_buffered = 0;
            
            for (int i = 0; i < comm_size; i++){
                int n = send_per_rank[i];
                if (n == 0)continue;
                T* in_buff = &out[n_buffered];
                MPI_Request req;
                MPI_Isend(in_buff,n*sizeof(T),MPI_BYTE,i,0,comm,&req);
                MPI_Request_free(&req);
                printf("rank %d sending to rank %d\n",comm_rank,i);
                n_buffered += n;
            }

            n_buffered = local_start;
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                int in_idx = local_sources[i];
                T* in_buff = &in[in_idx];
                int out_idx = local_idxs[i];
                T* out_buff = &out[n_buffered];
                //if(!comm_rank)printf("rank %d moving local %d to %d, size %d\n",comm_rank,in_idx,out_idx,n);
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }
            
            n_buffered = 0;

            for (int i = 0; i < comm_size; i++){
                int n = recv_per_rank[i];
                if (n==0)continue;
                T* out_buff = &in[n_buffered];
                MPI_Request req;
                printf("rank %d recv from rank %d\n",comm_rank,i);
                MPI_Irecv(out_buff,n*sizeof(T),MPI_BYTE,i,0,comm,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);
                n_buffered += n;
            }

            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                //int in_idx = local_sources[i];
                T* in_buff = &out[local_start];
                //int out_idx = local_idxs[i];
                T* out_buff = &in[n_buffered];
                //if(!comm_rank)printf("rank %d moving local %d to %d, size %d\n",comm_rank,in_idx,out_idx,n);
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
                local_start += n;
            }

            n_buffered = 0;
            for (int i = 0; i < total_recvs; i++){
                int n = recv_counts[i];
                int out_idx = recv_idxs[i];
                T* in_buff = &in[n_buffered];
                T* out_buff = &out[out_idx];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }

            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                //int in_idx = local_sources[i];
                T* in_buff = &in[n_buffered];
                int out_idx = local_idxs[i];
                T* out_buff = &out[out_idx];
                //if(!comm_rank)printf("rank %d moving local %d to %d, size %d\n",comm_rank,in_idx,out_idx,n);
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }

            /*for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                int in_idx = local_sources[i];
                T* in_buff = &out[n_buffered];
                int out_idx = local_idxs[i];
                T* out_buff = &in[out_idx];
                //if(!comm_rank)printf("rank %d moving local %d to %d, size %d\n",comm_rank,in_idx,out_idx,n);
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }*/

            /*for (int i = 0; i < total_sends; i++){
                int n = send_counts[i];
                int in_idx = send_idxs[i];
                T* in_buff = &in[in_idx];
                int dest = send_ranks[i];
                int tag = send_tags[i];
                //if(!comm_rank)printf("rank %d sending %d from %d to rank %d (tag = %d)\n",comm_rank,n,in_idx,dest,tag);
                MPI_Request req;
                MPI_Isend(in_buff,n * sizeof(T),MPI_BYTE,dest,tag,comm,&req);
                MPI_Request_free(&req);
                //if(!comm_rank)printf("rank %d sent %d from %d to rank %d (tag = %d)\n",comm_rank,n,in_idx,dest,tag);
            }

            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                int in_idx = local_sources[i];
                T* in_buff = &in[in_idx];
                int out_idx = local_idxs[i];
                T* out_buff = &out[out_idx];
                //if(!comm_rank)printf("rank %d moving local %d to %d, size %d\n",comm_rank,in_idx,out_idx,n);
                memcpy(out_buff,in_buff,n*sizeof(T));
            }

            for (int i = 0; i < total_recvs; i++){
                int n = recv_counts[i];
                int out_idx = recv_idxs[i];
                T* out_buff = &out[out_idx];
                int src = recv_ranks[i];
                int tag = recv_tags[i];
                //if(!comm_rank)printf("rank %d recieving %d to %d from rank %d (tag = %d)\n",comm_rank,n,out_idx,src,tag);
                MPI_Request req;
                MPI_Irecv(out_buff,n*sizeof(T),MPI_BYTE,src,tag,comm,&req);
                MPI_Wait(&req,MPI_STATUS_IGNORE);
                //if(!comm_rank)printf("rank %d recieved %d to %d from rank %d (tag = %d)\n",comm_rank,n,out_idx,src,tag);
            }

            for (int i = 0; i < total_secondary; i++){
                int n = secondary_counts[i];
                int source_idx = secondary_sources[i];
                int dest_idx = secondary_idxs[i];
                //if(!comm_rank)printf("rank %d mapping %d items from %d to %d\n",comm_rank,n,source_idx,dest_idx);
                memcpy(&out[dest_idx],&out[source_idx],n*sizeof(T));
                //if(!comm_rank)printf("rank %d mapped %d items from %d to %d\n",comm_rank,n,source_idx,dest_idx);
            }*/
            
        }

};

#endif