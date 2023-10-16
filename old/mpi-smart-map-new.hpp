#ifndef _SMAR_MAP_
#define _SMAR_MAP_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "hvector.hpp"
#include <vector>
#include <string.h>

class ContigGet{
    private:
        int rank;
        int get_start;
        int n;
        int dest_start;
        bool init;
    
    public:
        inline ContigGet(){

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

        int* recv_ranks;
        int* recv_idxs;
        int* recv_counts;

        int* secondary_idxs;
        int* secondary_counts;
        int* secondary_sources;
        
        int* local_idxs;
        int* local_counts;
        int* local_sources;

        int* send_per_rank;
        int* recv_per_rank;

        int total_recv_size;
        int total_send_size;

    public:
        inline void reduce(std::vector<ContigGet>& gets, std::vector<SecondaryGet>& sgets, std::vector<ContigGet>& local_gets){
            bool reset_local = true;
            bool reset_global = false;
            for (int i = 0; i < n; i++){
                if (i % (n/100) == 0){
                    if (!comm_rank){
                        printf("%g done\n",(double)i/(double)n);
                    }
                }
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

            int* s_origin_ptr = s_origin;
            int* s_n_ptr = s_n;
            int* send_counts_ptr = send_counts;
            int* send_idxs_ptr = send_idxs;
            int* send_ranks_ptr = send_ranks;

            for (int i = 0; i < comm_size; i++){
                if (nsends[i] == 0)continue;
                for (int j = 0; j < nsends[i]; j++){
                    *send_counts_ptr++ = *s_n_ptr++;
                    *send_idxs_ptr++ = *s_origin_ptr++;
                    *send_ranks_ptr++ = i;
                }
            }
        }

        inline void fill_recv_info(int r_dest[], int r_n[], int nrecvs[]){
            recv_counts = (int*)malloc(sizeof(int)*total_recvs);
            recv_idxs = (int*)malloc(sizeof(int)*total_recvs);
            recv_ranks = (int*)malloc(sizeof(int)*total_recvs);

            int* r_dest_ptr = r_dest;
            int* r_n_ptr = r_n;
            int* recv_counts_ptr = recv_counts;
            int* recv_idxs_ptr = recv_idxs;
            int* recv_ranks_ptr = recv_ranks;

            for (int i = 0; i < comm_size; i++){
                if (nrecvs[i] == 0)continue;
                for (int j = 0; j < nrecvs[i]; j++){
                    *recv_counts_ptr++ = *r_n_ptr++;
                    *recv_idxs_ptr++ = *r_dest_ptr++;
                    *recv_ranks_ptr++ = i;
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

            if(!comm_rank)printf("starting reduce\n");
            reduce(gets,sgets,local_gets);

            if(!comm_rank)printf("starting fill local info\n");
            fill_local_info(local_gets);

            if(!comm_rank)printf("starting fill secondary info\n");
            fill_secondary_info(sgets);
            
            if(!comm_rank)printf("starting unify\n");
            unify(gets,unified_gets);

            if(!comm_rank)printf("starting fill nrecvs\n");
            int nrecvs[comm_size];
            int nsends[comm_size];

            total_recvs = fill_nrecvs(nrecvs,unified_gets);

            int r_origin[total_recvs];
            int r_n[total_recvs];
            int r_dests[total_recvs];

            if(!comm_rank)printf("starting fill origin\n");
            fill_origin(r_origin,r_n,r_dests,nrecvs,unified_gets);

            if(!comm_rank)printf("starting fill recv info\n");
            fill_recv_info(r_dests,r_n,nrecvs);

            if(!comm_rank)printf("starting fill nsends\n");
            total_sends = fill_nsends(nrecvs,nsends);

            int s_origin[total_sends];
            int s_n[total_sends];

            if(!comm_rank)printf("starting send origin\n");
            send_origin(r_origin,r_n,s_origin,s_n,nrecvs,nsends);

            if(!comm_rank)printf("starting fill send info\n");
            fill_send_info(s_origin,s_n,nsends);

            if(!comm_rank)printf("starting finalizing\n");

            for (int i = 0; i < comm_size; i++){
                send_per_rank[i] = 0;
                recv_per_rank[i] = 0;
            }

            total_send_size = 0;
            for (int i = 0; i < total_sends; i++){
                int count = send_counts[i];
                int rank = send_ranks[i];
                send_per_rank[rank] += count;
                total_send_size += count;
            }
            
            total_recv_size = 0;
            for (int i = 0; i < total_recvs; i++){
                int count = recv_counts[i];
                int rank = recv_ranks[i];
                recv_per_rank[rank] += count;
                total_recv_size += count;
            }

            /*for (int i = 0; i < comm_size; i++){
                printf("rank %d send %d to   rank %d\n",comm_rank,send_per_rank[i],i);
                printf("rank %d recv %d from rank %d\n",comm_rank,recv_per_rank[i],i);
            }*/

            /*for (int i = 0; i < total_sends; i++){
                printf("rank %d sending %d from %d to rank %d\n",comm_rank,send_counts[i],send_idxs[i],send_ranks[i]);
            }

            for (int i = 0; i < total_recvs; i++){
                printf("rank %d recieving %d to %d from rank %d\n",comm_rank,recv_counts[i],recv_idxs[i],recv_ranks[i]);
            }

            for (int i = 0; i < total_secondary; i++){
                printf("rank %d mapping %d items from %d to %d\n",comm_rank,secondary_counts[i],secondary_sources[i],secondary_idxs[i]);
            }*/
            
        }

        inline SmartMap(MPI_Comm comm_, int n_, MapT map_) : comm(comm_), n(n_), map(map_),send_counts(NULL),send_idxs(NULL),send_ranks(NULL),
                                                                        recv_counts(NULL),recv_idxs(NULL),recv_ranks(NULL),
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
            if(recv_counts)free(recv_counts);
            if(recv_idxs)free(recv_idxs);
            if(recv_ranks)free(recv_ranks);
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

            //We might need to wait before we move the local stuff over
            MPI_Request send_reqs[comm_size];
            bool wait_on[comm_size];
            for (int i = 0; i < comm_size; i++){
                wait_on[i] = false;
            }

            int n_buffered = 0;
            int current_rank_start = 0;
            int current_n = 0;

            //for every batch we need to send
            for (int i = 0; i < total_sends; i++){
                int n = send_counts[i];
                int in_idx = send_idxs[i];
                T* in_buff = &in[in_idx];
                int dest = send_ranks[i];
                T* out_buff = &out[n_buffered];
                //move it so it is contiguous
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;

                //make sure we wait if we need to wait (I think this works but idk)
                if ((in_idx + n) >= total_recv_size){
                    wait_on[dest] = true;
                }

                current_n += n;
                //If we have packed everything that we need to send to dest, then send it all to dest
                if (current_n == send_per_rank[dest]){
                    T* in_buff = &out[current_rank_start];
                    MPI_Isend(in_buff,current_n*sizeof(T),MPI_BYTE,dest,0,comm,&send_reqs[dest]);
                    current_rank_start = n_buffered;
                    current_n = 0;
                }
            }
            
            //Save this for later, we need to know where we moved the local data to
            int local_start = n_buffered;

            //Now for every local transfer (i.e. one local to this rank), buffer it all in out
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                int in_idx = local_sources[i];
                T* in_buff = &in[in_idx];
                int out_idx = local_idxs[i];
                T* out_buff = &out[n_buffered];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }
            
            n_buffered = 0;

            MPI_Request recv_reqs[comm_size];
            //Begin our recv requests
            for (int i = 0; i < comm_size; i++){
                int n = recv_per_rank[i];
                if (n==0)continue;
                T* out_buff = &in[n_buffered];
                //printf("rank %d recv from rank %d\n",comm_rank,i);
                MPI_Irecv(out_buff,n*sizeof(T),MPI_BYTE,i,0,comm,&recv_reqs[i]);
                n_buffered += n;
            }

            //we might have multiple sub requests from each rank, but we only want to wait once
            int wait_count[comm_size];
            for (int i = 0; i < comm_size; i++){
                wait_count[i] = 0;
            }

            //Wait if we need to for the local moves
            for (int i = 0; i < comm_size; i++){
                if (send_per_rank[i] == 0)continue;
                if (wait_on[i])
                    MPI_Wait(&send_reqs[i],MPI_STATUS_IGNORE);
            }

            //now we move the local stuff back into in, leaving space for the transfered stuff
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                T* in_buff = &out[local_start];
                T* out_buff = &in[n_buffered];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
                local_start += n;
            }

            //and then move it straight back
            //we have to do this because we might overrite the stuff if we do it in one go
            n_buffered = total_recv_size;
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                T* in_buff = &in[n_buffered];
                int out_idx = local_idxs[i];
                T* out_buff = &out[out_idx];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }

            //wait on the rest of the send requests, to make sure we dont overrite anything
            for (int i = 0; i < comm_size; i++){
                if (send_per_rank[i] == 0)continue;
                if (!wait_on[i])
                    MPI_Wait(&send_reqs[i],MPI_STATUS_IGNORE);
            }

            //and finally, we wait on each request once per rank, and move it to out
            n_buffered = 0;
            for (int i = 0; i < total_recvs; i++){
                int r_rank = recv_ranks[i];
                if(wait_count[r_rank] == 0){
                    MPI_Wait(&recv_reqs[r_rank],MPI_STATUS_IGNORE);
                    wait_count[r_rank]++;
                }
                int n = recv_counts[i];
                int out_idx = recv_idxs[i];
                T* in_buff = &in[n_buffered];
                T* out_buff = &out[out_idx];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }
            
        }

        template<class T>
        inline void backward(T* in, T* out){

            //We might need to wait before we move the local stuff over
            MPI_Request send_reqs[comm_size];
            bool wait_on[comm_size];
            for (int i = 0; i < comm_size; i++){
                wait_on[i] = false;
            }

            int n_buffered = 0;
            int current_rank_start = 0;
            int current_n = 0;

            //for every batch we need to send
            for (int i = 0; i < total_recvs; i++){
                int n = recv_counts[i];
                int in_idx = recv_idxs[i];
                T* in_buff = &in[in_idx];
                int dest = recv_ranks[i];
                T* out_buff = &out[n_buffered];
                //move it so it is contiguous
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;

                //make sure we wait if we need to wait (I think this works but idk)
                if ((in_idx + n) >= total_send_size){
                    wait_on[dest] = true;
                }

                current_n += n;
                //If we have packed everything that we need to send to dest, then send it all to dest
                if (current_n == recv_per_rank[dest]){
                    T* in_buff = &out[current_rank_start];
                    MPI_Isend(in_buff,current_n*sizeof(T),MPI_BYTE,dest,0,comm,&send_reqs[dest]);
                    current_rank_start = n_buffered;
                    current_n = 0;
                }
            }
            
            //Save this for later, we need to know where we moved the local data to
            int local_start = n_buffered;

            //Now for every local transfer (i.e. one local to this rank), buffer it all in out
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                //int in_idx = local_sources[i];
                int out_idx = local_idxs[i];
                T* in_buff = &in[out_idx];
                
                T* out_buff = &out[n_buffered];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }
            
            n_buffered = 0;

            MPI_Request recv_reqs[comm_size];
            //Begin our recv requests
            for (int i = 0; i < comm_size; i++){
                int n = send_per_rank[i];
                if (n==0)continue;
                T* out_buff = &in[n_buffered];
                //printf("rank %d recv from rank %d\n",comm_rank,i);
                MPI_Irecv(out_buff,n*sizeof(T),MPI_BYTE,i,0,comm,&recv_reqs[i]);
                n_buffered += n;
            }

            //we might have multiple sub requests from each rank, but we only want to wait once
            int wait_count[comm_size];
            for (int i = 0; i < comm_size; i++){
                wait_count[i] = 0;
            }

            //Wait if we need to for the local moves
            for (int i = 0; i < comm_size; i++){
                if (recv_per_rank[i] == 0)continue;
                if (wait_on[i])
                    MPI_Wait(&send_reqs[i],MPI_STATUS_IGNORE);
            }

            //now we move the local stuff back into in, leaving space for the transfered stuff
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                T* in_buff = &out[local_start];
                T* out_buff = &in[n_buffered];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
                local_start += n;
            }

            //and then move it straight back
            //we have to do this because we might overrite the stuff if we do it in one go
            n_buffered = total_send_size;
            for (int i = 0; i < total_local; i++){
                int n = local_counts[i];
                T* in_buff = &in[n_buffered];
                int out_idx = local_sources[i];
                T* out_buff = &out[out_idx];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }

            //wait on the rest of the send requests, to make sure we dont overrite anything
            for (int i = 0; i < comm_size; i++){
                if (recv_per_rank[i] == 0)continue;
                if (!wait_on[i])
                    MPI_Wait(&send_reqs[i],MPI_STATUS_IGNORE);
            }

            //and finally, we wait on each request once per rank, and move it to out
            n_buffered = 0;
            for (int i = 0; i < total_sends; i++){
                int r_rank = send_ranks[i];
                if(wait_count[r_rank] == 0){
                    MPI_Wait(&recv_reqs[r_rank],MPI_STATUS_IGNORE);
                    wait_count[r_rank]++;
                }
                int n = send_counts[i];
                int out_idx = send_idxs[i];
                T* in_buff = &in[n_buffered];
                T* out_buff = &out[out_idx];
                memcpy(out_buff,in_buff,n*sizeof(T));
                n_buffered += n;
            }
            
        }

};

#endif