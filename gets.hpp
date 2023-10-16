#include <stdlib.h>
#include <stdio.h>

struct map_return_t{
    int rank;
    int src_idx;
    int dest_idx;
};

inline map_return_t make_map_return(int rank, int src_idx, int dest_idx){
    map_return_t out;
    out.rank = rank;
    out.src_idx = src_idx;
    out.dest_idx = dest_idx;
    return out;
}

struct sget_t{
    int src;
    int dest;
    int stride;
    int n;
    sget_t* next;
};

inline sget_t* make_sget(int src, int dest, int stride, int n){
    sget_t* out = (sget_t*)malloc(sizeof(sget_t));
    out->src = src;
    out->dest = dest;
    out->stride = stride;
    out->n = n;
    return out;
}

inline void free_sget(sget_t* out){
    sget_t* tmp = out;
    while (tmp->next){
        sget_t* prev = tmp;
        tmp = tmp->next;
        free(prev);
    }
    free(tmp);
}

inline void append_sget(sget_t* base, sget_t* add){
    sget_t* tmp = base;
    while (tmp->next){
        tmp = tmp->next;
    }
    tmp->next = add;
}

inline sget_t* make_empty_sget(){
    sget_t* out = (sget_t*)malloc(sizeof(sget_t));
    out->src = -1;
    out->dest = -1;
    out->stride = -1;
    out->n = -1;
    return out;
}

struct get_t{
    int rank;
    int src;
    int dest;
    int stride;
    int n;
    get_t* next;
};

inline get_t* make_get(int rank, int src, int dest, int stride, int n){
    get_t* out = (get_t*)malloc(sizeof(get_t));
    out->rank = rank;
    out->src = src;
    out->dest = dest;
    out->stride = stride;
    out->n = n;
    out->next = NULL;
    return out;
}

inline get_t* make_empty_get(){
    get_t* out = (get_t*)malloc(sizeof(get_t));
    out->rank = -1;
    out->next = NULL;
    return out;
}

inline void free_get(get_t* out){
    get_t* tmp = out;
    while (tmp->next){
        get_t* prev = tmp;
        tmp = tmp->next;
        free(prev);
    }
    free(tmp);
}

inline void append_get(get_t* base, get_t* add){
    get_t* tmp = base;
    while (tmp->next){
        tmp = tmp->next;
    }
    tmp->next = add;
}

inline get_t* add_get(get_t* base, map_return_t ret){
    if ((base)->rank == ret.rank){
        if ((base)->n == 1){
            if ((base)->src < ret.src_idx){
                (base)->stride = ret.src_idx - (base)->src;
                (base)->n++;
                return base;
            }
        }
        if (((base)->n * (base)->stride + (base)->src) == ret.src_idx){
            (base)->n++;
            return base;
        }
    }
    
    if(base->rank != -1){
        get_t* new_get = make_get(ret.rank,ret.src_idx,ret.dest_idx,1,1);
        append_get(base,new_get);
        return new_get;
    }
    base->rank = ret.rank;
    base->n = 1;
    base->next = NULL;
    base->dest = ret.dest_idx;
    base->stride = 1;
    base->src = ret.src_idx;
    return base;
}

template<class Map>
inline get_t* find_gets(int n, Map& map){
    get_t* init = make_empty_get();
    get_t* cur = init;
    for (int i = 0; i < n; i++){
        map_return_t ret = map.map(i);
        cur = add_get(cur,ret);
    }
    return init;
}