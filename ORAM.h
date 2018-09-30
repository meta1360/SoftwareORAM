//
// Created by meysam on 7/18/18.
//

#ifndef SOFTWARE_ORAM_ORAM_H
#define SOFTWARE_ORAM_ORAM_H

#include <vector>
#include <unordered_map>
#include <random>
#include<list>

#include <random>
#include <cmath>
#include <time.h>
#include "Util.h"


//const uint32_t ORAM_Z= 4; // the number of cache lines per bucket
const uint32_t STASH_SIZE = 200; // number of cache lines, the capacity of the stash
const uint32_t STASH_THSHOLD = (STASH_SIZE * 4)/5; // if exceeds, Overflow is imminent.
const unsigned  long long int OFFSET = 0; // offset of address

const bool STASH_VERBOSE = 0;
const int MAX_HISTOGRAM = 40;

const int background_evict_off = 0;
const int EPOC = 1000000;
const int STASH_DRAIN_RATE = 3;
const int STASH_DRAIN_LENGTH = 1000;

//ORAM_Z numbers of nods are in one bucket
typedef struct Node {
    long long int addr;
    int leaf_id;
    bool dummy;
}Node_t;


class ONC {
private:
    std::unordered_map<long long int, std::pair<bool, int>> ONC_umap;
    long long int miss;
    long long int hit;


public:
    ONC(int cached_height);
    ~ONC();
    bool find(long long int tg);
    void insert(long long int tg, int loc);
    void delete_node(long long int tg);
    void write_cache();
    bool Is_empty();
    void hit_update();
    void miss_update();
    void ONC_stat_write();
    long long int get_hit();
    long long int get_miss();
    void reset_hit_miss();
    int get_location(long long int addr);
    void write_ONC();

};





class ORAM : public ONC {
private:
    std::vector<Node_t> buckets;
    std::list<std::pair<int,  long long int>> stash;
    uint32_t * pos_map;
    int hight_cached;
    uint32_t ORAM_Z;
    uint32_t num_leaf;
    uint32_t height;
    unsigned long long int CYCLE_VAL;
    Node_t * sub_buffer;
    long long int histogram[MAX_HISTOGRAM];
    long long int access_BE;
    long long int *leaf_id_count;

    ONC my_ONC;

    bool VERBOSE ;
    bool VERBOSE1 ;
    bool VERBOSE2 ;
    bool VERBOSE_BAEC;
    long long int found_in_stash;
    long long int miss_in_stash;
    unsigned long long int data_cache;
    unsigned long long int data_memory;





public:
    ORAM(const uint32_t num_node, uint32_t zi, int cached_hight,int utilization);
    ~ORAM();


    int access(const char &op, const long long int & addr);
    bool IS_Stash_Overflown();
    bool look_in_stash(const  long long int& addr );
    void write_buckets();
    void write_stash();
    unsigned int get_size_stash();
    void stash_eviction();
    void SnapShot();
    void update_histogram(int level, long long int value);
    void update_histogram_PE(int level, long long int value);
    void set_VERBOSE();
    void reset_VERBOSE();
    long long int get_found_in_stash();
    void write_ONC_stat();
    void write_histogran();
    long long int get_miss_in_stash();
    void background_eviction();
    unsigned long long int get_time();
    unsigned  long long int get_data_cache();
    unsigned  long long int get_data_memory();
    void dummy_data_stat();
    void write_ONC();

    //static size_t rand_int(size_t n);


    //{
    //      return std::uniform_int_distribution<size_t>(0, n - 1)(gen);
    //return rand() % n;
    //  }







};



#endif //SOFTWARE_ORAM_ORAM_H
