//
// Created by meysam on 7/18/18.
//

#include <cmath>
#include "ORAM.h"
#include <random>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <list>
//#include "Util.h"


ONC::ONC(int cached_height){
    long long int size_cache = (1 << cached_height) - 1;
    hit = 0;
    miss = 0;
    /*for (int i=0; i<size_cache; i++)
    {
        ONC_umap.insert(std::make_pair(0,false));
    }*/
}
ONC::~ONC(){
    std::unordered_map<long long int, std::pair< bool, int>>::iterator it;
    for (it= ONC_umap.begin(); it != ONC_umap.end();)
    {
        it = ONC_umap.erase(it);
    }
    ONC_stat_write();
    std::cout << " ONC is distryed\n";

}
void ONC::write_ONC() {
    std::unordered_map<long long int, std::pair< bool, int>>::iterator it;

    for (it = ONC_umap.begin(); it != ONC_umap.end(); it++)
    {
        std::cout << std::hex << it->first  <<" " << std::dec << it->second.first << " " << it->second.second <<"\n";
    }
}
bool ONC::find(long long int tg) {
    auto found = ONC_umap.find(tg);
    if (found == ONC_umap.end())
    {
        //std::cout << "not found\n";
        return (false);
    }
    else
    if (found->second.first == true) {
        //  std::cout << "found\n";
        return (true);
    }
    else {
        //    std::cout << "found invalid\n";
        return (false);
    }
}
void ONC::insert(long long int tg, int loc ) {
    assert(find(tg) == false);
    ONC_umap.insert(std::make_pair(tg, std::make_pair(true, loc)));
}
void ONC::delete_node(long long int tg) {
    //std::cout <<std::hex<< tg << " in the function --> to be deleted..\n";
    if (find(tg)) {
        ONC_umap.erase(tg);
        //  std::cout <<std::hex<< tg << " in the function --> deleted..\n";
    }
    else
        assert(0);
}
void ONC::write_cache() {
    std::unordered_map<long long int, std::pair<bool, int>>::iterator it;
    std::cout << "----------------ONC-------------------\n";
    for (it = ONC_umap.begin(); it != ONC_umap.end(); it++)
        if (it->second.first == true)
            std::cout << it->first <<"\n";
    std::cout << "----------------END-------------------\n";
}

bool ONC::Is_empty() {
    return(ONC_umap.empty());
}

void ONC::hit_update() {
    hit++;
}
void ONC::miss_update() {
    miss++;
}
void ONC::ONC_stat_write() {
    std::cout << "hit miss = " <<std::dec<< hit << " "<< miss <<"\n";
    //output_file << "hit miss = " << hit << " "<< miss <<"\n";
}

long long int ONC::get_hit() {
    return (hit);
}

long long int ONC::get_miss() {
    return (miss);
}

void ONC::reset_hit_miss() {
    miss = 0;
    hit = 0;
}

int ONC::get_location(long long int addr) {

    std::unordered_map<long long int, std::pair< bool, int>>::iterator it;
    for (it= ONC_umap.begin(); it != ONC_umap.end(); it++) {
        if (it->first == addr)
            return (it->second.second);
    }
    assert(0);

}




ORAM::ORAM(const uint32_t num_node, uint32_t zi, int cached_height, int utilization) : ONC(cached_height), my_ONC(cached_height) {
    Node_t Oram_N;
    //srand(time(NULL));
    ORAM_Z = zi;
    height = (uint32_t)ceil(log2((double) num_node));
    hight_cached = cached_height;
    num_leaf = (uint32_t)1 << (height - 1);
    //(uint32_t) (1 << (height-1));
    leaf_id_count = new long long int[num_leaf];
    for (int i=0; i<num_leaf; i++)
    {
        leaf_id_count[i] = 0;
    }
    for (int i=0; i<height; i++)
    {
        histogram[i] = 0;
    }



    // to generate random number
    //std::random_device rd;  //Will be used to obtain a seed for the random number engine

    //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    //std::uniform_int_distribution<> rnd(0, num_leaf-1);
    pos_map = new uint32_t [num_node*ORAM_Z];
    long long int num_block = 0;
    for (uint32_t i =0; i< ((1 << height) - 1); i++)
    {
        for (uint32_t j=0; j<ORAM_Z; j++)
        {
            Oram_N.dummy = 1; // it is a dummy node
            Oram_N.addr = -1; //OFFSET + ((i * ORAM_Z + j)<< 6); // cache line address
            Oram_N.leaf_id = -1; //rnd(gen);
            pos_map[(i * ORAM_Z + j)]= -1; //Oram_N.leaf_id; //fill the posmap
            buckets.push_back(Oram_N); // insert the node into the ORAM_node vector
            num_block++;

        }
    }
    // initilize the pos-map with random number for entire memory
    for (int i =0; i<(num_node*ORAM_Z* utilization/100); i++ )
    {
        pos_map[i] = Util::rand_int(num_leaf);
    }

    VERBOSE = 0;
    VERBOSE1 = 0;
    VERBOSE2 = 0;
    VERBOSE_BAEC = 0;

    found_in_stash = 0;
    miss_in_stash = 0;
    CYCLE_VAL = 0;
    data_cache = 0;
    data_memory = 0;
    access_BE = 0;

    std::cout << "ORAM and its pos map are initilized \n";
    std::cout << "height ............................ " << height <<"\n";
    std::cout << "treetop level ..................... " << hight_cached <<"\n";
    std::cout << "num of leaf ....................... " << num_leaf << "\n";
    std::cout << "num of blocks ..................... " << num_block << "\n";
    std::cout << "ORAM_Z ............................ " << ORAM_Z << "\n";
    std::cout << "STASH_SIZE ........................ " << STASH_SIZE << "\n";
    std::cout << "STASH_THSHOLD ..................... " << STASH_THSHOLD << "\n";
    std::cout << "OFFSET ............................ " << OFFSET << "\n";


}

ORAM::~ORAM() {

    //write_buckets();
    buckets.clear();
    free (pos_map);
    write_histogran();


    std::cout<< "hit = : "<< my_ONC.get_hit() << "\n";
    std::cout<< "miss = : "<< my_ONC.get_miss() << "\n";

    std::cout << "found_in_stash = " <<  found_in_stash << "\n";
    std::cout << "missed_in_stash = " <<  miss_in_stash << "\n";

    std::cout << "access BE = " <<  access_BE << "\n";
    // draw the histogram for count
    long long int max, min;
    max = 0;
    min = 1000000;
    for (int i=0; i<num_leaf; i++)
    {
        if (max < leaf_id_count[i])
            max = leaf_id_count[i];
        if (leaf_id_count[i] < min)
            min = leaf_id_count[i];
    }
    const int lev = 100;
    int length_levels = (max - min) / lev;
    /*int hist[lev];
    for (int i=0; i<lev; i++) hist[i] = 0;
    for (int i = 0; i<num_leaf; i++)
    {
        hist[((leaf_id_count[i] - min)/length_levels)]++;
    }*/
    //  for (int i=0; i<num_leaf; i++)
    //      std::cout << "leaf_id_count[" << i << " ] = " << leaf_id_count[i] <<"\n";


    /*for (int i=0; i<lev; i++)
    {
        std::cout << "hist[" << i << " ] = " << hist[i] <<"\n";
    }*/



    std::cout << "ORAM destroyed memory back\n";

}

void ORAM::set_VERBOSE() {
    VERBOSE = 1;
    VERBOSE1 = 1;
    VERBOSE2 = 1;

}


void ORAM::reset_VERBOSE() {
    VERBOSE = 0;
    VERBOSE1 = 0;
    VERBOSE2 = 0;

}

void ORAM::write_ONC_stat() {
    my_ONC.ONC_stat_write();
}




bool ORAM::look_in_stash(const  long long int &addr) {
    bool found = 0;
    std::list<std::pair<int,  long long int>>::iterator  it;
    for (it=stash.begin(); it!=stash.end(); it++)
    {
        if (it->second == addr) {
            found = 1;
            break;
        }

    }
    return (found);

}

long long int ORAM::get_found_in_stash() {
    return (found_in_stash);
}

void ORAM::write_buckets() {
    std::vector<Node_t>::iterator it;
    long long int num = 0;

    for (it=buckets.begin(); it!=buckets.end(); ++it)
    {
        if (it->dummy == 0)
            std::cout <<std::dec<< num << "  [ " <<std::dec<< it->dummy << " " <<std::dec<< it->leaf_id << " " <<std::hex<< it->addr << " ]" << "\n";
        num++;
    }
}

void ORAM::background_eviction() {
    // VERBOSE_BAEC = 0;
    int lf_id = -1;
    //while (lf_id == -1)
    //lf_id = pos_map[rand()%1024][rand()%((ORAM_Z*((1<<(height))-1))/1024)];
    lf_id = Util::rand_int(num_leaf);//rand()%num_leaf;
    //assert(lf_id != -1);
    //if (VERBOSE_BAEC)
    //{
    std::cerr << " ==========================  BackGraound Eviction is undergone for ld_id = " << std::dec << lf_id << " stash size = " << std::dec << get_size_stash() <<"\n";
    //}
    access_BE++;
    int Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
    int level_node_in_path = height - 1;

    while (Node_in_Path >= 0) {
        CYCLE_VAL++;

        if (VERBOSE_BAEC)
            std::cout << "Node_in_Path= " <<std::dec<< Node_in_Path << "\n";
        //starting from leaf go up to the root
        for (uint32_t j = 0; j < ORAM_Z; j++) {
            if (VERBOSE_BAEC)
                std::cout << "j= " << j << " bucket[ " <<std::dec<< Node_in_Path * ORAM_Z + j << " ] " <<std::hex<<  buckets[Node_in_Path * ORAM_Z + j].addr;

            if (level_node_in_path < hight_cached)
            {
                if (buckets[Node_in_Path * ORAM_Z + j].dummy == 0)
                    my_ONC.delete_node(buckets[Node_in_Path * ORAM_Z + j].addr);
            }

            // this is not our target either it is a dummy throw it a way!
            // or it is a data put it in the stash
            if (buckets[Node_in_Path * ORAM_Z + j].dummy != 1) {
                // data
                if (VERBOSE_BAEC) {
                    std::cout << "--> not target but data taken to stash \n";

                }
                if (VERBOSE_BAEC)
                    std::cout <<std::dec<<"bucket[" <<Node_in_Path * ORAM_Z + j <<" ] "<< std::hex<<buckets[Node_in_Path * ORAM_Z + j].addr<<"\n";
                if (look_in_stash( buckets[Node_in_Path * ORAM_Z + j].addr)) {
                    std::cout << " --> addr = " <<  buckets[Node_in_Path * ORAM_Z + j].addr << "\n";
                    // write_stash();
                }
                assert(!look_in_stash(buckets[Node_in_Path * ORAM_Z + j].addr));
                stash.push_back(std::make_pair(buckets[Node_in_Path * ORAM_Z + j].leaf_id,buckets[Node_in_Path * ORAM_Z + j].addr));


                //if (VERBOSE_BAEC)
                //    std::cout << "+++stash size is " << stash[leaf_id_2_stash_id(buckets[Node_in_Path * ORAM_Z + j].leaf_id)].size()<< "\n";
            } else {
                // dummy
                if (VERBOSE_BAEC)
                    std::cout << "--> dummy \n";
            }
            // whether or not it is dummy, it should be converted to a  dummy one because it is evacuated into the stash or deleted
            buckets[Node_in_Path * ORAM_Z + j].dummy = 1;
            buckets[Node_in_Path * ORAM_Z + j].addr = -1;
            buckets[Node_in_Path * ORAM_Z + j].leaf_id = -1;

        }// for cache lines per bucket
        if (Node_in_Path == 0) break;
        Node_in_Path = (Node_in_Path - 1) / 2;
        level_node_in_path--;
        // go to the parent node up to the root
    } // while bucket

    if (VERBOSE_BAEC)
        std::cout << " Back_Evict: eviction from stash back to ORAM starting .....\n";

    Node_in_Path = (1 << (height - 1)) - 1 + lf_id;

    int level = height - 1; // the maximum level is for leaves in the cache

    std::list<std::pair<int,  long long int>>::iterator it_j;

    int lv = 0;  //the max is for the root

    while (Node_in_Path >= 0)
    {
        CYCLE_VAL++;
        // find the correct cachelines belonging to correct leaf that can be located in this bucket
        if(VERBOSE2)
        {
            std::cout << "Node_in path = " << Node_in_Path << "\n";
        }

        int num_cache_per_bucket = 0;
        it_j = stash.begin();
        for (; it_j != stash.end();){
            //  std::cout << std::dec << it_j->first <<" "<< lf_id <<" " << (it_j->first ^ lf_id) <<" " << ((1<< lv)-1)<<"\n";
            //  std::cout << std::dec << it_j->first <<" " <<std::hex<<" " << it_j->second <<"\n";


            if ((it_j->first ^ lf_id) <= ((1<< lv)-1))
            {
                if(VERBOSE2)
                    std::cout <<std::dec << it_j->first <<" "<< lf_id <<" " << lv <<" "<< ((1<< lv)-1) <<"\n";
                if(VERBOSE2)
                {
                    std::cout << "buckets[ " <<std::dec<< Node_in_Path*ORAM_Z+num_cache_per_bucket <<" ]" << "is updated ";
                }

                // this is good for this node
                buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].dummy = 0;
                buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].leaf_id = it_j->first;
                buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].addr = it_j->second;
                if(VERBOSE2)
                {
                    std::cout << "with " <<std::dec << it_j->first << std::hex <<"  "<< it_j->second << "\n";
                }
                if (level < hight_cached)
                {
                    my_ONC.insert(buckets[Node_in_Path * ORAM_Z + num_cache_per_bucket].addr, Node_in_Path * ORAM_Z + num_cache_per_bucket);
                    //std::cout << std::dec<< "bucket[" << Node_in_Path * ORAM_Z + num_cache_per_bucket << "] = " <<std::hex <<
                    //          buckets[Node_in_Path * ORAM_Z + num_cache_per_bucket].addr << " inserted in cache \n";
                }
                it_j = stash.erase(it_j);
                num_cache_per_bucket++;
                if (num_cache_per_bucket >= ORAM_Z)
                    break;
            }
            else ++it_j;
        }
        //std::cout <<"===========\n";
        //getchar();

        if (num_cache_per_bucket < ORAM_Z)
        {
            // not found enough fit canidates go for dummies!
            for (int j=(num_cache_per_bucket); j<ORAM_Z; j++)
            {
                Node_t tmp;
                tmp.dummy = 1;
                tmp.addr = -2;
                tmp.leaf_id = -2;
                buckets[Node_in_Path * ORAM_Z + j] = tmp;
                if(VERBOSE2)
                {
                    std::cout << "bucket[ " <<std::dec<<  Node_in_Path * ORAM_Z + j << " ] is dummy " << "\n";
                }
            }
        }


        /*           for (int i=0; i<(1<<lv); i++)
                   {


                       int leaf_look = ((1<<lv)*(Node_in_Path - ((1<< level)-1)))+i;
                       // go to the stash and find the cacheline with this leaf_look
                       if(VERBOSE)
                       {
                           std::cout << "leaf_look = " << leaf_look << "\n";
                       }

                       it_j = stash.begin();
                       for (; it_j != stash.end();)
                       {
                           if (it_j->first  == leaf_look)
                           {
                               if(VERBOSE)
                               {
                                   std::cout << "buckets[ " << Node_in_Path*ORAM_Z+num_cache_per_bucket <<" ]" << "is updated ";
                               }
                               // this cacheline should be evicted from stash to this node.
                               buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].dummy = 0;
                               buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].leaf_id = leaf_look;
                               buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].addr = it_j->second;
                               if(VERBOSE)
                               {
                                   std::cout << "with " <<std::dec << leaf_look << std::hex <<"  "<< it_j->second << "\n";
                               }
                              // std::cout << " before erase ...\n";
                               stash.erase(it_j);
                               if (STASH_VERBOSE)
                                   std::cout << "+++stash size is " << stash.size()<< "\n";
                              // std::cout << " After erase ...\n";
                               num_cache_per_bucket++;
                               if (num_cache_per_bucket == ORAM_Z)
                               {
                                   // ths node is full
                                   break;
                               }

                           }
                           else ++it_j;
                       }
                       if (num_cache_per_bucket ==ORAM_Z )
                           break;

                       if ((i == (1<<lv)-1) && (num_cache_per_bucket < ORAM_Z))
                       {
                           if(VERBOSE)
                           {
                               std::cout << "dummy --> num_cache_per_bucket = " << num_cache_per_bucket << "\n";
                           }
                           // we could not find enough non-dummy cache lines go for dummies
                           if (VERBOSE)
                                std::cout << "num_cache_per_bucket = " << num_cache_per_bucket << "\n";
                           for (int j=(num_cache_per_bucket); j<ORAM_Z; j++)
                           {
                               Node_t tmp;
                               tmp.dummy = 1;
                               tmp.addr = -2;
                               tmp.leaf_id = -2;
                               buckets[Node_in_Path * ORAM_Z + j] = tmp;
                               if(VERBOSE)
                               {
                                   std::cout << "bucket[ " <<  Node_in_Path * ORAM_Z + j << " ] is dummy " << "\n";
                               }
                           }
                       }
                   }*/
        if (Node_in_Path == 0) break;
        Node_in_Path = (Node_in_Path - 1) / 2;
        lv++;
        level--;
    }

    /*  while (Node_in_Path >= 0)
      {
          CYCLE_VAL++;
          // find the correct cachelines belonging to correct leaf that can be located in this bucket
          if(VERBOSE_BAEC)
          {
              std::cout << "Node_in path = " <<std::dec<< Node_in_Path << "\n";
          }
          int num_per_buck = 0;

          int num_group = (1<<(height-level-1)); // howmany leaf id may fit in to this level of lf_id
          int num_cache_per_bucket = 0;
          for (int i=0; i<num_group; i++)
          {

              // this leaf can be fit into this level of lf_id
              int cand_leaf_id = num_group * (lf_id / num_group) + i;

              std::list<std::pair<int, long long int>>::iterator it_st;


              //for (int j = 0; j < ORAM_Z; j++) Z_top[j].top_pair.second.first = -1;
              //int min_index = 0; // index of minimum in Z_top
              //num_cache_per_bucket = 0;
              //if (VERBOSE_BAEC)
              //    write_a_stash(leaf_id_2_stash_id(cand_leaf_id));


              for (it_st = stash.begin(); it_st != stash.end();) {
                  if (it_st->first == cand_leaf_id) {
                      if (VERBOSE_BAEC) {
                          std::cout << "buckets[ " << std::dec << Node_in_Path * ORAM_Z + num_per_buck << " ]"
                                    << "is updated ";
                      }
                      buckets[Node_in_Path * ORAM_Z + num_per_buck].dummy = 0;
                      buckets[Node_in_Path * ORAM_Z + num_per_buck].leaf_id = it_st->first;
                      buckets[Node_in_Path * ORAM_Z + num_per_buck].addr = it_st->second;

                      if (level < hight_cached)
                          my_ONC.insert(buckets[Node_in_Path * ORAM_Z + num_per_buck].addr);

                      if (VERBOSE_BAEC) {
                          std::cout << "with " << std::dec << it_st->first << std::hex << "  " << it_st->second
                                    << "  from memory" << "\n";
                      }
                      // delete from stash
                      it_st = stash.erase(it_st);
                      num_per_buck++;
                      if (num_per_buck == ORAM_Z)
                          break;
                  }
                  else
                      it_st++;

              }
              if (num_per_buck == ORAM_Z)
                  break;

          }
          if (num_per_buck <ORAM_Z)
          {
              // go for dummy
              for (int ii=num_per_buck; ii<ORAM_Z;ii++)
              {
                  buckets[Node_in_Path * ORAM_Z + ii].dummy = 1;
                  buckets[Node_in_Path * ORAM_Z + ii].leaf_id = -1;
                  buckets[Node_in_Path * ORAM_Z + ii].addr = -1;

                  if(VERBOSE_BAEC)
                  {
                      std::cout << "bucket[ " <<std::dec<<  Node_in_Path * ORAM_Z + ii << " ] is dummy " << "\n";
                  }

              }
          }

          if (Node_in_Path == 0) break;
          Node_in_Path = (Node_in_Path - 1) / 2;
          // lv++;
          level--;
      } // while bucket for eviction
  */
    // if (VERBOSE_BAEC)
    std::cerr << " Back eviction has endded !!!! ====================== " << " stash size = " << std::dec <<get_size_stash() <<"\n";

}


void ORAM::write_ONC() {
    my_ONC.write_ONC();

}

int ORAM::access(const char &op, const long long int &addr) {

    if (VERBOSE) {
        if (op == 'r')
            std::cout << "addr = " <<std::hex<< addr << "  =========ORAM access started for a read!======================\n";
        else
            std::cout << "addr = " <<std::hex<< addr << "  =========ORAM access started for a write!======================\n";
    }
    // determine whether or not it already exists in stash
    if (look_in_stash(addr)) {
        found_in_stash++;
        return 0;

    }

    miss_in_stash++;

    if (!my_ONC.Is_empty())
    {
        if (my_ONC.find(addr) == true)
        {
            my_ONC.hit_update();
            return 0;
        }
        else
            my_ONC.miss_update();
    }
    else
        my_ONC.miss_update();

    // determine which path (leaf_id) should be fetched
    int lf_id = pos_map[(addr - OFFSET)/64];
    /*if (op == 'r')
    {
        //assert(pos_map[(addr - OFFSET) / 64] != -1);
        lf_id = pos_map[(addr - OFFSET) / 64];
        if (VERBOSE)
            std::cout << "lf_id = " <<std::dec<< lf_id << "\n";

    }
    else
    {
        //std::cout <<"wrrrrrrrrite\n";
        //write
        if (pos_map[(addr - OFFSET) / 64] == -1)
        {
            // this is the first time to write
            pos_map[(addr - OFFSET) / 64] = Util::rand_int(num_leaf); //rand()%num_leaf; //rnd(gen);
            leaf_id_count[pos_map[(addr - OFFSET) / 64]]++;
            lf_id = pos_map[(addr - OFFSET) / 64];
            // std::cout <<"rand_int = " << lf_id <<"\n";
            // getchar();

        }
        else {
            lf_id = pos_map[(addr - OFFSET) / 64];
        }
    }*/

    std::vector<std::pair<int,  long long int>> stash_buffer; // to tore candidates to evict there is a triple (int,int, long long int)
    // to generate random number

    // fetch all the blocks along with the path
    int Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
    if (VERBOSE) {
        std::cout << "lf_id = " <<std::dec<< lf_id << " Node_in_Path " << Node_in_Path << "\n";
    }
    bool found = 0;
    int level_node_in_path = height - 1;


    while (Node_in_Path >= 0)
    {
        CYCLE_VAL++;
        if (VERBOSE)
            std::cout << "Node_in_Path= " <<std::dec<< Node_in_Path << "\n";
        //starting from leaf go up to the root
        for (uint32_t j = 0; j < ORAM_Z; j++)
        {
            if (VERBOSE)
                std::cout << "j= " << j << " bucket[ " <<std::dec<< Node_in_Path * ORAM_Z + j << " ] " <<std::hex<<  buckets[Node_in_Path * ORAM_Z + j].addr;

            if (level_node_in_path < hight_cached) // if it exists in treetop
            {
                if (buckets[Node_in_Path * ORAM_Z + j].dummy == 0)
                {
                    // std::cout << "\nbuckets[" <<std::dec<< Node_in_Path * ORAM_Z + j << "]" << "deleted from the cache\n";
                    my_ONC.delete_node(buckets[Node_in_Path * ORAM_Z + j].addr);
                }
            }
            if (buckets[Node_in_Path * ORAM_Z + j].addr == addr)
            {
                found =1;
                //update_histogram(level_node_in_path, histogram[level_node_in_path]+ 1 );
                histogram[level_node_in_path]++;
                // if it is a read it should be valid
                if (op == 'r')
                    assert(buckets[Node_in_Path * ORAM_Z + j].dummy != 1);
                // this is our target
                // change the leaf id and put in the stash
                int new_leaf_id = Util::rand_int(num_leaf); //rand()%num_leaf; //rnd(gen);
                leaf_id_count[new_leaf_id]++;
                stash.push_back(std::make_pair(new_leaf_id,addr)); // (new_leaf_id, addr)
                if (VERBOSE)
                    std::cout << "\n stash size.................. " <<std::dec<< get_size_stash() <<"\n";
                pos_map[(addr - OFFSET) / 64] = new_leaf_id;
                if (VERBOSE)
                    std::cout << "--> our target taken to stash " << " new leaf = " << new_leaf_id << "\n";
            } else
                {
                // this is not our target either it is a dummy throw it a way!
                // or it is a data put it in the stash
                    if (buckets[Node_in_Path * ORAM_Z + j].dummy != 1)
                    {
                        // data
                        if (VERBOSE) {
                            std::cout << "--> not target but data taken to stash \n";

                        }
                        // its exact node will go to the stash because it is not
                        // our target while it is not dummy either
                        stash.push_back(std::make_pair(buckets[Node_in_Path * ORAM_Z + j].leaf_id, buckets[Node_in_Path * ORAM_Z + j].addr ));
                        if (STASH_VERBOSE)
                            std::cout << "+++stash size is " << stash.size()<< "\n";
                        if (VERBOSE){
                            std::cout << "stash size.................. " <<std::dec<< get_size_stash() <<"\n";
                            write_stash();
                        }
                    }
                    else
                    {
                        // dummy
                        if (VERBOSE)
                            std::cout << "--> dummy \n";
                    }
                }
            // whether or not it is dummy, it should be converted to a  dummy one because it is evacuated into the stash or deleted
            buckets[Node_in_Path * ORAM_Z + j].dummy = 1;
            buckets[Node_in_Path * ORAM_Z + j].addr = -1;
            buckets[Node_in_Path * ORAM_Z + j].leaf_id = -1;
        }// for cache lines per bucket
        if (Node_in_Path == 0) break;
        Node_in_Path = (Node_in_Path - 1) / 2;
        level_node_in_path--;
        // go to the parent node up to the root
    } // while bucket

    if (op == 'r')
    {
        if (found == 0)
            std::cout << " address " <<std::hex<< addr << " not found in read!\n";
        assert(found == 1);
    }
    else
    {
        if (found == 0)
        {
            // this was a write into a new address
            // we should change its leaf id and store in the stash
            int new_leaf_id = Util::rand_int(num_leaf); //rand()%num_leaf; //rnd(gen);
            leaf_id_count[new_leaf_id]++;
            pos_map[(addr-OFFSET)/64] = new_leaf_id;
            stash.push_back(std::make_pair(new_leaf_id, addr));
            if (STASH_VERBOSE)
                std::cout << "+++stash size is " << stash.size()<< "\n";
            if (VERBOSE)
                std::cout << " addr = " << addr << " with leaf id = " <<std::dec << new_leaf_id << "--> stash \n";

        }

    }


    // evict appropriate cache lines back to the tree in a same path
    if (VERBOSE2)
        std::cout << " eviction from stash back to ORAM starting .....\n";
    int found_to_evict = 0;
    std::list<std::pair<int,  long long int>>::iterator it_j;
    Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
    if (VERBOSE2)
        std::cout << " Node_in_Path = " <<std::dec<<Node_in_Path<<" lf_id= "<< lf_id<<"\n";
    int lv = 0;  //the max is for the root
    int level = height - 1; // the maximum level is for leaves

    if (!stash.empty()) {
        if (VERBOSE2) {
            std::cout << "stash is NOT empty ..\n";
            write_stash();
        }

        // this should happen for every node (bucket)
        while (Node_in_Path >= 0)
        {
            CYCLE_VAL++;
            // find the correct cachelines belonging to correct leaf that can be located in this bucket
            if(VERBOSE2)
            {
                std::cout << "Node_in path = " <<std::dec<< Node_in_Path << "\n";
            }

            int num_cache_per_bucket = 0;
            it_j = stash.begin();
            for (; it_j != stash.end();)
            {
                if ((it_j->first ^ lf_id) <= ((1<< lv)-1))
                    //if ((it_j->first >> lv) == (lf_id >> lv))
                {
                    if(VERBOSE2)
                        std::cout <<std::dec << it_j->first <<" "<< lf_id <<" " << lv <<" "<< ((1<< lv)-1) <<"\n";
                    if(VERBOSE2)
                    {
                        std::cout << "buckets[ " <<std::dec<< Node_in_Path*ORAM_Z+num_cache_per_bucket <<" ]" << "is updated ";
                    }

                    // this is good for this node
                    buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].dummy = 0;
                    buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].leaf_id = it_j->first;
                    buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].addr = it_j->second;
                    if(VERBOSE2)
                    {
                        std::cout << "with " <<std::dec << it_j->first << std::hex <<"  "<< it_j->second << "\n";
                    }
                    if (level < hight_cached)
                    {
                        my_ONC.insert(buckets[Node_in_Path * ORAM_Z + num_cache_per_bucket].addr, Node_in_Path * ORAM_Z + num_cache_per_bucket);
                        //std::cout << std::dec<< "bucket[" << Node_in_Path * ORAM_Z + num_cache_per_bucket << "] = " <<std::hex <<
                        //          buckets[Node_in_Path * ORAM_Z + num_cache_per_bucket].addr << " inserted in cache \n";
                    }
                    it_j = stash.erase(it_j);
                    num_cache_per_bucket++;
                    if (num_cache_per_bucket >= ORAM_Z)
                        break;
                }
                else ++it_j;
            }
            if (num_cache_per_bucket < ORAM_Z)
            {
                // not found enough fit canidates go for dummies!
                for (int j=(num_cache_per_bucket); j<ORAM_Z; j++)
                {
                    Node_t tmp;
                    tmp.dummy = 1;
                    tmp.addr = -1;
                    tmp.leaf_id = -1;
                    buckets[Node_in_Path * ORAM_Z + j] = tmp;
                    if(VERBOSE2)
                    {
                        std::cout << "bucket[ " <<std::dec<<  Node_in_Path * ORAM_Z + j << " ] is dummy " << "\n";
                    }
                }
            }


            /*           for (int i=0; i<(1<<lv); i++)
                       {


                           int leaf_look = ((1<<lv)*(Node_in_Path - ((1<< level)-1)))+i;
                           // go to the stash and find the cacheline with this leaf_look
                           if(VERBOSE)
                           {
                               std::cout << "leaf_look = " << leaf_look << "\n";
                           }

                           it_j = stash.begin();
                           for (; it_j != stash.end();)
                           {
                               if (it_j->first  == leaf_look)
                               {
                                   if(VERBOSE)
                                   {
                                       std::cout << "buckets[ " << Node_in_Path*ORAM_Z+num_cache_per_bucket <<" ]" << "is updated ";
                                   }
                                   // this cacheline should be evicted from stash to this node.
                                   buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].dummy = 0;
                                   buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].leaf_id = leaf_look;
                                   buckets[Node_in_Path*ORAM_Z+num_cache_per_bucket].addr = it_j->second;
                                   if(VERBOSE)
                                   {
                                       std::cout << "with " <<std::dec << leaf_look << std::hex <<"  "<< it_j->second << "\n";
                                   }
                                  // std::cout << " before erase ...\n";
                                   stash.erase(it_j);
                                   if (STASH_VERBOSE)
                                       std::cout << "+++stash size is " << stash.size()<< "\n";
                                  // std::cout << " After erase ...\n";
                                   num_cache_per_bucket++;
                                   if (num_cache_per_bucket == ORAM_Z)
                                   {
                                       // ths node is full
                                       break;
                                   }

                               }
                               else ++it_j;
                           }
                           if (num_cache_per_bucket ==ORAM_Z )
                               break;

                           if ((i == (1<<lv)-1) && (num_cache_per_bucket < ORAM_Z))
                           {
                               if(VERBOSE)
                               {
                                   std::cout << "dummy --> num_cache_per_bucket = " << num_cache_per_bucket << "\n";
                               }
                               // we could not find enough non-dummy cache lines go for dummies
                               if (VERBOSE)
                                    std::cout << "num_cache_per_bucket = " << num_cache_per_bucket << "\n";
                               for (int j=(num_cache_per_bucket); j<ORAM_Z; j++)
                               {
                                   Node_t tmp;
                                   tmp.dummy = 1;
                                   tmp.addr = -2;
                                   tmp.leaf_id = -2;
                                   buckets[Node_in_Path * ORAM_Z + j] = tmp;
                                   if(VERBOSE)
                                   {
                                       std::cout << "bucket[ " <<  Node_in_Path * ORAM_Z + j << " ] is dummy " << "\n";
                                   }
                               }
                           }
                       }*/
            if (Node_in_Path == 0) break;
            Node_in_Path = (Node_in_Path - 1) / 2;
            lv++;
            level--;
        }

        /*int Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
        int leveling = height - 1;


        data_cache = 0; data_memory = 0;
        while (Node_in_Path >= 0)
        {
            for (int iii=0; iii<ORAM_Z; iii++)
            {
                if (buckets[Node_in_Path * ORAM_Z + iii].dummy == 0)
                {
                    if (leveling < hight_cached)
                        data_cache++;
                    else
                        data_memory++;

                }

            }

            if (Node_in_Path == 0) break;
            Node_in_Path = (Node_in_Path - 1) / 2;
            leveling--;
        }*/

        /*it_j = stash.begin();
        while (it_j != stash.end()) {
            if (VERBOSE) {
                std::cout << "[ " << it_j->second << " ," << it_j->first << " ]" << " is in  the stash\n";
            }
            if (it_j->first == lf_id) {
                // this is a good candate for eviction
                if (VERBOSE)
                    std::cout << "inserting into stash_buffer  "<<std::dec<< it_j->first <<" "<<std::hex<< it_j->second <<"\n";
                stash_buffer.push_back(std::make_pair(it_j->first, it_j->second));
                 //std::pair<int,  long long int> hh = stash_buffer.back();
                 //std::cout <<"  " << hh.first << " " <<hh.second;
                stash.erase(it_j);
                //std::cout << "stash size.................. " << get_size_stash() <<"\n";
                if (VERBOSE) {
                    std::cout << it_j->second << " added to stash_buufer\n";
                }
            }
            it_j++;
            ///std::cout << "it_j added \n";

        }*/
        ///std::cout << " Exit from while ..  " << stash_buffer.size() <<"\n";
    } else if (VERBOSE)
        std::cout << "stash is empty ..\n";
    /*int num_candidates;
    if (stash_buffer.empty()) {
        num_candidates = 0;
        if (VERBOSE)
            std::cout << "stash_buffer is empty ..\n";
    } else
        num_candidates = stash_buffer.size();

    Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
    Node_t tmp;
    std::pair<int,  long long int> pair_tmp;

    while (num_candidates > 0 && Node_in_Path >= 0) {
        if (VERBOSE) {
            std::cout << "num_candidates= " << num_candidates << " Node_in_Path= " << Node_in_Path << "\n";
        }

        pair_tmp = stash_buffer.back(); // the last elements stored
        stash_buffer.pop_back(); // and deleted
        for (uint32_t j = 0; j < ORAM_Z; j++) {
            tmp.dummy = 0; // data
            tmp.addr = pair_tmp.second + j * 64;
            tmp.leaf_id = pair_tmp.first;
            if (VERBOSE)
                std::cout << "j= " << std::dec << j <<" "<<  tmp.addr <<" " << tmp.leaf_id << "\n";
            buckets[Node_in_Path * ORAM_Z + j] = tmp;
            num_candidates--;
            if (VERBOSE) {
                std::cout << "REAL bucket[  " << Node_in_Path * ORAM_Z + j << " ]" << "\n";
                std::cout << buckets[Node_in_Path * ORAM_Z + j].addr << " " << buckets[Node_in_Path * ORAM_Z + j].dummy << " " <<buckets[Node_in_Path * ORAM_Z + j].leaf_id <<"\n";
            }

        }
        if (Node_in_Path == 0) break;
        Node_in_Path = (Node_in_Path - 1) / 2; // go to its parents along with the path up to the root
    }*/

    //write_buckets();


    // we need to fill the rest of buckets with dummies!
    /*if (Node_in_Path != 0)
    {
        while (Node_in_Path >= 0) {
            if (VERBOSE) {
                std::cout << "Node_in_Path= " << Node_in_Path << "\n";
            }
            for (uint32_t j = 0; j < ORAM_Z; j++) {
                tmp.dummy = 1; // dummy
                tmp.addr = -2;
                tmp.leaf_id = -2;
                if (VERBOSE) {
                    std::cout << "bucket[  " << Node_in_Path * ORAM_Z + j << " ]" << "\n";
                }
                buckets[Node_in_Path * ORAM_Z + j] = tmp;

            }
            if (Node_in_Path == 0) break;
            Node_in_Path = (Node_in_Path - 1) / 2; // go to its parents along with the path up to the root


        } //while

    }//if
    */

}

bool ORAM::IS_Stash_Overflown() {
    return (stash.size() > STASH_SIZE);
}

unsigned int ORAM::get_size_stash(){
    return (stash.size());
}

void ORAM::write_stash(){
    std::list<std::pair<int,  long long int>>::iterator  it;
    std::cout  <<("Stash------------------------------\n");
    for (it=stash.begin(); it !=stash.end(); it++)
    {
        std::cout <<std::dec<< it->first <<" " <<std::hex<< it->second <<"\n";
    }

}

void ORAM::SnapShot() {
    long long int num_block = 0;
    std::vector<Node_t>::iterator it;

    for (it=buckets.begin(); it != buckets.end(); it++)
    {
        std::cout << " BUCKET " << it->dummy << " " <<std::dec<< it->leaf_id << " " <<std::hex<< it->addr << "\n";
        num_block++;

    }

    std::list<std::pair<int,  long long int>>::iterator it_stash;

    for (it_stash= stash.begin(); it_stash!=stash.end(); it_stash++)
    {
        std::cout << " STASH " << std::dec<< it_stash->first <<" " <<std::hex<< it_stash->second << "\n";
    }

    for (long long int i=0; i<num_block; i++)
    {
        std::cout << " POSMAP " <<std::dec<< pos_map[i] << "\n";
    }





    std::cout << "ORAM and its pos map are initilized \n";
    std::cout << "height ............................ " << height <<"\n";
    std::cout << "num of leaf ....................... " << num_leaf << "\n";
    std::cout << "num of blocks ..................... " << num_block << "\n";
    std::cout << "ORAM_Z ............................ " << ORAM_Z << "\n";
    std::cout << "STASH_SIZE ........................ " << STASH_SIZE << "\n";
    std::cout << "STASH_THSHOLD ..................... " << STASH_THSHOLD << "\n";
    std::cout << "OFFSET ............................ " << OFFSET << "\n";

}
/*
void ORAM::stash_eviction() {

    std::cout << "stash eviction is startin ....\n";

    std::unordered_map<int,  long long int>::iterator  it;
    std::unordered_map<int,  int> hm;
    std::unordered_map<int,  int>::iterator  it_hm;
    std::cout << " stash size: " << get_size_stash() <<"\n";
    // build a histogram
    for (it=stash.begin(); it != stash.end(); it++)
    {
        hm[it->first]++;
    }

    int max_hm = 0;
    int i_max;
    int ntimes = 0;
    // find the most beneficial leaf in ORAM to evict
    for (it_hm=hm.begin(); it_hm!=hm.end();it_hm++)
    {
      //  if (it_hm->second != 1)
        std::cout << it_hm->first << " " << it_hm->second << "\n";

        if (it_hm->second >= max_hm)
        {
            max_hm = it_hm->second;
            i_max = it_hm->first;
        }

    }

    std::cout << "i_max = " << i_max;

    // access to leaf = i_max and
    for (it=stash.begin(); it != stash.end(); it++)
    {
        if (it->first == i_max)
        {
            std::cout << "--->" << it->first << " " << it->second << "\n";
            //access('w', it->second);
            // here we should read this path and in write time nmake sure that
            // the one which is in the stash will be written back

            int lf_id = i_max; // read / write this path


            std::vector<std::pair<int,  long long int>> stash_buffer; // to tore candidates to evict
            // to generate random number

            // fetch all the blocks along with the path
            int Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
            if (VERBOSE) {
                std::cout << "lf_id = " << lf_id << " Node_in_Path " << Node_in_Path << "\n";
            }
            bool found = 0;

            while (Node_in_Path >= 0) {
                if (VERBOSE)
                    std::cout << "Node_in_Path= " << Node_in_Path << "\n";
                //starting from leaf go up to the root
                for (uint32_t j = 0; j < ORAM_Z; j++) {
                    if (VERBOSE)
                        std::cout << "j= " << j << " bucket[ " << Node_in_Path * ORAM_Z + j << " ] " <<  buckets[Node_in_Path * ORAM_Z + j].addr;


                    // this is not our target either it is a dummy throw it a way!
                    // or it is a data put it in the stash
                    if (buckets[Node_in_Path * ORAM_Z + j].dummy != 1) {
                        // data
                        if (VERBOSE) {
                            std::cout << "--> not target but data taken to stash \n";
                            write_stash();
                        }
                        // its exact node will go to the stash because it is not
                        // our target while it is not dummy either
                        stash[buckets[Node_in_Path * ORAM_Z + j].leaf_id] = buckets[Node_in_Path * ORAM_Z + j].addr;
                        if (VERBOSE)
                            std::cout << "stash size.................. " << get_size_stash() <<"\n";
                    } else {
                        // dummy
                        if (VERBOSE)
                            std::cout << "--> dummy \n";
                    }




                    // whether or not it is dummy, it should be converted to a  dummy one because it is evacuated into the stash or deleted
                   ////////// buckets[Node_in_Path * ORAM_Z + j].dummy = 1;
                }// for cache lines per bucket
                ///////////////if (Node_in_Path == 0) break;
                Node_in_Path = (Node_in_Path - 1) / 2;
                // go to the parent node up to the root
            } // while bucket




            // evict appropriate cache lines back to the tree in a same path
            if (VERBOSE)
                std::cout << " eviction from stash back to ORAM starting .....\n";
            int found_to_evict = 0;
            std::unordered_map<int,  long long int>::iterator it_j;

            if (!stash.empty()) {
                if (VERBOSE) {
                    std::cout << "stash is NOT empty ..\n";
                    write_stash();
                }


                it_j = stash.begin();
                while (it_j != stash.end()) {
                    if (VERBOSE) {
                        std::cout << "[ " << it_j->second << " ," << it_j->first << " ]" << " is in  the stash\n";
                    }
                    if (it_j->first == lf_id) {
                        // this is a good candate for eviction
                        stash_buffer.push_back(std::make_pair(it_j->first, it_j->second));
                        stash.erase(it_j);
                        if (VERBOSE){
                            std::cout << "stash size.................. " << get_size_stash() <<"\n";
                            std::cout << it_j->second << " added to stash_buufer\n";
                        }
                    }
                    it_j++;

                }
                std::cout << " Exit from while ..\n";
            } else if (VERBOSE)
                std::cout << "stash is empty ..\n";
            int num_candidates;
            if (stash_buffer.empty()) {
                num_candidates = 0;
                if (VERBOSE)
                    std::cout << "stash_buffer is empty ..\n";
            } else
                num_candidates = stash_buffer.size();

            Node_in_Path = (1 << (height - 1)) - 1 + lf_id;
            Node_t tmp;
            std::pair<int,  long long int> pair_tmp;

            while (num_candidates > 0 && Node_in_Path >= 0) {
                if (VERBOSE) {
                    std::cout << "num_candidates= " << num_candidates << " Node_in_Path= " << Node_in_Path << "\n";
                }

                for (uint32_t j = 0; j < ORAM_Z; j++) {
                    tmp.dummy = 0; // data
                    pair_tmp = stash_buffer.back(); // the last elements stored
                    stash_buffer.pop_back(); // and deleted
                    tmp.addr = pair_tmp.second;
                    tmp.leaf_id = pair_tmp.first;
                    if (VERBOSE)
                        std::cout <<  tmp.addr <<" " << tmp.leaf_id << "\n";
                    buckets[Node_in_Path * ORAM_Z + j] = tmp;
                    num_candidates--;
                    if (VERBOSE) {
                        std::cout << "REAL bucket[  " << Node_in_Path * ORAM_Z + j << " ]" << "\n";
                        std::cout << buckets[Node_in_Path * ORAM_Z + j].addr << " " << buckets[Node_in_Path * ORAM_Z + j].dummy << " " <<buckets[Node_in_Path * ORAM_Z + j].leaf_id <<"\n";
                    }

                }
                if (Node_in_Path == 0) break;
                Node_in_Path = (Node_in_Path - 1) / 2; // go to its parents along with the path up to the root
            }

            //write_buckets();


            // we need to fill the rest of buckets with dummies!
            if (Node_in_Path != 0)
            {
                while (Node_in_Path >= 0) {
                    if (VERBOSE) {
                        std::cout << "Node_in_Path= " << Node_in_Path << "\n";
                    }
                    for (uint32_t j = 0; j < ORAM_Z; j++) {
                        tmp.dummy = 1; // dummy
                        tmp.addr = -2;
                        tmp.leaf_id = -2;
                        if (VERBOSE) {
                            std::cout << "bucket[  " << Node_in_Path * ORAM_Z + j << " ]" << "\n";
                        }
                        buckets[Node_in_Path * ORAM_Z + j] = tmp;

                    }
                    if (Node_in_Path == 0) break;
                    Node_in_Path = (Node_in_Path - 1) / 2; // go to its parents along with the path up to the root


                } //while

            }//if












        }
        std::cout << " stash size: after  " << ntimes<<"   "<< get_size_stash() <<"\n";
    }

    std::cout << " Final stash size: " << get_size_stash() <<"\n";

}*/
void ORAM::update_histogram(int level, long long int value) {
    histogram[level] = value;
}

void ORAM::dummy_data_stat() {
    data_cache = 0;
    data_memory = 0;
    for (uint32_t i =0; i< ((1 << height) - 1); i++) {
        for (uint32_t j = 0; j < ORAM_Z; j++) {
            if (buckets[i * ORAM_Z + j].dummy != 1)
            {
                if (i < ((1 << hight_cached) - 1))
                {
                    // in cache
                    data_cache ++;
                }
                else
                {
                    data_memory++;
                }
            }

        }
    }
    //std::cout << data_cache <<" " << data_memory <<"\n";
}

void ORAM::write_histogran() {
    long long int total = 0;
    for (int i = 0; i< height; i++) {
        std::cout << i << " " <<std::dec<< histogram[i] << "\n";
        total = total + histogram[i];
    }

    std ::cout <<"=====================%%%%===================\n";

    for (int i = 0; i< height; i++) {
        std::cout << i << " " << (double)(100*histogram[i])/total  << "\n";
    }
}

long long int ORAM::get_miss_in_stash() {
    return (miss_in_stash);
}


unsigned long long int ORAM::get_time() {
    return (CYCLE_VAL);
}

unsigned  long long int ORAM::get_data_cache() {
    return data_cache;
}

unsigned  long long int ORAM::get_data_memory() {
    return data_memory;
}