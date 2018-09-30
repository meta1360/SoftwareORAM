#include <iostream>

#include "ORAM.h"
#include "time.h"
#include <fstream>
#include <istream>
#include <ostream>
#include <iostream>
#include <assert.h>
#include <string>
#include <sstream>      // std::istringstream
#include <cstring>


const bool M_VERBOSE = 1;
const unsigned long long int NUM_ACC = 60000000;

std::ifstream myfile;
//const unsigned long long int NUM_NODE = 64*1024*1024;
//const int DATA_PERCENT = 80;
//const unsigned long long int INIT_CACH_LINE = ((NUM_NODE * ORAM_Z * DATA_PERCENT)/100);


//g++ -std=c++11 -o main Util.cpp ORAM.cpp main.cpp
// to initiate the ORAM with some randome number
// size_ORAM shows how big our addresses are
// to initialize the ORAM
//size_ is the number of cache lines gonna be initilized.
void initialize_ORAM(ORAM & my_ORAM, const unsigned long long int & size_)
{
    std::cout << "initialization process starting .....\n";
    //unsigned long long int num_block = size_/64;
    for (long long int i=0; i<size_; i=i+64)
    {
        //std::cout << "write addr " << i <<"\n";
        if (i%(1000000) == 0) {
            std::cout << std::dec << time(NULL) << std::dec << "  " << i << " has been accessed to ORAM \n";
            std::cerr << std::dec << time(NULL) << std::dec << "  " << i << " has been accessed to ORAM \n";
        }
        // std::cout << "----------------starting time = " <<std::dec<< time(NULL) << "\n";
        my_ORAM.access('w',i);
        // std::cout << "----------------end time = " <<std::dec<< time(NULL) << "\n";

        // std::cout <<std::dec<<time(NULL)<< " STASH --> " << my_ORAM.get_size_stash() << "\n";


    }
    std::cout << "size oram = " << size_ << " done \n";

    my_ORAM.write_buckets();

}

void read_ORAM(ORAM &my_ORAM, const unsigned long long int & size_)
{
    ///std::cout << "reading process starting .....\n";

    for (long long int i=0; i<size_; i=i+64)
    {
        std::cout << "read addr " << i <<"\n";
        if (i%10000 == 0)
            std::cout << i <<" has been read from ORAM \n";
        my_ORAM.access('r',i);


    }
}


/*void read_write_ORAM(ORAM & my_ORAM)
{

}*/

// ./main num_node zi data_percent
int main(int argc, char *argv[]) {
    long long int n = (long long int) atoll(argv[1]);
    int ORAM_HIGHT = (uint32_t)ceil(log2((double) n));
    int num_node = (1 << ORAM_HIGHT)- 1;
    uint32_t  zi =  (uint32_t)atoll(argv[2]);
    int utilization = (int)atoll(argv[3]);

    double avg_data_memory = 0;
    double avg_data_cache = 0;
    unsigned long long int total_data_memory = 0;
    unsigned long long int total_data_cache = 0;


    unsigned  long long int prv_end_time = 0;
    long long int current_stash_size = 0;


    std::string prefix = "/uusoc/scratch/res/arch/students/meysam/input/";
    std::string file_n = prefix + argv[4];

    long long int INIT_CACH_LINE = (((zi * num_node * utilization)/100));

    std::cout << num_node <<" " <<INIT_CACH_LINE <<"\n";
    int CACHE_HIGHT = (ORAM_HIGHT > 4) ? (ORAM_HIGHT - 4) : 2;//(log2((((1<< ORAM_HIGHT) * 1)/16)) == 0) ? 2 : log2((((1<< ORAM_HIGHT) * 1)/16)) ;
    ORAM my_ORAM(num_node, zi, CACHE_HIGHT, utilization); //64*1024 = num_node total size X*2*64 = 256X Byte
    std::cout << "starting time = " <<std::dec<< time(NULL) << "\n";

    initialize_ORAM(my_ORAM,INIT_CACH_LINE*64);  // 128*1024 is the total number of cache line
    std::cout << "ending time = " <<std::dec<< time(NULL) << "\n";

    myfile.open(file_n, std::ios::in);
    assert(myfile.is_open());


    /// my_ORAM.SnapShot();
/*
    if (M_VERBOSE) {
        std::cout << "==========================> bucket <=========================\n";
        my_ORAM.write_buckets();
        std::cout << "=============>   stash size = " << my_ORAM.get_size_stash()<<  "   <=====================" << "\n";
        my_ORAM.write_stash();
        std::cout << "=============================================================\n";
        std::cout << "after stash size=  " << my_ORAM.get_size_stash() << "\n";
    }

    */
    //read_ORAM(my_ORAM,30*64);
    std::string line;
    long long int addr = INIT_CACH_LINE - 1;
    long long int instruction_count = 0;
    //for (unsigned long long int ii = 0; ii< NUM_ACC; ii++)
    long long int stash_size_prev = 1;


    while ((!myfile.eof()) && (instruction_count < NUM_ACC))
    {
        if (instruction_count == 10000000)
            my_ORAM.reset_hit_miss();
        getline(myfile, line);
        //std::cout << line << "\n";
        std::istringstream iss (line);
        //(line[2],gg);
        std::string word;
        int count = 0;
        bool W_R = 0;
        while (iss >> word )
        {
            if (count == 1)
            {
                if (word == "W")
                    W_R = true;
                else if (word == "R")
                    W_R = false;
                else{
                    std::cout << "trace is screwed up!\n";
                    assert(0);
                }
            }
            else
            if (count == 2)
            {
                //long long int raw_address = (rand()%INIT_CACH_LINE)*64; //std::stoll(word,nullptr,16);
                long long int address = (rand()%INIT_CACH_LINE)*64; //std::stoll(word,nullptr,16);
                // std::cout << "$$$$$$$$$$$$$$$$" << raw_address <<"\n";
                //getchar();
                //int num_shift = ((double)address_bit);
                //std::cout << "gggg " << num_shift << " " <<std::hex<< (((long long int)1<<(num_shift)) - 1) <<"\n";
                //long long int address = (((((((long long int)1<<(num_shift)) - 1) & raw_address) >> 6)<<6)) & (((long long int)1 << ((long long int )log2(INIT_CACH_LINE*64))) -1 );
                //std::cout << "raw address = " <<std::hex<< raw_address << " " <<address  <<std::dec<< "  address bit = " << address_bit  <<" " <<
                //                                                                           ((long long int) (log2(INIT_CACH_LINE*64)))<<" "<< (INIT_CACH_LINE*64) << "\n";

                // std::cout << "raw address = " <<std::hex<< raw_address << " " <<address <<"\n";

                if (W_R)
                {
                    // write
                    my_ORAM.access('w', address);
                    //std::cout <<" considered as a write " <<std::hex<< address <<"\n";
                }
                else
                {
                    // read
                    my_ORAM.access('r', address);
                    //std::cout <<" considered as a read " <<std::hex<< address <<"\n";


                }






            }
            count++;
        }


        if (instruction_count%1000 == 0)
        {
            my_ORAM.dummy_data_stat();

            total_data_cache = total_data_cache  + my_ORAM.get_data_cache();

            // std::cout << "ggggg  " <<total_data_cache <<" "<< my_ORAM.get_data_cache()<<"\n";

            // getchar();

            total_data_memory = total_data_memory  + my_ORAM.get_data_memory();

            //std::cout << "hhhh  " <<total_data_ <<" "<< my_ORAM.get_data_cache()<<"\n";

            if (instruction_count%1000001 == 0) {

                std::cout << instruction_count <<"\n";

                //assert(0);

                std::cout << "--> " << my_ORAM.get_data_cache() << " " << my_ORAM.get_data_memory() << " "
                          << my_ORAM.get_data_cache() + my_ORAM.get_data_memory() << " " << "stash size = "
                          << my_ORAM.get_size_stash() << " " <<INIT_CACH_LINE << "\n";
                assert((my_ORAM.get_data_cache() + my_ORAM.get_data_memory() + my_ORAM.get_size_stash() == INIT_CACH_LINE ) ||
                       (my_ORAM.get_data_cache() + my_ORAM.get_data_memory() + my_ORAM.get_size_stash() == INIT_CACH_LINE ));
            }

        }


        if (instruction_count%EPOC == 0) {
            current_stash_size = my_ORAM.get_size_stash();
            //  rate_stash_length = current_stash_size - prev_stash_size;
            //  prev_stash_size = current_stash_size;

            //if (rate_stash_length > STASH_DRAIN_RATE) {

            if (!background_evict_off) {
                unsigned long long int start_time = my_ORAM.get_time();
                std::cout <<std::dec<< "start time = " << my_ORAM.get_time() <<" ";
                unsigned long long int  total_time = my_ORAM.get_time() - prv_end_time;
                while (my_ORAM.get_size_stash() > 1.2*STASH_DRAIN_LENGTH) {
                    // std::cout << "ALARM: rate exceeds the threshold ....  " << std::dec
                    //          << my_ORAM.get_total_size_stash() << "\n";
                    my_ORAM.background_eviction();
                    // std::cout << " Background eviction happened ...  " << my_ORAM.get_total_size_stash() << "\n";
                }
                std::cout <<std::dec << "end time = " << my_ORAM.get_time() <<" " << (double)(my_ORAM.get_time() - start_time)/total_time  <<"\n";
                prv_end_time = my_ORAM.get_time();

            }
        }

        //std::cout << total_data_cache <<" "<< my_ORAM.get_data_cache()<<" "<< instruction_count <<"\n";
        //std::cout << total_data_memory <<" "<< my_ORAM.get_data_memory()<<" "<< instruction_count <<"\n";

        avg_data_memory = (double)total_data_memory/(instruction_count/1000);
        //if (instruction_count%1000000 == 0)

        avg_data_cache = (double)total_data_cache/(instruction_count/1000);



        // std::cout << total_data_cache <<" "<< my_ORAM.get_data_cache()<<" "<< instruction_count << " " <<avg_data_cache<< "\n";

        // if (instruction_count > 1000)
        //     getchar();

        if ((instruction_count%1000000 == 0)  && (instruction_count > 14000000 ))
        {
            std::cout <<std::dec<<time(NULL)<< std::dec <<" " << instruction_count << " read accesses done" <<  " stash size is " << std::dec << my_ORAM.get_size_stash()
                      << "  " << std::dec << " " <<  " found in stash  = " << my_ORAM.get_found_in_stash() << " " << "miss stash = " <<
                      my_ORAM.get_miss_in_stash() << " avg data memory " << avg_data_memory << " avg data cache " << avg_data_cache << "\n";
            my_ORAM.write_ONC_stat();
            std::cerr <<std::dec<<time(NULL)<< std::dec <<" " << instruction_count << " read accesses done" <<  " stash size is " << std::dec << my_ORAM.get_size_stash() <<
                      " found in stash = " << my_ORAM.get_found_in_stash()<<"\n";

            //my_ORAM.write_stash();
            //  std::cout << "_____________________________________________\n";
            //  my_ORAM.write_buckets();
        }
        if (instruction_count > 14000000 ) {
            std::cout << std::dec << instruction_count << "  stash --> " << my_ORAM.get_size_stash() << "\n";
            if (my_ORAM.get_size_stash() > 10)
                my_ORAM.write_ONC_stat();
        }
        /*if (stash_size_prev < my_ORAM.get_size_stash() )
        {
            stash_size_prev = my_ORAM.get_size_stash();
             std::cout <<std::dec<<time(NULL)<< std::dec <<" " << instruction_count << " stash size alarm  " <<  " stash size is " << std::dec << my_ORAM.get_size_stash() <<"\n";
            std::cerr <<std::dec<<time(NULL)<< std::dec <<" " << instruction_count << "  stash size alarm " <<  " stash size is " << std::dec << my_ORAM.get_size_stash() <<"\n";
            my_ORAM.write_ONC_stat();

           // my_ORAM.set_VERBOSE();
            //my_ORAM.write_stash();
           // std::cout << "_____________________________________________\n";

            //my_ORAM.write_buckets();


        }*/
        if ((instruction_count > 14931581) && (instruction_count < 14932291))
        {
            if (instruction_count == 76001)
                std::cout <<"Start monitoring \n";
            my_ORAM.set_VERBOSE();
            if (instruction_count%100 == 0)
            {
                my_ORAM.write_stash();
                std::cout <<"=====================BUCKET  ===================\n";
                my_ORAM.write_buckets();
                std::cout <<"=====================ONC ====================\n";
                my_ORAM.write_ONC();
            }
        }
        else {
            if (instruction_count > 14932291) {
                my_ORAM.reset_VERBOSE();
                std::cerr << "monitoring done \n";
                getchar();
            }
        }

        instruction_count++;
    }

    //my_ORAM.stash_eviction();
    return 0;
}