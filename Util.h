//
// Created by meysam on 9/13/18.
//

#ifndef SOFTWARE_ORAM_UTIL_H
#define SOFTWARE_ORAM_UTIL_H


#include <iostream>
#include <assert.h>
#include <time.h>
#include <list>
#include <random>

class Util{
private:
public:


    static size_t rand_int(size_t n) {
        // std::mt19937 ORAM::gen;
        return std::uniform_int_distribution<size_t>(0, n - 1)(gen);
        //  return rand() % n;
    }


    static std::random_device rd;
    static std::mt19937 gen;

};








#endif //SOFTWARE_ORAM_UTIL_H
