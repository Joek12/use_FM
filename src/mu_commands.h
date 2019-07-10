//
// Created by Joseph Kang on 2019-07-10.
//

#ifndef USE_FM_MU_COMMANDS_H
#define USE_FM_MU_COMMANDS_H

#include <iostream>
#include <deque>
#include "../resources/FMIndex.h"
#include "pbar.h"


std::deque<std::vector<int>>* stitch_n(std::deque<std::vector<int>> * mu){

    size_t length = mu->size();


    int bundle_start = 0;
    int bundle_end = 0;

    for (size_t i = 0; i < length; i++){
        auto vec = mu->at(0);
        mu->pop_front();
        int start = vec.at(0);
        int end = vec.at(1);

        if (start > bundle_end){
            std::vector<int> rvec = {bundle_start, bundle_end};
            mu->push_back(rvec);
            bundle_start = start;
            bundle_end = start + end;
        }

        else{
            bundle_end = start + end;
        }
    }

    return mu;
}


void stitch(std::deque<std::vector<int>> * mu){

    size_t length = mu->size();

    pbar pb(length, "stitching");
    int bundle_start = 0;
    int bundle_end = 0;

    for (size_t i = 0; i < length; i++){
        auto vec = mu->at(0);
        mu->pop_front();
        int start = vec.at(0);
        int end = vec.at(1);

        if (start > bundle_end){
            std::vector<int> rvec = {bundle_start, bundle_end};
            mu->push_back(rvec);
            bundle_start = start;
            bundle_end = start + end;
        }

        else{
            bundle_end = start + end;
        }
        pb.update();
        pb.show();
    }

}

struct SNP{
    int pos;
    char s;
};

std::deque<SNP> read_snp_file(std::string pos, std::string ch){
    std::deque<SNP> reads;
    std::ifstream f1 (pos);
    std::ifstream f2 (ch);

    std::string line1;
    std::string line2;

    if (f1.is_open() and f2.is_open()){
        while(getline(f1, line1) and getline(f2, line2)) {
            SNP mys;
            mys.pos = std::stoi(line1);
            assert(line2.size() == 1);
            mys.s = line2[0];
            reads.push_back(mys);
        }

    }

    return reads;
}

#endif

