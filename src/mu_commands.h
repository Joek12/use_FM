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

void distance(std::deque<std::vector<int>> * mu, int d = 300){

    auto vec = mu->at(0);
    size_t length = mu->size();
    int past_start = vec.at(0);
    int past_end = (vec.at(1) > past_start ? vec.at(1) : (past_start + vec.at(1)));
    bool flag = true;

    pbar pb (length, "finding within distance");

    for (size_t i = 0; i < length; i++){
        vec = mu->at(0);
        mu->pop_front();
        int start = vec.at(0);

        if (start <= past_end + d){
            if (!flag){
                std::vector<int> v = {past_start, past_end};
                mu->push_back(v);
            }
            flag = true;
            mu->push_back(vec);
        }

        else flag = false;

        past_end = vec.at(1);
        past_start = start;

        pb.update();
        pb.show();
    }

    std::cout << "number of valids after distance()" << mu->size();
}

void perfect_stitch(std::deque<std::vector<int>> * mu, int bot = 0, int top = 0){


    auto vec = mu->at(0);
    size_t length = mu->size();
    int past_start = vec.at(0);
    int past_end = (vec.at(1) > past_start ? vec.at(1) : (past_start + vec.at(1)));
    int begin = past_start;
    pbar pb(length, "perfect stitching");
    int count_mu = 0;

    try {

        for (size_t _ = 0; _ < length; _++) {

            pb.update();
            pb.show();

            //vec = mu->at(0);
            if (mu->empty())
                break;
            mu->pop_front();
            int start = vec.at(0);
            count_mu++;
            int end = (vec.at(1) > start ? vec.at(1) : vec.at(1) + start);


            if (start != past_start + 1 or end != past_end + 1 or count_mu > past_end - past_start + 1) {
                if (past_start != begin and ((bot == 0 and top == 0) or bot < past_end - begin + 1 < top) and
                    count_mu == past_end - past_start + 1) {
                    std::vector<int> mini_v = {begin, past_end};
                    mu->push_back(mini_v);
                }

                begin = start;
                past_end = end;
                count_mu = 0;
            }

            past_start = start;

        }

        std::cout<< "number of perfect stitches: " << mu->size() << std::endl;
    }

    catch(const std::out_of_range& oor){
        std::cout << "size of mu until exception: " << mu->size() << std::endl;
        return;
    }


}



#endif


