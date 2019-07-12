//
// Created by Joseph Kang on 2019-07-02.
//

#include <iostream>
#include "../resources/BitVector.h"
#include "../resources/openbwt.h"
#include "../resources/BWT.c"
#include "../resources/FMIndex.h"
#include "../resources/WaveletTree.h"
#include "../resources/misc.h"
#include "../resources/serializing.h"
#include <sstream>
#include <string>
#include "../resources/FMIndex.cpp"
#include "../resources/WaveletTree.cpp"
#include "../resources/BitVector.cpp"
#include "fastaReader.h"
#include "pbar.h"
#include <deque>
#include <queue>
#include "mu_commands.h"

#ifndef USE_FM_FMCOMMANDS_H
#define USE_FM_FMCOMMANDS_H


void make_fm(std::string inFile, std::string outFile){

    fastaReader fr;
    //std::string inFile = "/Users/josephkang/PycharmProjects/minimal_uniquemer/data/fake_geno_1/fake_geno.fa";

    // std::string read = fr.read(inFile);

    std::string unams = fr.read(inFile);

    std::cout << "Finished reading chars from : " << inFile << "\n";

    auto fmi = new FMIndex(unams);

    std::cout << "Finished building FM Index" << "\n";

    fmi->serialize_to_file(outFile);

    std::cout << "Finished serializing FM Index to file: " << outFile << "\n";


    // cout << fmi;

    // cout << read << endl;
    // cout << unams << endl;

    //std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;

    //cout << fmi->find(matches, std::string("CGGCGCGATCGGGGGCGCGATCCCCCCCAGAGAGA")) << endl;

}

FMIndex * make_fm(std::string inFile){

    fastaReader fr;
    //std::string inFile = "/Users/josephkang/PycharmProjects/minimal_uniquemer/data/fake_geno_1/fake_geno.fa";

    // std::string read = fr.read(inFile);

    std::string reads = fr.read(inFile);

    auto fmi = new FMIndex(reads);

    return fmi;

}

std::vector<int> find_matches(std::vector<std::string> seqs, FMIndex * fmi){

    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;
    std::vector<int> hits;
    for (std::string seq : seqs){
        hits.push_back(fmi->find(matches, seq));
    }
    return hits;


}

std::vector<int> find_matches(std::string filename, FMIndex * fmi, std::string geno){
    // find matches from json file

    // read the file (should be one continuous string read
    std::ifstream myf (filename);
    std::string line;
    if (myf.is_open()){
        getline(myf, line);
    }

    // parse the string
    // "{\"num0\": num1, \"num2\": num3, ... }"
    line = line.substr(1,line.length() - 2);

    // "\"num0\": num1, \"num2\": num3, ..."
    std::string delimiter = ",";
    std::string mini_delimiter = ":";
    size_t pos = 0;
    size_t mini_pos = 0;
    std::string token;

    int num0, num1;

    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;

    std::vector<int>hits;
    std::string seq;

    while ((pos = line.find(delimiter)) != std::string::npos){

        std::vector<int> nums;
        token = line.substr(0, pos);

        mini_pos = token.find(mini_delimiter);
        num0 = std::stoi(token.substr(1,mini_pos));
        token.erase(0, mini_pos + mini_delimiter.length() + 1); // 2 added to account for ": "

        num1 = std::stoi(token);

        line.erase(0, pos + delimiter.length() + 1);
        seq = geno.substr(num0, num1-num0);
        hits.push_back(fmi->find(matches, seq));

    }

    return hits;

}

std::vector<int> find_matches(std::vector<std::vector<int>> seqs, FMIndex * fmi, std::string geno){
    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;
    std::vector<int> hits;

    for (std::vector<int> mini : seqs){
        std::string seq = geno.substr(mini[0], mini[1] - mini[0]);
        hits.push_back(fmi->find(matches, seq));
    }

    return hits;

}

FMIndex* readFMFile(const std::string filename){

    FMIndex* fmi;
    fmi = fmi->new_from_serialized_file(filename);

    return fmi;

}

void simple_fm(){

    // FM Index
    const std::string ins{"mississippi"};

    FMIndex *fmi = new FMIndex(ins);

    auto result = fmi->findn("miss");

    std::cout << result << '\n';
}

void make_bwt(){

    // BWT
    std::string s = "mississippi";
    char t[s.size()];

    int pidx = BWT(reinterpret_cast<const unsigned char *>(s.c_str()),
                   reinterpret_cast<unsigned char *>(t),
                   static_cast<int>(s.size()));

    std::cout << t << '\n';
}

std::vector<int> * get_hits(std::deque<std::vector<int>> * reads, std::string * geno, FMIndex * fmi){
    std::vector<int> * hits = new std::vector<int>;
    pbar pb(reads->size(), "getting hits");
    for (std::vector<int> vec : *reads){
        int start = vec.at(0);
        int end = vec.at(1);
        std::string subs = geno ->substr(start, (end > start ? end-start : end));
        int hit = 0;

        subs.empty() ? NULL : hit = fmi->findn(subs);
        assert(hit <= 1);
        hits->emplace_back(hit);

        pb.update();
        pb.show();
    }
    return hits;

}


void assert_unique(std::deque<std::vector<int>> * reads, std::string* geno, FMIndex * fmi){
    size_t length = reads->size();
    int not_unique = 0;
    pbar pb(length, "asserting uniqueness");
    for (size_t i = 0; i < length; i++){
        auto vec = reads->at(0);
        reads->pop_front();
        int start = vec.at(0);
        int end = vec.at(1);

        // if end > start, most likely means that actual end position was stored
        // else, most likely means the length of kmer stored
        std::string subs = geno -> substr(start, (end > start ? end-start : end));

        int hit = 0;
        subs.empty() ? NULL : hit = fmi -> findn(subs);

        if (hit == 1) reads->emplace_back(vec);
        else ++not_unique;

        pb.update();
        pb.show();

    }

    std::cout << "not uniques: " << not_unique <<"\n";
    std::cout << "total number of original sequences: " << length << std::endl;

}

void write_finds(std::vector<std::vector<int>> * reads, std::string * geno, FMIndex * fmi, std::string * fn){

    std::ofstream myf (*fn);
    if (myf.is_open()){
        for (std::vector<int> vec : *reads){
            int start = vec.at(0);
            myf<<fmi->findn(geno->substr(start, vec.at(1) - start)) << ", ";
        }
    }
}

void write_start_end(std::deque<std::vector<int>> * reads, std::string start_file, std::string end_file){
    std::ofstream starts (start_file);
    std::ofstream ends (end_file);
    if (starts.is_open() and ends.is_open()){
        size_t length = reads->size();
        for(size_t i = 0; i < length; i++){
            auto vec = reads->at(i);
            starts << vec.at(0) << "\n";
            ends << vec.at(1) << "\n";
        }
    }
}

std::deque<std::vector<int>> * read_start_end(const std::string start_file, const std::string end_file){
    std::ifstream starts_if (start_file);
    std::ifstream ends_if (end_file);
    std::string start_line, end_line;

    if (starts_if.is_open() && ends_if.is_open()){
        // get lengths
        getline(starts_if, start_line);
        const std::string delim = "length";
        start_line.erase(0, delim.length() + 2);

        getline(ends_if, end_line);
        end_line.erase(0, delim.length() + 2);
        assert(std::stoi(start_line) == std::stoi(end_line));

        pbar pb (std::stoi(start_line), "reading starts and end");

        auto *reads = new std::deque<std::vector<int>>;
        std::vector<int> mini;

        while( getline(starts_if, start_line) && getline(ends_if, end_line)){
            mini.emplace_back(std::stoi(start_line));
            mini.emplace_back(std::stoi(end_line));
            reads->emplace_back(mini);
            mini.clear();

            pb.update();
            pb.show();
        }

        starts_if.close();
        ends_if.close();
        return reads;
    }
    return nullptr;

}

struct valid_SNP{
    int start_pos;
    int end_pos;
    int pos;
    char s;
};

std::queue<valid_SNP> check_snp_unique(std::deque<std::vector<int>> * reads, std::string * geno,  std::deque<SNP>* snps, FMIndex * fmi){
    // assumes that reads contains start as pos and ends as length

    // phase 1: popping unnecessary k-mers
    // phase 2: reading necessary k-mers

    size_t snp_arr_length = snps->size();
    pbar pb (snp_arr_length, "checking uniqueness after snp");
    size_t count = 0;

    int past_end = 0;
    int past_start = 0;

    std::queue<valid_SNP> v_snp;

    for (size_t i = 0; i < snp_arr_length; i++){

        auto snp = snps->at(0);
        snps->pop_front();

        // pop until acquiring first important kmer
        while(past_end < snp.pos){
            auto r_vec = reads->at(0);
            reads->pop_front();
            int start = r_vec.at(0);
            int end = r_vec.at(1);
            past_end = (end > start? end : start + end);
        }

        // reading phase
        while(past_start <= snp.pos){
            auto r_vec = reads->at(0);
            reads->pop_front();
            int start = r_vec.at(0);
            int end = r_vec.at(1);
            past_start = start;
            end = end > start? end - start : end;
            if (start <= snp.pos <= end + start){

                auto subs = geno->substr(start, end);
                if(fmi->findn(subs) == 1){
                    //++count;
                    valid_SNP sn = {start, end, snp.pos, snp.s};
                    v_snp.push(sn);
                }

            }
            else if (start > snp.pos) reads->push_front(r_vec);
        }


        pb.update();
        pb.show();

        /*
        if (holder.empty()) {
            int past_end = 0;
            while(past_end < snp.pos){
                auto r_vec = reads->at(0);
                int start = r_vec.at(0);
                int end = r_vec.at(1);
                past_end = (end > start? end : start + end);
                reads->pop_front();
            }
        }
         */




    }
    return v_snp;
}



void write_q_file(std::queue<valid_SNP> q, std::string fn){
    std::ofstream myf (fn);

    while(not q.empty()){

    }

}

#endif
