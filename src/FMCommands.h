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

#ifndef USE_FM_FMCOMMANDS_H
#define USE_FM_FMCOMMANDS_H

#endif //USE_FM_FMCOMMANDS_H





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

std::vector<int> * get_hits(std::vector<std::vector<int>> * reads, std::string * geno, FMIndex * fmi){
    std::vector<int> * hits;
    for (std::vector<int> vec : *reads){
        int start = vec.at(0);
        hits->emplace_back(fmi->findn(geno->substr(start, vec.at(1) - start)));
    }
    return hits;

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
