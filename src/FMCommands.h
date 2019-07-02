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

    std::string unams = fr.unam(inFile);

    std::cout << "Finished reading unams from : " << inFile << "\n";

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

std::vector<int> find_matches(std::vector<std::string> seqs, FMIndex * fmi){

    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;
    std::vector<int> hits;
    for (std::string seq : seqs){
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
