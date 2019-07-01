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

using namespace std;


int main() {

    // BWT
    string s = "mississippi";
    char t[s.size()];

    int pidx = BWT(reinterpret_cast<const unsigned char *>(s.c_str()),
                   reinterpret_cast<unsigned char *>(t),
                    static_cast<int>(s.size()));

    cout << t << endl;


    // FM Index
    const string ins{"mississippi"};

    FMIndex *fmi = new FMIndex(ins);

    auto result = fmi->findn("miss");

    // cout << result << endl;

    fastaReader fr;
    //std::string fn = "/Users/josephkang/PycharmProjects/minimal_uniquemer/data/fake_geno_1/fake_geno.fa";
    //std::string fn = "/Users/josephkang/Documents/uniquekmer/data/22.fa";
    std::string fn = "/work/bioinformatics/s187520/uniquekmer/data/genome/fa";

    // std::string read = fr.read(fn);

    std::string unams = fr.unam(fn);

    cout << "Finished reading unams from : " << fn << endl;

    fmi = new FMIndex(unams);

    cout << "Finished building FM Index" << endl;

    //std::string output_file = "/Users/josephkang/CLionProjects/use_FM/output/FM_index";
    std::string output_file = "/work/bioinformatics/s187520/use_FM/output/genome_FM_index";
    fmi->serialize_to_file(output_file);

    cout << "Finished serializing FM Index to file: " << output_file << endl;


    // cout << fmi;

    // cout << read << endl;
    // cout << unams << endl;

    // std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;

    // cout << fmi->find(matches, std::string("CGGCGCGATCGGGGGCGCGATCCCCCCCAGAGAGA")) << endl;



    return 0;

}
