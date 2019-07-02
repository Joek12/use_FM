#include <iostream>

#include "python.hpp"

#include "FMCommands.h"

#define PY_SSIZE_T_CLEAN
#include "Python/Python.h"


using namespace std;






int main() {

    const std::string inFile = "/Users/josephkang/Documents/uniquekmer/data/22.fa";
    const std::string outFile = "/Users/josephkang/CLionProjects/use_FM/output/22_FM";

    // make_fm(inFile, outFile);
    FMIndex* c22 =readFMFile(outFile);

    std::vector<std::string> sequences;
    sequences.emplace_back("GAATTCTTGTGTTT");
    sequences.emplace_back("AATTCTTGTGTTT");

    auto hits = find_matches(sequences, c22);




    return 0;

}
