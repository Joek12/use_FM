//
// Created by Joseph Kang on 2019-06-27.
//

#ifndef USE_FM_FASTAREADER_H
#define USE_FM_FASTAREADER_H
#include <string>

class fastaReader {

    std::string filename;
    std::string sequence;

public:

    fastaReader(std::string fn);
    fastaReader();
    std::string read(std::string fn);
    std::string unam(std::string fn);
    std::string unam();

};


#endif //USE_FM_FASTAREADER_H
