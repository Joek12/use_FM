//
// Created by Joseph Kang on 2019-07-12.
//

#ifndef USE_FM_PRIMER_DESIGN_H
#define USE_FM_PRIMER_DESIGN_H
#include <iostream>


class primer_design {
    struct range{
        double min;
        double max;
    };

    range rec_length;
    range rec_temp;

public:

    primer_design();
    primer_design(int min_length, int max_length, double min_temp, double max_temp);
    bool dimer(const std::string * primer);

    long nucToInt(std::string kmer);


};


#endif //USE_FM_PRIMER_DESIGN_H
