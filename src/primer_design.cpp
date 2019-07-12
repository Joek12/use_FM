//
// Created by Joseph Kang on 2019-07-12.
//

#include "primer_design.h"
#include <iostream>
#include <math.h>

primer_design::primer_design() {
    // default recommended values from https://ocw.mit.edu/courses/biology/7-02-experimental-biology-communication-spring-2005/recitations/rdm6_allrecnotes.pdf
    rec_length.max = 25;
    rec_length.min = 18;
    rec_temp.max = 65;
    rec_temp.min = 55;

}

primer_design::primer_design(int min_length, int max_length, double min_temp, double max_temp) {
    rec_length.max = max_length;
    rec_length.min = min_length;
    rec_temp.max = max_temp;
    rec_temp.min = min_temp;
}

bool primer_design::dimer(const std::string * primer) {

}

long primer_design::nucToInt(const std::string kmer) {
/*
 * 'A' = 1
 * 'C' = 2
 * 'G' = 3
 * 'T' = 4
 */
    int temp = 0;
    long return_me = 0;
    for (size_t i = 0 ; i < kmer.size(); i++){
        switch(kmer[i]){
            case 'A':
                temp = 1;
                break;
            case 'C':
                temp = 2;
                break;
            case 'G':
                temp = 3;
                break;
            case 'T':
                temp = 4;
                break;
        }
        return_me += long(pow(10, kmer.size() - i - 1)) * temp;
    }
    return return_me;
}
