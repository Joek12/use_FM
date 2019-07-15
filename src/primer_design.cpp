//
// Created by Joseph Kang on 2019-07-12.
//

#include "primer_design.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <sstream>

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


unsigned long primer_design::nucToInt(const std::string kmer, bool gc_at_switch) {
/*
 * 'A' = 1
 * 'C' = 2
 * 'G' = 3
 * 'T' = 4
 */

    return (gc_at_switch? _inner_nuc2int(kmer, 2) : _inner_nuc2int(kmer, 4));

}

bool primer_design::length_check(const std::string kmer) {
    return (rec_length.min <= kmer.size() and kmer.size() <= rec_length.max);
}


double primer_design::gc_content(const std::string kmer) {
    // distribution() implicitly returns info on gc content!

    size_t g_count = std::count(kmer.begin(), kmer.end(), 'G');
    size_t c_count = std::count(kmer.begin(), kmer.end(), 'C');
    return double(g_count + c_count) / kmer.size();
}


double primer_design::distribution(const std::string kmer) {
    // convert the input k-mer into a string with '1' representing A/T chars and '2' representing G/C
    /*
    std::stringstream ss;
    std::string k_str;
    ss << _inner_nuc2int(kmer, 2);
    ss >> k_str;
     */

    /*check 2-mers. The value for number of "12" and "21" corresponds to how well distributed the k-mer is
     * a perfectly distributed k-mer (ex. ACACACACACAC) should return a value of 1.0
     *
     * the higher presence of "22" or "11" should lower the value returned
     *  (i.e. the substring is less distributed)
    */

    // count the number of "12" or "21" pairs
    int count = 0;
    char past = ' ';
    for (char c : kmer){
        if (past != ' '){
            // 21 condition
            if ((c == 'A' or c == 'T') and (past == 'G' or past == 'C')) ++count;
            // 12 condition
            else if ((c == 'G' or c == 'C') and (past == 'A' or past == 'T')) ++ count;
        }
        past = c;
    }

    // total number of possible 2-mers from the input is the number of chars - 1
    return double(count) / (kmer.size() - 1);

}


bool primer_design::good_3_end(std::string kmer) {
    return (kmer[kmer.size() - 1] != 'T');
}


unsigned long primer_design::_inner_nuc2int(std::string kmer, int nucs) {

    int temp = 0;
    long return_me = 0;

    if (nucs == 2){
        for (size_t i = 0; i < kmer.size(); i++){
            switch(kmer[i]){
                case 'A':
                case 'T':
                    temp = 1;
                    break;
                case 'C':
                case 'G':
                    temp = 2;
                    break;
            }
            return_me += long(pow(10, kmer.size() - i - 1)) * temp;
        }
        return return_me;
    }
    else if (nucs == 4) {
        for (size_t i = 0; i < kmer.size(); i++) {
            switch (kmer[i]) {
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

    else{
        std::cerr << "Unknown quantity entered for _inner_nuc2int()";
        return 0;
    }
}

bool primer_design::good_last_5(std::string kmer) {
    if (kmer.size() < 5) return false;
    std::string mini_s = kmer.substr(kmer.size() - 5, kmer.size());
    size_t count = std::count(mini_s.begin(), mini_s.end(), 'G') + std::count(mini_s.begin(), mini_s.end(), 'C');
    return (count <= 2 and (mini_s[mini_s.size()-1] == 'G' or mini_s[mini_s.size()-1] == 'C'));
}

bool primer_design::optimal(const std::string kmer, const bool verbose) {
    // run through all the tests and return whether this kmer is optimal as a primer
    bool flag = true;

    if(not length_check(kmer)){

        if (verbose) {
            std::cout << kmer << " is not optimal as a primer because of the length." << '\n';
            flag = false;
        }

        else return false;

    }

    auto gc_c = gc_content(kmer);
    if (gc_c < 0.4 or gc_c > 0.6){

        if (verbose){
            std::cout << kmer << " is not optimal as a primer because of its G/C content (" << gc_c << ")."<< '\n';
            flag = false;
        }

        else return false;
    }

    // will choose 0.5 as the minimum distribution score
    if (distribution(kmer) < 0.5){

        if (verbose){
            std::cout << kmer << " is not optimal as a primer because of its distribution." << '\n';
            flag = false;
        }

        else return false;
    }

    if (not good_3_end(kmer)) {

        if (verbose){
            std::cout << kmer << " is not optimal as a primer because of the presence of 'T' in its 3' end." << '\n';
            flag = false;
        }

        else return false;
    }

    if (not good_last_5(kmer)){

        if (verbose){
            if (kmer.size() >= 5)
                std::cout << kmer << " is not optimal as a primer because of the last 5 bp" << " (" << kmer.substr(kmer.size() - 5, kmer.size())
                << ") " << '\n';

            flag = false;
        }

        else return false;
    }

    return flag;
}


