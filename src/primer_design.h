//
// Created by Joseph Kang on 2019-07-12.
//

#ifndef USE_FM_PRIMER_DESIGN_H
#define USE_FM_PRIMER_DESIGN_H
#include <iostream>


// General class to provide functions that check the following primer design guidelines:

/*
 * From "PCR Primer Design Basic Guidelines"
 *
 * Optimal Primer Characteristics:
 * 1. usually between 18-25 bp long
 * 2. Primer GC content ~50%
 * 3. Amplicon should not have G/C content in excess of 60%
 *     - (Ideally 40 - 60%)
 *     - High GC content may not denature well during cycling and also susceptible to non-specific interactions
 * 4. Balanced distribution of G/C and A/T bases
 * 5. 3' terminal position needs to be a G or C, but no more than 2 G or C in last 5 bases is best
 * 6. Avoid 3' terminal T
 * 7. No runs of single bases greater than three bases especially for G or C
 * 8. Tm of primer 55-65 C and pairs should not differ by more than 1-2 C
 * 9. Match primers for delta G more than Tm. Delta G no more than -10 Kcal/mol
 */

/*
 * From MIT Department of Biology's
 * 7.02 Experimental Biology and Communication, Spring 2005
 * "Primer Design"
 *
 * 1. Primers should **flank** the DNA that you want to amplify (i.e. one on either side), such that the exponentially
 *      amplified product consists of the primer sequences and everything in between them
 * 2. Primers are generally between 18-25 bp long
 * 3. Each primer should have a Tm between 55-65 C and G/C content of 50-60%
 * 4. Each primer should have a 3' end that hybridizes very well to the template, but the 5' end can be initially
 *      less complementary (or non-complementary) to the template
 * 5. The primers should not form "primer dimers" or "hairpins"
 */


class primer_design {
    struct range{
        double min;
        double max;
    };

    range rec_length;
    range rec_temp;

    unsigned long _inner_nuc2int(std::string kmer, int nucs);



public:
    bool length_check(std::string kmer);
    double gc_content(std::string kmer);
    double distribution(std::string kmer);
    bool good_3_end(std::string kmer);
    bool good_last_5(std::string kmer);


    primer_design();
    primer_design(int min_length, int max_length, double min_temp, double max_temp);
    bool dimer(std::string primer); // checks hairpins as well

    unsigned long nucToInt(std::string kmer, bool gc_at_switch);

    bool optimal(std::string kmer, bool verbose);



};


#endif //USE_FM_PRIMER_DESIGN_H
