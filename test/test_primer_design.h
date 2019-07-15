//
// Created by Joseph Kang on 2019-07-15.
//

#include "../src/primer_design.h"
using namespace std;

void test_distribution(){

    primer_design pd;
    vector<string> vs;
    vector<double> expected;

    string test1 = "ACACACACACA"; // perfect score: 1
    vs.emplace_back(test1);
    expected.emplace_back(1);

    string test2 = "AAAAAAAAAAA"; // worst score: 0
    vs.emplace_back(test2);
    expected.emplace_back(0);

    string test3 = "ACGTACGTACGT"; // 6 out of 11 -> 0.545
    vs.emplace_back(test3);
    expected.emplace_back(double(6)/11);

    string test4 = "ATATATATATAT"; // worst score: 0
    vs.emplace_back(test4);
    expected.emplace_back(0);

    string test5 = "AAAAAGGGGGG"; // 1 out of 10: 0.1
    vs.emplace_back(test5);
    expected.emplace_back(double(1)/10);

    size_t len = vs.size();

    for (size_t i = 0; i < len; i++){

        std::string ms = vs.at(i);

        auto result = pd.distribution(ms);

        assert(expected.at(i) == result);

        cout << "check distribution of " << ms << ": " << result << "\n";


    }


}

void test_nucToInt(){
    primer_design pd;

    //test case with 4 nucs
    string kmer = "ACGTACGTACGT";
    auto result = pd.nucToInt(kmer, false);
    cout << "nucToInt of " << kmer << ": " << result << '\n';
    assert(result == 123412341234);

    result = pd.nucToInt(kmer, true);
    cout << "nucToInt of " << kmer << ": " << result << '\n';
    assert(result == 122112211221);

    kmer = "AAAAAAAAAAAA";
    result = pd.nucToInt(kmer, false);
    cout << "nucToInt of " << kmer << ": " << result << '\n';
    assert(result == 111111111111);


}

void test_gc_content(){
    primer_design pd;

    string kmer = "ACGTACGTACGT";
    auto result = pd.gc_content(kmer);
    cout << "gc_content() of " << kmer << ": " << result << '\n';


    kmer = "AAAAAAAAAAAA";
    result = pd.gc_content(kmer);
    cout << "gc_content() of " << kmer << ": " << result << '\n';

    kmer = "GCGCGCGCGCGCGCG";
    result = pd.gc_content(kmer);
    cout << "gc_content() of " << kmer << ": " << result << '\n';

    kmer = "GGGGGGGGGGGGGG";
    result = pd.gc_content(kmer);
    cout << "gc_content() of " << kmer << ": " << result << '\n';

    kmer = "CCCCCCCCCCCCCC";
    result = pd.gc_content(kmer);
    cout << "gc_content() of " << kmer << ": " << result << '\n';
}

void test_optimal(){
    primer_design pd;

    vector<string> kmers = {
            "AG",
            "ACGTACGTACGT",
            "AAAAAAAAAAAA",
            "GTACGATAAGCTAAGTA",
            "GTACGATAAGCTGGGGGG",
            "GAGCCTGAGTAGATAG"
    };

    vector <bool> results;

    cout << "VERBOSE MODE: \n_________________\n\n\n";

    for (string s : kmers){
        auto result = pd.optimal(s, true);
        cout << s << " is " << (result? "a " : "not a ") << "good primer.\n\n";
    }


    cout << "\n\nQUIET MODE: \n_________________\n\n\n";

    for (string s : kmers){
        cout << s << " is " << (pd.optimal(s, false)? "a " : "not a ") << "good primer." << '\n';
    }

}







