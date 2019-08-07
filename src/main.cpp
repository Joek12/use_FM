#include <iostream>
#include <fstream>
#include "FMCommands.h"
#include "fastaReader.h"
#include "pbar.h"
#include "mu_commands.h"
#include "primer_design.h"
//#include "../test/test_primer_design.h"
#include <queue>

auto json_reader(std::string filename){

    // read the file (should be one continuous string read
    std::ifstream myf (filename);
    std::string line;
    if (myf.is_open()){
        getline(myf, line);
    }

    // parse the string
    // "{\"num0\": num1, \"num2\": num3, ... }"
    std::vector<std::vector<int>> *reads;
    line = line.substr(1,line.length() - 2);

    // "\"num0\": num1, \"num2\": num3, ..."
    std::string delimiter = ",";
    std::string mini_delimiter = "\"";
    size_t pos = 0;
    size_t mini_pos = 0;
    std::string token;


    while ((pos = line.find(delimiter)) != std::string::npos){

        std::vector<int> nums;
        token = line.substr(0, pos);

        token.erase(0, token.find(mini_delimiter) + mini_delimiter.length());

        mini_pos = token.find(mini_delimiter);
        nums.emplace_back(std::stoi(token.substr(0,mini_pos)));
        token.erase(0, mini_pos + mini_delimiter.length() + 2); // 2 added to account for ": "

        nums.emplace_back(std::stoi(token));

        reads->emplace_back(nums);
        line.erase(0, pos + delimiter.length());
    }

    return reads;

}

auto just_read_it(std::string filename){

    std::ifstream myf (filename);
    std::string line;
    std::vector<int> reads;
    if (myf.is_open()) {
        while (getline(myf, line)){
            reads.emplace_back(std::stoi(line));
        }
    }
    return reads;

}

auto just_dump_read(std::string filename){

    // read the file (should be one continuous string read)
    std::ifstream myf (filename);
    std::string line;
    if (myf.is_open())
        getline(myf, line);

    //parse string
    std::vector<std::vector<int>> reads;

    std::string delimiter = ",";
    std::string mini_delim = " ";
    size_t pos, mini_p = 0;
    //size_t pos = 0;
    std::string token;

    while (( pos = line.find(delimiter)) != std::string::npos){
        std::vector<int> nums;
        token = line[0] == ' ' ? line.substr(1, pos) : line.substr(0, pos);

        mini_p = token.find(mini_delim);
        nums.emplace_back(std::stoi(token.substr(0, mini_p)));
        token.erase(0, token.find(mini_delim) + mini_delim.length());

        nums.emplace_back(std::stoi(token));

        //token.erase(0, token.find(delimiter) + delimiter.length());

        reads.emplace_back(nums);
        line.erase(0, pos + delimiter.length());

    }

    return reads;
}

void write_hits_file(std::string fn, std::vector<int> * hits){
    std::ofstream myf (fn);
    if (myf.is_open()){
        for (int hit: *hits)
            myf<<hit<<", ";
    }
}

void test_fmi(){
    FMIndex fmi = FMIndex("mississippi");
    std::cout << "siss hits: " << fmi.findn("siss") << "\n";
    std::cout << "issis hits: " << fmi.findn("issis") << "\n";
    std::cout << "ssiss hits: " << fmi.findn("ssiss") << "\n";
    std::cout << "sis hits: " << fmi.findn("sis") << "\n";
    std::cout << "si hits: " << fmi.findn("si") << "\n";

}

#ifndef CHECK_AS_PRIMERS
#define CHECK_AS_PRIMERS

void check_as_primers(std::deque<std::vector<int>> * reads, const std::string * geno){
    try {
        primer_design pd;
        int start = 0;
        int end = 0;
        int count = 0;

        size_t length = reads->size();

        pbar pb(length, "checking valid primers: ");
        for (size_t i = 0; i < length; i++) {

            auto mini = reads->at(0);

            start = mini.at(0);

            end = (mini.at(1) < start ? mini.at(1) : mini.at(1) - start);

            reads->pop_front();

            auto kmer = geno->substr(start, end);

            if (pd.optimal(kmer, false)) {
                ++count;
                reads->emplace_back(mini);
            }

            pb.update();
            pb.show();
        }

        std::cout << "number of valid kmers as primers: " << count << '\n';
    }
    catch(std::out_of_range){
        std::cout<<"error caught in check_as_primers()";
    }

}

#endif


void work(){

    //iMac prefix
    const std::string prefix = "/Users/josephkang/Documents/uniquekmer";

    //22.fa
    const std::string chr22 = prefix + "/data/22.fa";

    //genome.fa
    const std::string genome = prefix + "/data/genome.fa";

    // use FM prefix
    const std::string fm_prefix = "/Users/josephkang/CLionProjects/use_FM";

    // fm indices
    const std::string chr22_fm = fm_prefix + "/output/22_FM";
    const std::string geno_fm = fm_prefix + "/output/genome_FM_index";

    // MU start and end
    const std::string mu_start = prefix + "/src/c22_mu_starts_0709";
    const std::string mu_end = prefix + "/src/c22_mu_ends_0709";

    // SNP
    const std::string snp_pos = prefix + "/output/ez_c22_SNP_pos";
    const std::string snp_char = prefix + "/output/ez_c22_SNP_char";

    // read the MU's from the mu start and ends
    auto reads = std::deque<std::vector<int>>();
    read_start_end(&reads, mu_start, mu_end);
    //stitch(reads);

    std::cout << "number of mu read: " << reads.size() << '\n';

    // stitch the reads
    //stitch(reads);
    //perfect_stitch(reads);
    //distance(reads);


    fastaReader fr;

    // read the genome file
    std::string geno = fr.read(chr22);
    //std::string geno = fr.read(genome);

    // read the SNPs
    auto snps = std::vector<SNP>();

    read_snp_file( &snps, snp_pos, snp_char);
    std::cout << "number of snps read: " << snps.size() << '\n';

    FMIndex *c22 = readFMFile(chr22_fm);
    //FMIndex *c22 = readFMFile(geno_fm);

    // change the mu's to only the unique ones
    assert_unique(&reads, &geno, c22);

    // write the start and end of mu after assert
    write_start_end(&reads, "mu_unique_starts", "mu_unique_ends");


    //check_as_primers(&reads, &geno);

    //std::cout << "number of valid primers received: " << reads.size() << '\n';


    auto v_snp = std::queue<valid_SNP>();
    // check if kmers are still unique after SNP
    auto reads_copy = reads;
    check_snp_unique(&v_snp, &reads_copy, &geno, &snps, c22);


    check_as_primers(&reads, &geno);

    std::cout << "number of valid snps: " << v_snp.size() << '\n';
    std::cout << "size of reads afterwards: " << reads.size() << '\n';

}





int main() {

    work();


    //test_distribution();
    //test_nucToInt();
    //test_gc_content();
    //test_optimal();

    return 0;
}
