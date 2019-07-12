#include <iostream>
#include <fstream>
#include "FMCommands.h"
#include "fastaReader.h"
#include "pbar.h"
#include "mu_commands.h"
#include "primer_design.h"

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

int main() {

    const std::string inFile = "/Users/josephkang/Documents/uniquekmer/data/22.fa";
    const std::string outFile = "/Users/josephkang/CLionProjects/use_FM/output/22_FM";
    const std::string json_file = "/Users/josephkang/Documents/uniquekmer/output/c22_mu";

    const std::string mu_start = "/Users/josephkang/Documents/uniquekmer/src/c22_mu_starts_0709";
    const std::string mu_end = "/Users/josephkang/Documents/uniquekmer/src/c22_mu_ends_0709";

    const std::string snp_pos = "/Users/josephkang/Documents/uniquekmer/output/ez_c22_SNP_pos";
    const std::string snp_char = "/Users/josephkang/Documents/uniquekmer/output/ez_c22_SNP_char";

    primer_design pd;
    std::string dum = "ATTCA";
    std::cout << "A + G " << 'A' + 'G' << " " << std::endl;
    std::cout << "find G in ATTCA " << (dum.find('G') != std::string::npos ? "true" : "false") << std::endl;
    std::cout << "n_2_i of ACGT: " << pd.nucToInt("ACGT") << std::endl;
    auto k2 = "ACGTTTAAATGTGC";
    std::cout << "n_2_i of " << k2 << ": " << pd.nucToInt(k2) << std::endl;

    return 0;

    auto reads = read_start_end(mu_start, mu_end);
    stitch(reads);
    //perfect_stitch(reads);
    distance(reads);

    std::cout << "start and end files read\n";

    fastaReader fr;

    //std::string geno = fr.read(inFile);
    std::string geno = fr.read(inFile);
    std::cout << "genome read\n";
    //FMIndex * c22 = make_fm(inFile);

    auto snps = read_snp_file(snp_pos, snp_char);
    std::cout << "number of snps: " << snps.size();

    FMIndex *c22 = readFMFile(outFile);
    //auto hits = get_hits(reads, &geno, c22);
    assert_unique(reads, &geno, c22);
    std::cout << "hits get";

    write_start_end(reads, "mu_unique_starts", "mu_unique_ends");

    auto v_snp = check_snp_unique(reads, &geno, &snps, c22);

    std::cout << "The number of valid SNPs: " << v_snp.size() << std::endl;


    //write_hits_file("c22_mu_hits_with_ambs_0710", hits);

}
/*
    // check hits and output to file
    const std::string output = "mu_hits";
    std::ofstream ofs (output);
    if(ofs.is_open()){
        for (std::vector<int> arr : *reads){

        }
    }


    // make_fm(inFile, outFile);
    // FMIndex * c22 = make_fm(inFile);
    FMIndex *c22 = readFMFile(outFile);
    fastaReader fr;
    std::string geno = fr.unam(inFile);
    //std::string geno = fr.read(inFile);
    std::cout << "genome read\n";

    //std::vector<std::vector<int>> * reads = json_reader("/Users/josephkang/Documents/uniquekmer/output/c22_mu");
    //std::vector<std::vector<int>> *reads = just_dump_read(
    //        "/Users/josephkang/Documents/uniquekmer/src/c22_mu_just_dump");
    //std::cout << "json file read\n";
    */

    /*
    // std::vector<std::string> sequences;
    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;
    std::vector<int> hits;
    //hits.emplace_back(c22->find(matches, "GAATTCTTGTGT"));//0 to 10
    //hits.emplace_back(c22->find(matches, "ATTCTTGTGTTTAT"));//2 to 15
    //hits.emplace_back(c22->find(matches, "AATTCTTGTGTTT"));//1 to 14
    hits.emplace_back(c22->find(matches, "AATTCTTGTGTTTA"));//1 to 15
    assert("AATTCTTGTGTTTA" == geno.substr(1,15-1));
    hits.emplace_back(c22->find(matches, geno.substr(2,16-2)));
    hits.emplace_back(c22->find(matches, geno.substr(7,21-7)));
    hits.emplace_back(c22->find(matches, geno.substr(13,26-13)));
    hits.emplace_back(c22->find(matches, geno.substr(16,29-16)));
    //hits.emplace_back(c22->find(matches, "TTCTTGTGTTTAT"));//3 to 16

    //hits.emplace_back(c22->find(matches, "TTCTTGTGTTTATA"));//3 to 17
    //hits.emplace_back(c22->findn("GAATTCTTGTG"));
    //auto hits = find_matches(json_file, c22, geno);
    //std::cout << "hits counted\n";
    //std::cout << "length of hits vector: " << hits.size();
     */

/*

    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;
    std::vector<int> hits;
    //std::cout << geno.substr(10510000, 13) << ": ";

    auto starts = just_read_it("/Users/josephkang/Documents/uniquekmer/src/c22_mu_starts");
    auto ends = just_read_it("/Users/josephkang/Documents/uniquekmer/src/c22_mu_ends");

    assert(starts.size() == ends.size());
    size_t length = starts.size();

    //progress bar
    float progress = 0.0;

    std::vector<int> hits_checked;

    while(progress < 1.0){
        std::cout << "Finding all the hits:"<<std::endl;
        for (size_t i = 0; i < length; i++){

            int barWidth = 70;
            std::cout << "[";
            int pos = int(barWidth * progress);
            for (size_t j = 0; j < barWidth; ++j){
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }

            std::cout << "] " << int (progress * 100.0)  << " %\r";
            std::cout.flush();

            progress += 1.0/length;

            // end progress bar

            //assert(c22->find(matches, geno.substr(starts[i], ends[i]-starts[i])) == 1);
            hits_checked.emplace_back(c22->find(matches, geno.substr(starts[i], ends[i] - starts[i])));

        }
    }
    */
/*
    //hits.emplace_back(c22->find(matches, geno.substr(10510000, 13)));
    std::cout << geno.substr(4,17-4) << "\n";
    hits.emplace_back(c22->find(matches, geno.substr(0, 13-0)));
    hits.emplace_back(c22->find(matches, geno.substr(1, 15-1)));
    hits.emplace_back(c22->find(matches, geno.substr(2, 16-2)));
    hits.emplace_back(c22->find(matches, geno.substr(4, 17-4)));
    hits.emplace_back(c22->find(matches, geno.substr(6, 20-6)));
    hits.emplace_back(c22->find(matches, geno.substr(7, 21-7)));



    //auto hits = get_hits(reads, &geno, c22);
    //write_hits_file("c22_fmi_hits", hits);
    //std::string fn = "/Users/josephkang/CLionProjects/use_FM/output/c22_hits";
    //write_finds(reads, &geno, c22, &fn);
    for (int hit: hits) {
        std::cout << hit << "\n";

    }

    return 0;
}
*/

