#include <iostream>
#include <fstream>
#include "FMCommands.h"
#include "fastaReader.h"
#include "pbar.h"
#include "mu_commands.h"
#include "primer_design.h"
#include "../test/test_primer_design.h"

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
    //test_distribution();
    //test_nucToInt();
    //test_gc_content();
    test_optimal();
    return 0;
}
