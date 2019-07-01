//
// Created by Joseph Kang on 2019-06-27.
//

#include <iostream>
#include "fastaReader.h"
#include <fstream>

fastaReader::fastaReader() {
    this->filename = "";
    this->sequence = "";
}

fastaReader::fastaReader(std::string fn) {
    this->filename = fn;
    this->sequence = this->read(fn);
}

std::string fastaReader::read(std::string fn) {

    this->filename = filename.empty() ? fn : filename;
    std::ifstream input (fn);
    std::string line;

    std::string DNA_sequence = "";

    while (std::getline(input, line)) {

        // line may be empty so you *must* ignore blank lines
        // or you have a crash waiting to happen with line[0]
        if (line[0] == '>' || line.empty())
            continue;

        else
            DNA_sequence += line;

    }

    return DNA_sequence;

}

std::string fastaReader::unam(std::string fn){

    this->filename = filename.empty() ? fn : filename;
    std::ifstream input(fn);
    std::string line;

    std::string unams = "";

    while (std::getline(input, line)){

        // line may be empty so you *must* ignore blank lines
        // or you have a crash waiting to happen with line[0]
        if (line[0] == '>' || line.empty())
            continue;

        // handle ambs
        else if (line.find('N') != std::string::npos){

            for (char c : line){

                // traverse each char in the line read and only append unams
                if (c != 'N' && c != '\n')
                    unams += c;
            }
        }

        else
            unams += line;
    }

    return unams;
}


std::string fastaReader::unam() {
    assert(!this->sequence.empty());

    std::string unams = "";

    if (this->sequence.find('N') != std::string::npos){
        for (char c : this->sequence){
            if (c != 'N')
                unams += c;
        }

        return unams;
    }

    else
        return this->sequence;

}
