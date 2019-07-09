//
// Created by Joseph Kang on 2019-07-09.
//

#include "pbar.h"
#include <iostream>



pbar::pbar(size_t length, const std::string description) {
    this->current_progress="[";
    this->length = length;
    this->progress = 0.0;
    this->pos = 0;
    this->barWidth = 70;

    this->desc = description;
}

void pbar::update(){
    progress += 1.0 / length;
    pos = int(barWidth * progress);
    current_progress = "[";
    for (size_t j = 0; j < barWidth; ++j){
        if (j < pos) current_progress += "=";
        else if (j == pos) current_progress += ">";
        else current_progress += " ";
    }

}

void pbar::show(){
    // print description
    if (!this->desc.empty()) std::cout << this->desc << " ";

    std::cout << current_progress << "]" << int(progress * 100.0) << "%\r";
    std::cout.flush();

}

