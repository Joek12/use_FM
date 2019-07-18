//
// Created by Joseph Kang on 2019-07-09.
//

#include "pbar.h"
#include <iostream>



pbar::pbar(size_t length, const std::string description) {
    this->current_progress="|";
    this->length = length;
    this->progress = 0.0;
    this->pos = -100;
    this->barWidth = 70;
    this->count = 0;

    this->desc = description;

    for (size_t j = 0; j < barWidth; ++j){
        current_progress += " ";
    }

}

void pbar::update(){
    progress += 1.0 / length;
    ++count;
    auto temp = int(barWidth * progress);
    pos = temp;

    /*
    if (temp != pos){
        if (pos >= 0) {
            current_progress = "|";
        }
        pos = temp;
    }
     */
    current_progress = "|";
    for (size_t j = 0; j < barWidth; ++j) {
        if (j < pos) current_progress += "█";
            //else if (j == pos) current_progress += "█";
        else current_progress += " ";
    }

}

void pbar::show(){
    // print description
    if (!this->desc.empty()) std::cout << this->desc << " ";

    std::cout << int(progress * 100.0) << "%" << current_progress << " | " << count << '/' << length << " \r";
    std::cout.flush();

}

