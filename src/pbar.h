//
// Created by Joseph Kang on 2019-07-09.
//
// progress bar

#ifndef USE_FM_PBAR_H
#define USE_FM_PBAR_H
#include <iostream>



class pbar {

    std::string current_progress;
    size_t length;
    float progress;
    int pos;
    int barWidth;
    std::string  desc;
    int count;



public:
    pbar(size_t length, std::string  description);
    void update();
    void show();
};


#endif //USE_FM_PBAR_H
