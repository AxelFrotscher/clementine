//
// Created by afrotscher on 9/26/18.
//

#pragma once

#include <vector>
#include <zconf.h>

class setting{
public:
    setting(int i){
        loadnumbers(i);
    }
    void loadnumbers(int i);
    std::vector<uint> goodruns;
    std::vector<uint> transmissionrun;
    std::vector<uint> emptyrun;
    uint analysedfile = 0;

private:

};
