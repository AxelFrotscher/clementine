//
// Created by afrotscher on 7/27/18.
//

#pragma once
#include "treereader.hh"       // contains already many classes

class triggercut {
public:
    std::vector<std::atomic<bool>> &goodevents;

    void innerloop(treereader* tree, std::vector<std::atomic<bool>> &goodevents,
                   std::vector<uint> range);
    void analyse(std::vector<std::string> input);
    triggercut(const std::vector<std::string> input, std::vector<std::atomic<bool>> &goodevents_)
            : goodevents(goodevents_){
        analyse(input);
    };

private:
    //TFile *output;
    const int badtrg = 6; // triggerbit to exclude (does not contain F7DS)
    const int threads = 10;
};
