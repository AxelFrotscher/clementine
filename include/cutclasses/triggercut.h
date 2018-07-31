//
// Created by afrotscher on 7/27/18.
//

#pragma once
#include "treereader.hh"
#include "TFile.h"
#include "histogram_cuts.hh"
#include "MakeAllTree_78Ni.hh"
#include <atomic>
#include <thread>


class triggercut {
public:
    std::vector<std::atomic<bool>> &goodevents;

    std::thread innerloop(treereader* tree, std::vector<std::atomic<bool>> &goodevents,
                   std::vector<int> range);
    void analyse(std::vector<treereader*> tree);
    triggercut(std::vector<treereader*> tree, std::vector<std::atomic<bool>> &goodevents_)
            : goodevents(goodevents_){
        analyse(tree);
    };

private:
    TFile *output;
    const int badtrg = 6; // triggerbit to exclude (does not contain F7DS)
    int threads;
};
