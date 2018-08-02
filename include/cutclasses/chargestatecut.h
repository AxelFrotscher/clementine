//
// Created by afrotscher on 8/1/18.
//

#pragma once
#include "treereader.hh"
#include "TFile.h"
#include "histogram_cuts.hh"
#include "MakeAllTree_78Ni.hh"
#include <atomic>
#include <thread>
#include <string>

class ccsc{
public:
    std::thread innerloop(treereader *tree, std::vector<std::atomic<bool>> &goodevents,
                          std::vector<int> range);
    void analyse(std::vector<treereader*> tree, TFile* output);
    ccsc(std::vector<treereader*> tree, std::vector<std::atomic<bool>> &goodevents_,
            TFile* output):goodevents(goodevents_){
            analyse(tree, output);
    };

    std::vector<std::atomic<bool>> &goodevents;

private:
    int threads;
    std::vector<TH2D> cschist;

    std::vector<TCutG*> mycut;
    std::mutex unitemutex;

};
