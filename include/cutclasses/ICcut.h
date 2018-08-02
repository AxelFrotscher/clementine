//
// Created by afrotscher on 8/2/18.
//

#pragma once
#include "treereader.hh"
#include "TFile.h"
#include "histogram_cuts.hh"
#include "MakeAllTree_78Ni.hh"
#include <atomic>
#include <thread>
#include <string>

class iccut{
public:
    std::thread innerloop(treereader *tree, std::vector<std::atomic<bool>>
    &goodevents, std::vector<int> range);
    void analyse(std::vector<treereader*> tree, TFile* output);
    iccut(std::vector<treereader*> tree, std::vector<std::atomic<bool>>
            &goodevents_, TFile* output):goodevents(goodevents_){
            analyse(tree, output);
    };

    std::vector<std::atomic<bool>> &goodevents;

private:
    int threads;
    std:: mutex unitemutex;

    const int numchannel = 6; // There are 6 channel per IC
    const int numic = 2; // Number of Ionisation Chambers
    std::vector<std::vector<TH2D>> comparediag;

};

