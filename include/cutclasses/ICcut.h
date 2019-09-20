//
// Created by afrotscher on 8/2/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH2D.h"

class iccut{
public:
    void innerloop(treereader &tree, std::vector<int> range);
    void analyse(const std::vector<std::string> &input, TFile* output);
    iccut(const std::vector<std::string> &input, std::vector<std::vector<std::atomic<bool>>>
            &goodevents_, TFile* output):goodevents(goodevents_){
            analyse(input, output);
    };

    std::vector<std::vector<std::atomic<bool>>> &goodevents;

private:
    int threads = std::min(25, std::max((int)sqrt(goodevents.size())/750,2));
    std::mutex unitemutex;

    const int numchannel = 6; // There are 6 channel per IC
    //const int numic = 2; // Number of Ionisation Chambers
    std::vector<std::vector<TH2I>> comparediag;

};

