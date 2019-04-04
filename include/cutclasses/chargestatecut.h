//
// Created by afrotscher on 8/1/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH2D.h"
#include "TCutG.h"
#include "progress.h"

class ccsc{
public:
    void innerloop(treereader *tree, std::vector<std::vector<std::atomic<bool>>> &goodevents,
                          std::vector<uint> range);
    void analyse(std::vector<std::string> input, TFile* output);
    ccsc(const std::vector<std::string> input, std::vector<std::vector<std::atomic<bool>>> &goodevents_,
            TFile* output):goodevents(goodevents_){
            analyse(input, output);
    };

    std::vector<std::vector<std::atomic<bool>>> &goodevents;

private:
    //const int threads = 7;
    int threads = std::min(25, std::max((int)sqrt(goodevents.size())/750,2));
    std::vector<TH2I> cschist;

    std::vector<TCutG*> mycut;
    std::mutex unitemutex;

};
