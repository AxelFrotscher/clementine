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

class plasticcut {
public:
    void innerloop(treereader *tree, std::vector<std::atomic<bool>>
                          &goodevents, std::vector<int> range);
    void analyse(std::vector<std::string> input, TFile* output);
    plasticcut(const std::vector<std::string> input, std::vector<std::atomic<bool>>
               &goodevents_, TFile* output):goodevents(goodevents_){
            analyse(input, output);
    };

    std::vector<std::atomic<bool>> &goodevents;

private:
    const int threads = 7;
    std::mutex unitemutex;

    // generate output diagrams
    std::vector<TH2D> qcorr2D;  // Charge Correlation between 1 & 2
    std::vector<TH1D> qcorr;    // deposited charge at detector
    std::vector<TH2D> tqcorr2D; // Charge-Time correlation
    std::vector<std::vector<std::string>> arrayname, arraytitle;

    const int numplastic = 4;
    const int F7pos = 1;
    const int threshhold = 450; // random noise suppresion value

    // To avoid multihit-triggers we define a range of acceptance
    std::vector<std::vector<int>> acceptance_range{{480,620},{700,920},
                                                   {220,330},{270,1510}};

};

