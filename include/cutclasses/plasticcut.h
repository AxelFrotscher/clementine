//
// Created by afrotscher on 8/2/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "zdssetting.h"

class plasticcut {
public:
    void innerloop(treereader &tree, std::vector<int> range);
    void analyse(const std::vector<std::string> &input, TFile* output);
    plasticcut(const std::vector<std::string> &input, std::vector<std::vector<std::atomic<bool>>>
               &goodevents_, TFile* output):goodevents(goodevents_){
        acceptance_range = setting::getPlasticRange();
        analyse(input, output);
    };


private:
    std::vector<std::vector<std::atomic<bool>>> &goodevents;
    int threads = std::min(25, std::max((int)sqrt(goodevents.size())/750,2));
    std::mutex unitemutex;

    std::vector<std::vector<TCutG*>> mycut;

    // generate output diagrams
    std::vector<TH2I> qcorr2D;  // Charge Correlation between 1 & 2
    std::vector<TH1I> qcorr;    // deposited charge at detector
    std::vector<TH2I> tqcorr2D; // Charge-Time correlation
    std::vector<TH1I> qxF11;
    std::vector<std::vector<std::string>> arrayname, arraytitle;

    int numplastic = 4; // 4
    const int F7pos = 1;
    const int threshhold = 450; // random noise suppresion value

    // To avoid multihit-triggers we define a range of acceptance
    std::vector<std::vector<int>> acceptance_range;

};

