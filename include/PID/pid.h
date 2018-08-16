//
// Created by afrotscher on 8/7/18.
//

#pragma once
#include "treereader.hh"
#include "TH2D.h"
#include "atomic"

class PID {
public:
    void innerloop(treereader *tree, std::vector<std::atomic<bool>>
    &goodevents, std::vector<uint> range);
    void analyse(std::vector<std::string> input, TFile* output);
    double offctrans();
    void crosssection(double transmission);
    void reactionparameters();
    void histogramsetup();

    PID(std::vector<std::string> input, std::vector<std::atomic<bool>>
    &goodevents_, TFile* output, const std::string reaction_):
    goodevents(goodevents_),reaction(reaction_){
        analyse(input, output);
    };

private:
    const int threads= 20;
    std::string reaction = "";

    std::vector<std::atomic<bool>> &goodevents;

    std::vector<std::vector<TH2D>> PIDplot;
    std::vector<TH1D> reactF5;

    std::vector<double> incval; // cut on incoming particles (F7)
    std::vector<double> targetval; // second cut to detected particles (F11)

    //Setup Crossection trigger:
    std::atomic<int> reactionpid1{0};
    std::atomic<int> reactionpid2{0};

    std::vector<double> acceptancerange{-50,70}; // mm
    std::mutex unitemutex;
    int binning = 100;

};
