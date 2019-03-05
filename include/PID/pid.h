//
// Created by afrotscher on 8/7/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH2D.h"
#include <atomic>

class PID {
public:
    void innerloop(treereader *tree, treereader *minostree,
                   std::vector<std::vector<std::atomic<bool>>> &goodevents,
                   std::vector<uint> range);
    void analyse(const std::vector<std::string> &input, TFile* output);
    void offctrans();
    void crosssection();
    void reactionparameters();
    void histogramsetup();
    void brhoprojections();
    void chargestatecut();

    PID(const std::vector<std::string> &input, std::vector<std::vector<std::atomic<bool>>>
    &goodevents_, TFile* output, const std::string &reaction_):
    goodevents(goodevents_),reaction(reaction_){
        analyse(input, output);
    };

private:
    std::vector<std::vector<std::atomic<bool>>> &goodevents;

    int threads =1; //= std::min(25, std::max((int)sqrt(goodevents.size())/525,2));
    std::string reaction = "";

    std::vector<std::vector<TH2D>> PIDplot;
    std::vector<TH1D> reactF5;
    std::vector<std::vector<TH2D>> reactPPAC; // F7,F9,F11
    std::vector<std::vector<TH1D>> brhoprojection; // F5,7,9,11(Brho)
    std::vector<TH2D> chargestate;
    std::vector<TH2D> minosresults;

    std::vector<double> incval; // cut on incoming particles (F7)
    std::vector<double> targetval; // second cut to detected particles (F11)

    //Setup Crossection trigger:
    std::atomic<int> reactionpid1{0};
    std::atomic<int> reactionpid2{0};

    //Charge state victim count
    double chargestatevictims = 0;

    //std::vector<double> acceptancerange{-50,70}; // mm
    std::mutex unitemutex;
    int binning = 100;

    std::vector<TH2D> fitplot;
    std::vector<std::vector<TF1*>> fitstyle; // off-center 2D matrix fits

    double maxchisq = 1.15;
    TF1* bestfit = nullptr;

    double offcentertransmission = 1;
    double offcentertransmissionerror = 0;
};
