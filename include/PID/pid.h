//
// Created by afrotscher on 8/7/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH2D.h"
#include "TH3D.h"
#include <atomic>
#include "TGraphErrors.h"

class PID {
public:
    void innerloop(treereader &tree, treereader &minostree,
                   const std::vector<int> &range, bool minosanalyse);
    void analyse(const std::vector<std::string> &input, TFile* output);
    void offctrans();
    void crosssection();
    void reactionparameters();
    void histogramsetup();
    void brhoprojections();
    void chargestatecut();

    PID(const std::vector<std::string> &input, std::vector<std::vector<std::atomic<bool>>>
    &goodevents_, TFile* output, const std::string &reaction_, TGraphErrors &tcross_):
    goodevents(goodevents_),reaction(reaction_), tcross(tcross_){
        analyse(input, output);
    };

private:
    std::vector<std::vector<std::atomic<bool>>> &goodevents;

    int threads = std::min(25, std::max((int)sqrt(goodevents.size())/600,2));
    std::string reaction = "";

    std::vector<TH1F> reactF5;
    std::vector<std::vector<TH2S>> reactPPAC; // F7,F9,F11
    std::vector<std::vector<TH1F>> brhoprojection; // F5,7,9,11(Brho)
    std::vector<TH2I> chargestate;
    std::vector<decltype(chargestate)> PIDplot;  // "Uncorrected/corrected"
    std::vector<TH2I> minosresults;
    std::vector<TH3S> minos3dresults;
    std::vector<TH2C> minossingleevent; // Contains histograms for single events
    std::vector<TH1I> minos1dresults;

    std::array<double ,10> incval; // cut on incoming particles (F7)
    std::array<double, 4> targetval; // second cut to detected particles (F11)

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

    //TGraph to store all cross sections
    TGraphErrors &tcross;
    int projectileN = 0;
    int ncross =0;
};
