//
// Created by afrotscher on 8/2/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH2D.h"

class higherorder{
public:
    void innerloop(treereader *tree, std::vector<std::atomic<bool>>
    &goodevents, std::vector<uint> range);
    void analyse(std::vector<std::string> input, TFile* output);
    higherorder(const std::vector<std::string> &input, std::vector<std::atomic<bool>>
            &goodevents_, TFile* output):goodevents(goodevents_){
            analyse(input, output);
    };

    std::vector<std::atomic<bool>> &goodevents;

private:
    const int threads = 7;
    std:: mutex unitemutex;

    const std::vector<int> beam{0,4}; // Evaluate Beam F3-7 (1st element)
    const int corrcount = 3;     // Number of corrections
    const std::vector<std::vector<std::string>> arrname = {
            { "DepF3X","DepF3A","DepF5X", "DepF5A", "DepBetaF3-7"},
            { "DepF9X", "DepF9A","DepF11X", "DepF11A", "DepBetaF8-11"}};

    const std::vector<std::vector<std::string>> arrtitle = {
            {"Dependence of F3X vs AoQ", "Dependence of F3A vs AoQ",
                    "Dependence of F5X vs AoQ", "Dependence of F5A vs AoQ",
                    "Dependence of #beta (F7) vs AoQ"},
            {"Dependence of F9X vs AoQ", "Dependence of F9A vs AoQ",
                    "Dependence of F11X vs AoQ", "Dependece of F11A vs AoQ",
                    "Dependence of #beta (F11) vs AoQ"}};

    std::vector<std::vector<std::vector<TH2D>>> culpritdiag;
    std::vector<std::vector<double>> cutval;

    // Get linear fits of the projected means
    std::vector<std::vector<std::vector<TProfile *>>> projections;
    const double cutfrac = 0.8; // fraction to consider for fit

    const std::vector<std::vector<std::string>> folders{
            {"Corrections/Pre/Raw",
             "Corrections/Pre/F5Xcorr",
             "Corrections/Pre/F5XF5Acorr",
             "Corrections/Pre/F5XF5AF3X"},
            {"Corrections/Post/Raw",
             "Corrections/Post/F9Xcorr",
             "Corrections/Post/F9Acorr",
             "Corrections/Post/F11Acorr"}};
};
