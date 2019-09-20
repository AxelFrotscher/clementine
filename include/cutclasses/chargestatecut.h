//
// Created by afrotscher on 8/1/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH2D.h"
#include "TCutG.h"
#include "progress.h"

using std::vector, std::mutex, std::string, std::atomic;

class ccsc{
public:
    void innerloop(treereader &tree, vector<int> range);
    void analyse(const vector<string> &input, TFile* output);
    ccsc(const vector<string> &input,
         vector<vector<atomic<bool>>> &goodevents_,
         TFile* output):
         goodevents(goodevents_){
         analyse(input, output);
    };

    vector<vector<atomic<bool>>> &goodevents;

private:
    //const int threads = 7;
    int threads = std::min(25, std::max((int)sqrt(goodevents.size())/750,2));
    vector<TH2I> cschist;

    vector<TCutG*> mycut;
    mutex unitemutex;
};
