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

class targetcut{
public:
    std::thread innerloop(treereader *tree, std::vector<std::atomic<bool>>
    &goodevents, std::vector<int> range);
    void analyse(std::vector<treereader*> tree, TFile* output);
    targetcut(std::vector<treereader*> tree, std::vector<std::atomic<bool>>
            &goodevents_, TFile* output):goodevents(goodevents_){
            analyse(tree, output);
    };

    std::vector<std::atomic<bool>> &goodevents;
private:
    int threads;
    std::mutex unitemutex;

    const int ppacoffset = 18; // [18] is F8PPAC-1A
    const int f8ppac = 4;
    const int  magnum = -9999;
    const int f8offset = -2; // mm
    const int targetradius = 20; //mm

    // F8 PPAC positions relative to F8 0:PPAcno 1: z(x),z(y)
    const std::vector<std::vector<double>> f8z{{-1310.7,-1302.1},{-1273.3,-1281.9},
                                     {-810.7,-802.1},{-773.3,-781.9}};

    std::vector<TH2D> tarhist;
};