//
// Created by afrotscher on 7/30/18.
//

#pragma once
#include "helper/treereader.hh"
#include "TH1D.h"
#include "TH2D.h"

class ppaccut{
public:
    void innerloop(treereader &tree, std::vector<std::vector<std::atomic<bool>>> &goodevents,
                          std::vector<uint> range);
    void analyse(std::vector<std::string> input, TFile* output);
    ppaccut(const std::vector<std::string> input, std::vector<std::vector<std::atomic<bool>>> &goodevents_,
    TFile* output): goodevents(goodevents_){
        analyse(input, output);
    };

    std::vector<std::vector<std::atomic<bool>>> &goodevents;

private:
    int threads = std::min(25, std::max((int)sqrt(goodevents.size())/750,2));
    const int numplane = 36;
    const int pl11position = 3; //Plastic at F11 is fourth in array
    std::vector<TH1F> effPPAC;
    std::vector<std::vector<TH2I>> sumdiffppac; // 1 X,Y 2 Sum,diff (2D NoPPAC, Quantity)

    const std::vector<std::vector<int>> ppacplane{
        {4,5,6,7},{9,10,11,12},{14,15,16,17},{18,19,20,21},{22,23,24,25},{30,31,32,33}};

    const std::vector<std::vector<std::vector<int>>> ppacrange{ //{xlow,xup,ylow,yup}
        {{85,105,85,100},{85,100,85,100},{165,180,85,105},{170,185,90,105}},  //Fpl 3
        {{165,185,85,105},{165,185,90,105},{165,175,90,105},{165,180,95,110}},//Fpl 5
        {{160,175,95,110},{165,180,85,100},{75,90,90,100},{85,105,95,110}},   //Fpl 7
        {{170,185,99,110},{165,185,95,110},{165,185,95,110},{165,185,90,105}},//Fpl 8
        {{130,145,60,75},{135,150,60,75},{140,155,70,80},{130,145,50,65}},    //Fpl 9
        {{145,160,65,80},{135,150,55,70},{145,160,75,90},{135,150,60,75}}     //Fpl11
    };

    std::atomic<int> efficiencytotal{0};
    std::mutex unitemutex;

};

