//
// Created by axel on 05.06.18.
//

#pragma once
#include "TString.h"
#include "TStreamerElement.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TMath.h"
#include "TROOT.h"
#include "treereader.hh"
#include "TProfile.h"
#include "TF1.h"

struct calibpar{
    // Nice encapsulation of variables for the correction
    double F7absF5X = 0;
    double F7linF5X = 0;
    double F7linF5A = 0;
    double F7linF3X = 0;
    double F7absF5X0 = 0;

    double F11absF9X = 0;
    double F11linF9X = 0;
    double F11linF9A = 0;
    double F11linF11X = 0;
    double F11linF11A = 0;
    double F11absF9X0 = 0;
};

// Definition of global variables here. Keep as short as possible!
extern calibpar p1;

double linfit(double *x, double *par);
const bool closeness(const std::vector<double> &d, double sigma=0.1);
void plastics(treereader &tree, TFile &output, std::vector<bool> &goodevents);
void ppacs(treereader &tree, TFile &output, std::vector<bool> &goodevents);
void ionisationchamber(treereader &alt2dtree, TFile &output,
                       std::vector<bool> &goodevents);
void chargestatecut(treereader &tree, TFile &output,
                    std::vector<bool> &goodevents);