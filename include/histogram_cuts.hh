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
#include "TCutG.h"

struct calibpar{
    calibpar(){
    }
    calibpar(double d, double e, double f, double g, double h, double i,
             double j, double k, double l, double m, double n) :
            F7absF5X(d), F7linF5X(e), F7linF5A(f), F7linF3X(g), F7absF5X0(h),
            F11absF9X(i), F11linF9X(j),F11linF9A(k), F11linF11X(l), F11linF11A(m),
            F11absF9X0(n){
    }
    // Nice encapsulation of variables for the correction
    double F7absF5X;
    double F7linF5X;
    double F7linF5A;
    double F7linF3X;
    double F7absF5X0;

    double F11absF9X;
    double F11linF9X;
    double F11linF9A;
    double F11linF11X;
    double F11linF11A;
    double F11absF9X0;
};

// Definition of global variables here. Keep as short as possible!
extern calibpar p1;
extern std::mutex consolemutex;
extern std::mutex writemutex;

double linfit(double *x, double *par);
const bool closeness(const std::vector<double> &d, double sigma=0.1);
void plastics(treereader *tree, TFile *output, std::vector<std::atomic<bool>> &goodevents);
void ionisationchamber(treereader *alt2dtree, TFile *output,
                       std::vector<std::atomic<bool>> &goodevents);
void targetcut(treereader *tree, TFile *output, std::vector<std::atomic<bool>> &goodevents);