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
#include <atomic>

using namespace std;

void makehistograms(const std::vector<std::__cxx11::string> input);
void makepid(const vector<string> input, TFile* output,
             const std::vector<std::atomic<bool>> &goodevents);
void dalicalib(treereader *tree, TFile *output);
void pidth(treereader *tree, vector<vector<TH2D>> &PID, const vector<double> &cutval,
           const vector<double> &targetval, vector<uint> &reactionval,
           const vector<uint> &range, const vector<std::atomic<bool>> &goodevents,
           const int id, vector<TH1D> &reactf9);