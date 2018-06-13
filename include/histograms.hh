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

void makehistograms(const std::vector<std::__cxx11::string> input);
void makepid(treereader& datree, TFile* output,
             const std::vector<bool> &goodevents);
void highordercorrection(treereader &tree, TFile *output,
                         const std::vector<bool> &goodevents);
void dalicalib(treereader &tree, TFile *output);