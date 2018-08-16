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

void makehistograms(std::vector<std::__cxx11::string> input);
void dalicalib(treereader *tree, TFile *output);