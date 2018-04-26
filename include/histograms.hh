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

void makehistograms(const std::__cxx11::string input);
bool closeness(std::vector< double >& d, double sigma);
void makepid(const TTreeReader& datree, TFile& output);