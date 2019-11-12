#pragma once
#include "helper/treereader.hh"
#include "TGraphErrors.h"

void makehistograms(std::vector<std::__cxx11::string> input);
TGraphErrors nancycs(const int &setnumber);
void writeroot(const vector<std::string> &input, const std::string &out,
               const vector<vector<std::atomic<bool>>> &goodevents,
               const bool minosbool);