#pragma once
#include <string>
#include <vector>
#include <iostream>
#include "TArtEventStore.hh"

void generatetree(std::__cxx11::string infile, std::__cxx11::string output);
//void progressbar(int currevent, int totevent, int offset, int barwidth=30);

uint getset(TArtEventStore &estore);