//
// Created by afrotscher on 8/23/18.
//

#include <fstream>
#include "txtwriter.h"
#include <chrono>

std::string txtwriter::filename;
std::vector<std::string> txtwriter::outputbuffer;
std::chrono::system_clock::time_point txtwriter::begin;

void txtwriter::addline(std::string line) {
    // Perform checks on line

    //add it to the queue
    outputbuffer.push_back(line);
}

void txtwriter::writetofile(){
    // Generating output

    auto end = std::chrono::system_clock::now();

    std::ofstream out(filename);
    out << "Analysis Report SEASTAR2: "<< filename << std::endl << std::endl;
    out << "Number of Modes analysed: " << outputbuffer.size()-1 << std::endl;
    out << "Analyse Time: " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count()
        << "s "<< std::endl << std::endl;
    for(auto &i : outputbuffer) out << i << std::endl;
    out.close();
}