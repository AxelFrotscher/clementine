//
// Created by afrotscher on 8/23/18.
//

#pragma once
#include <iostream>
#include <vector>
#include <chrono>

class txtwriter{
public:
    explicit txtwriter(std::string filename_){
        filename = filename_;
        begin = std::chrono::system_clock::now();
    };

    explicit txtwriter() = default;
    void addline(std::string line);
    void writetofile();

private:
    static std::string filename;
    static std::vector<std::string> outputbuffer;

    static std::chrono::system_clock::time_point begin;
};