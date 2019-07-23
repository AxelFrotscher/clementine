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
        filename = std::move(filename_);
        begin = std::chrono::system_clock::now();
    };

    explicit txtwriter() = default;
    static void addline(const std::string &line);
    static void writetofile();

private:
    inline static std::string filename;
    inline static std::vector<std::string> outputbuffer;
    inline static std::chrono::system_clock::time_point begin;
};