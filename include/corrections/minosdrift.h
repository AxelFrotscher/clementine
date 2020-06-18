//
// Created by afrotscher on 12.04.19.
//

#pragma once

#include <vector>
#include <string>
#include <array>

class minosdrift{
public:
    explicit minosdrift(const std::vector<std::string> &input) {
        make_drift(input);
    }

private:
    void make_drift(const std::vector<std::string> &input);

    std::vector<double> gettimeborders(const std::string &i);

    const float tpclength = 300; //mm
};