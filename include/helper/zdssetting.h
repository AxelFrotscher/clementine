//
// Created by afrotscher on 9/26/18.
//

#pragma once

#include <vector>
#include <zconf.h>
#include "libconstant.h"
#include "TCutG.h"

class setting{
public:
    // This class aims to deliver the right numbers for each setting, by just
    // setting it once
    explicit setting(int i){
        settingnumber = i;
        loadnumbers(i);
    }
    explicit setting(){;};
    void loadnumbers(int i);
    void setcountno(int i){
        eventcounts = i;
        checkphysicsrun();
    };
    void checkphysicsrun();

    const std::vector<std::vector<double>> getHOcutval();
    const calibpar getHOparameters();
    TCutG* getbrhocut();
    const std::vector<std::vector<int>> getPlasticRange();
    const std::vector<uint> getZrange();
    const std::vector<std::string> getreactions();
    const bool isemptyortrans(){return isemptyrun || istransmissionrun;}
    const std::vector<double> getPIDincutvalue();
    const std::vector<double> getPIDoutcutvalue();
    const std::string getmodename();
    std::string getsetname(){return setname.at(settingnumber);};
    std::vector<uint> goodruns;
    std::vector<uint> transmissionrun;
    std::vector<uint> emptyrun;
    uint analysedfile = 0;

private:
    static int settingnumber;
    static int eventcounts;
    static bool istransmissionrun;
    static bool isemptyrun;
    static const std::vector<std::string> setname;
};
