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
    explicit setting() = default;
    void loadnumbers(int i);
    static void setcountno(int i){
        eventcounts = i;
        checkphysicsrun();
    };
    static void checkphysicsrun();

    static const std::vector<std::vector<double>> getHOcutval();
    static const calibpar getHOparameters();
    static TCutG* getbrhocut();
    static std::vector<TCutG*> getplasticcut();
    static const std::vector<std::vector<int>> getPlasticRange();
    static const std::vector<uint> getZrange();
    static const std::vector<std::string> getreactions();
    static bool isemptyortrans(){return isemptyrun || istransmissionrun;}
    static const std::array<double, 10> getPIDincutvalue();
    static const std::array<double, 4> getPIDoutcutvalue();
    static const std::string getmodename();
    static std::string getsetname(){return setname.at(settingnumber);};
    std::vector<uint> goodruns;
    std::vector<uint> transmissionrun;
    std::vector<uint> emptyrun;
    uint analysedfile = 0;
    static int getsetnumber(){return settingnumber;};
    static void setminos(bool minosbool){minos = minosbool;};
    static bool getminos(){return minos;};
    static const std::vector<double> geticlimits();

private:
    static inline int settingnumber = 2'000'000'000;
    static inline int eventcounts = 0;
    static inline bool istransmissionrun = false;
    static inline bool isemptyrun = false;
    static inline bool minos = false;
    static inline const std::vector<std::string> setname{
                      "110Nb", "88Ge", "94Se", "100Kr", "66Cr", "70Fe", "78Ni"};
};
