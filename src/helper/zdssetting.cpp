//
// Created by afrotscher on 9/26/18.
//

#include <iostream>
#include <TFile.h>
#include "zdssetting.h"
#include "libconstant.h"

using std::vector;

int setting::settingnumber = 2'000'000'000;
int setting::eventcounts = 0;
bool setting::istransmissionrun = false;
bool setting::isemptyrun = false;
const vector<std::string> setting::setname{"110Nb", "88Ge", "94Se", "100Kr", "66Cr", "70Fe", "78Ni"};

void setting::loadnumbers(int i) {
    // Loading the right numbers for th right setting
    std::cout << "Setting number " << i << std::endl;
    switch(i){
        case 0:{ // 110Nb
            analysedfile = 57; // 57 Index of analysed file (first, offset)
            goodruns = vector<uint>{1,3,4,5,7,8,9,10,11,13,14,19,26,29,32,33,35,
                                    37, 38,40,41,42,46,49,50,51,52,54,57,58,59};
            transmissionrun = vector<uint>{54};
            emptyrun = vector<uint>{52,53};
            break;
        }
        case 1:{ // 88Ge
            analysedfile = 119;
            goodruns = vector<uint>{0,3,8,9,10,11,12};
            transmissionrun = vector<uint>{118};
            emptyrun = vector<uint>{}; // 119 is not an empty run
            break;
        }
        case 2:{ // 94Se
            analysedfile = 131;
            goodruns = vector<uint>{27,28,29,30,31, 4,5,8,19,21,22,23,24,25,26,
                                    32,33,34,35,36,39,40,45,46,47}; // 48
            transmissionrun = vector<uint>{};
            emptyrun = vector<uint>{132};
            break;
        }
        case 3:{ // 100Kr
            analysedfile = 180;
            goodruns = vector<uint>{4,6,7,8,9,10,11,12,13,14,15,17,19,20,16,21,
                                    22,23,24,25,26,27,28,29,30,31,32,33,34,35,
                                    36,38};
            transmissionrun = vector<uint>{183};
            emptyrun = vector<uint>{182};
            break;
        }
        case 4:{ //66Cr
            analysedfile = 277;
            goodruns = vector<uint>{38,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
                                    54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,
                                    69};
            transmissionrun = vector<uint>{299};
            emptyrun = vector<uint>{279};
            break;
        }
        case 5:{ // 70Fe
            analysedfile = 346;
            goodruns = vector<uint>{7,8,9,11,12,13,14,15,18,19,20,21,22,23,24,25,
                                    26,27,28,30,31,32,33,34,35,36,37,39,40,41,42,
                                    43,44,45};
            transmissionrun = vector<uint>{};
            emptyrun = vector<uint>{};
            break;
        }
        case 6:{
            analysedfile = 395;
            goodruns = vector<uint>{8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
                                    24,25,26,28,29,30,31,34,35,36,37,38,39,40,41,
                                    42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,
                                    57,58,59,62,63,64,65,66,67,68,69,70,71,72,73,
                                    74,75,76,77,78,79,81,82,83,84,85,80,90,91,93,
                                    94,95,96,97,98,99,100,101,102,103,104,105,106,
                                    108};
            for(int i=110; i<156; i++) goodruns.push_back(i);
            transmissionrun = vector<uint>{394};
            emptyrun = vector<uint>{393};
            break;
        }
        default: std::__throw_invalid_argument("Chosen setting is not implemented.\n");
    }
}

void setting::checkphysicsrun() {
    // Determine type of analysed run
    if(eventcounts == runinfo::transsize.at(settingnumber)){
        std::cout << "Analysing Transmission run Setting "
                  << settingnumber << std::endl;
        istransmissionrun = true;
        return;
    }
    else if (eventcounts == runinfo::emptysize.at(settingnumber)){
        std::cout << "Analysing Empty Target run Setting "
                  << settingnumber << std::endl;
        isemptyrun = true;
        return;
    }
    else if(eventcounts == runinfo::fulldata.at(settingnumber)){
        std::cout << "Analysing full physics run Setting "
                  << settingnumber << std::endl;
        return;
    }

    std::__throw_invalid_argument("Run/Setting does not exist!\n");
}

const std::vector<std::vector<double>> setting::getHOcutval() {
    // Get the right cut values to perform a linear correction
    if (istransmissionrun) return nancytrans::cutval.at(settingnumber);
    if (isemptyrun) return nancyempty::cutval.at(settingnumber);
    return nancy::cutval.at(settingnumber);
}

const calibpar setting::getHOparameters() {
    // Get the right linear parameters
    if (istransmissionrun) return nancytrans::hoparame.at(settingnumber);
    if (isemptyrun) return nancyempty::hoparame.at(settingnumber);
    return nancy::hoparame.at(settingnumber);
}

TCutG* setting::getbrhocut() {
    // Preparing the right cut for the brho measurements CCSC
    auto nam = [](std::string str){
        return "brhocut" + str + setname.at(settingnumber);};

    // Get Cut
    TFile cutfile("config/cut.root");
    if(!(cutfile.IsOpen())) std::__throw_invalid_argument(
            "Could not open cut file at config/cut.root\n");

    if (istransmissionrun && ((TCutG*)cutfile.Get(nam("trans").c_str()) != nullptr)){
        return (TCutG*)cutfile.Get(nam("trans").c_str());
    }
    if (isemptyrun && ((TCutG*)cutfile.Get(nam("empty").c_str()) != nullptr)){
        return (TCutG*)cutfile.Get(nam("empty").c_str());
    }
    if ((TCutG*)cutfile.Get(nam("data").c_str()) != nullptr)
        return (TCutG*)cutfile.Get(nam("data").c_str());

    std::__throw_invalid_argument("Could not load cut from file!\n");
}

const vector<vector<int>> setting::getPlasticRange() {
    // Getting the right Parameter
    return runinfo::plasticrange.at(settingnumber);
}

const vector<uint> setting::getZrange() {
    // Getting the correct Z Range for the settings
    return runinfo::pidZrange.at(settingnumber);
}

const vector<std::string> setting::getreactions() {
    // Determine the reactions to calculate for the settings
    return runinfo::reactionmodes.at(settingnumber);
}

const vector<double> setting::getPIDincutvalue() {
    // return cut particle for incoming beam (trans/empty only)
    if(!setting::isemptyortrans()) std::__throw_invalid_argument("Run in not empty or transmision run!\n");

    if(istransmissionrun) return nancytrans::incval.at(settingnumber);
    if(isemptyrun) return nancyempty::incval.at(settingnumber);
}

const vector<double> setting::getPIDoutcutvalue() {
    // return cut particle for incoming beam (trans/empty only)
    if(!setting::isemptyortrans()) std::__throw_invalid_argument("Run in not empty or transmision run!\n");

    if(istransmissionrun) return nancytrans::targetval.at(settingnumber);
    if(isemptyrun) return nancyempty::targetval.at(settingnumber);
}

const std::string setting::getmodename() {
    if(istransmissionrun) return "transmission";
    if(isemptyrun) return "empty";
    return "data";
}