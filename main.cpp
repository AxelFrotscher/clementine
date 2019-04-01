#include "TH1.h"
#include "TF1.h"
#include <boost/algorithm/string/replace.hpp>
#include "histograms.hh"
#include "MakeAllTree_78Ni.hh"
#include "zdssetting.h"
#include <iostream>
#include <algorithm>
#include <thread>
#include <fstream>
#include <vector>
#include <Math/MinimizerOptions.h>
#include "libconstant.h"
#include "TArtStoreManager.hh"
#include "TMinuitMinimizer.h"

R__LOAD_LIBRARY(libanacore.so)

using std::vector, std::string, std::ifstream, std::cout, std::endl, std::cerr,
      std::cin, std::invalid_argument, std::__throw_invalid_argument;

vector<string> getlist(const char *instring){
    // read in the file list for all .rdif's
    ifstream infile(instring);
    if(!infile) cerr << "File list not found " << endl;
    
    vector<string> file_list;
    string temp;
    while(infile >> temp){
        file_list.push_back(temp);
    }
    printf("For this analysis %lu files are available!\n", file_list.size());
    return file_list;
}

int main(int argc, char**argv){
    // Make ROOT Thread-aware
    ROOT::EnableThreadSafety();
    // Do not store TObject's in ROOT static class
    ROOT::GetROOT()->SetObjectStat(false);
    // Do not store TH1 (And derived) in ROOT static class
    TH1::AddDirectory(kFALSE);
    // Do not store TF1 in ROOT static class
    TF1::DefaultAddToGlobalList(kFALSE);
    // Set Loglevel to warning (no canvas creation notice)
    gErrorIgnoreLevel = kWarning;
    // Set minimizer to Minuit2, successor of TMinuit, threadsafe
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    // Set TMinuit not to use static TMinuit instance by gMinuit
    TMinuitMinimizer::UseStaticMinuit(kFALSE);

    // The Scope of this project is to analyse the 2015 SEASTAR
    printf("Welcome to the analysis program for 2014/2015 SEASTAR DATA\n"
           "Would you like to analyse a .ridf file first? \n"
           " (1) yes or (0) no: ");
    int analyse_raw = 0;
    if(!(cin >> analyse_raw)) throw invalid_argument("WTF");

    // Step 1: analyse the raw data
    //const uint analysedfile = 57; // 57 Index of analysed file (first, offset)
    const vector<string> input = getlist("../config/minosridf.txt");

    printf("Which Setting do you want to analyse?\n"
           "[0] 110Nb \n[1] 88Ge \n[2] 94Se \n[3] 100Kr\n[4] 66Cr\n[5] 70Fe\n"
           "[6] 78Ni\n");
    int set = 0;
    if(!(cin >> set)) throw invalid_argument("That's no setting!");
    auto s = setting(set);

    vector<string> output;
    for(auto run : s.goodruns)
        output.push_back(boost::replace_all_copy(input.at(s.analysedfile+run),
                                                 "ridf", "root"));

    vector<string> transmissionout;
    for(auto run: s.transmissionrun)
        transmissionout.push_back(boost::replace_all_copy(input.at(run),
                                                          "ridf", "root"));

    vector<string> emptyout;

    for(auto run: s.emptyrun)
        emptyout.push_back(boost::replace_all_copy(input.at(run),"ridf", "root"));

    switch(analyse_raw){
        case 1:{ // Analyse SEASTAR-DATA
            cout << "Which file to analyse?"
                    " [0] empty, [1] trans, [2] file " << endl;
            uint i = 0;
            if(!(cin >> i)) throw invalid_argument("WTH no int");

            switch(i){
                case 0:{
                    cout << "Analyzing emptyrun [0-" << s.emptyrun.size()-1 << "] "
                         << endl;
                    for(int i=0; i<s.emptyrun.size(); i++){
                        //Delete inernal static storage of master class
                        auto man = TArtStoreManager::Instance();
                        delete man;
                        generatetree(input.at(s.emptyrun.at(i)),
                                     emptyout.at(i));
                    }
                    break;
                }
                case 1:{
                    cout << "(Trans) Analyzing SEASTAR:"
                         << input.at(s.transmissionrun.at(0)) << endl;
                    generatetree(input.at(s.transmissionrun.at(0)),
                                 transmissionout.at(0));
                    break;
                }
                case 2:{
                    cout << "Which physicsrun? [0-"<< s.goodruns.size()-1
                         <<"], [" << s.goodruns.size() << "] = all" << endl;
                    if(!(cin >> i)) throw invalid_argument("WTH");
                    if(i == s.goodruns.size()){
                        cout << "Analyze ALL the events" << endl;
                        for(int j=0; j<s.goodruns.size(); j++){
                            //Delete inernal static storage of master class
                            auto man = TArtStoreManager::Instance();
                            delete man;
                            generatetree(input.at(s.analysedfile+s.goodruns.at(j)),
                                      output.at(j));
                        }
                        break;
                    }
                    cout << "(Empty) Analyzing SEASTAR:"
                         << input.at(s.analysedfile+s.goodruns.at(i)) << endl;
                    {generatetree(input.at(s.analysedfile+s.goodruns.at(i)),
                                 output.at(i));}
                    break;
                }
                default: __throw_invalid_argument("Option not available\n");
            }
            break;
        }
        case 0:{
            int choose;
            printf("Now proceeding to make histograms...\n");
            cout << "Which setting? [0 physics, 1 trans, 2 empty] " << endl;
            if(!(cin >> choose)) throw invalid_argument("WTH no int");
            switch(choose){
                case 0:{
                    //Add Minos choice
                    bool minoschoice = false;
                    cout << "MINOS Analysis? [0] no, [1] yes " << endl;
                    if(!(cin >> minoschoice)) throw invalid_argument("WTH!");
                    s.setminos(minoschoice);
                    makehistograms(output);
                    break;
                }
                case 1:{
                    makehistograms(transmissionout);
                    break;
                }
                case 2:{
                    makehistograms(emptyout);
                    break;
                }
                default: __throw_invalid_argument("Option not available\n");
            }
            break;
        }
        default:
            __throw_invalid_argument("Option not available");
    }
    return 0;
}