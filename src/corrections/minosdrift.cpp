//
// Created by afrotscher on 12.04.19.
//

#include <iostream>
#include <fstream>
#include "minosdrift.h"
#include "TFile.h"
#include "TH1I.h"
#include "TF1.h"
#include "zdssetting.h"
#include <sstream>

using std::string, std::vector, std::cout, std::endl, std::to_string;

void minosdrift::make_drift(const std::vector<std::string> &input) {
    /// First get all associated root files
    auto convert = [](auto &i){
        return "output/MINOS/" + i.substr(34,14);
    };

    vector<string> minosfiles;
    minosfiles.reserve(input.size());
    for(auto &i:input) minosfiles.push_back(convert(i));

    ///Call fit method
    vector<vector<double>> driftresults;
    driftresults.reserve(minosfiles.size());

    for(auto &i: minosfiles) driftresults.push_back(gettimeborders(i));

    /// get mean start (only depending on electronics -> one value)
    double mean_t_start = 0;
    for(auto &i: driftresults) mean_t_start += i.at(0)/driftresults.size();

    /// Replace values with mean
    for(auto &i:driftresults){
        i.at(0) = mean_t_start;
        i.push_back(tpclength/(i.at(1)-i.at(0)));
    }

    ///Generate strings to write
    vector<string> resultio;
    for(int i=0; i<input.size(); i++){
        resultio.push_back(input.at(i).substr(34,9) + " " +
                           to_string(driftresults.at(i).at(0)) + " " +
                           to_string(driftresults.at(i).at(1)) + " " +
                           to_string(driftresults.at(i).at(2)));
        cout << resultio.back() << endl;
    }
    cout << "Done!" << endl;

    /// Read in all results generated so far
    vector<string> filelocation{  // 7 locations for the seven settings
        "../config/db/ConfigMINOSDrift_fixed.txt",
        "../config/db/ConfigMINOSDrift_fixed.txt",
        "../config/db/ConfigMINOSDrift_fixed.txt",
        "../config/db/ConfigMINOSDrift_fixed.txt",
        "../config/db2014/ConfigMINOSDrift.txt",
        "../config/db2014/ConfigMINOSDrift.txt",
        "../config/db2014/ConfigMINOSDrift.txt",
    };


    std::ifstream preresults(filelocation.at(setting::getsetnumber()));
    assert(preresults.is_open());
    string dummy;
    vector<string> prestrings;
    while(std::getline(preresults, dummy)) prestrings.push_back(dummy);
    preresults.close();

    /// Loop over all previous results
    for(auto &i: prestrings){
        std::stringstream ss(i);
        string temp;
        ss >> temp; // get the runnumber
        for(int j=0; j<minosfiles.size(); j++){ // loop over all generated results to see matchs
            if(minosfiles.at(j).find(temp) != string::npos){
                i = resultio.at(j);   // if it matches replace results, and delete them from vector
                cout << "Replacing Run: " << minosfiles.at(j).substr(13,9) << endl;
                resultio.erase(resultio.begin()+j);
                minosfiles.erase(minosfiles.begin()+j);
                break;
            }
        }
    }
    /// Append the new runs
    prestrings.insert(prestrings.end(), resultio.begin(), resultio.end());
    std::sort(prestrings.begin()+3, prestrings.end());

    /// Write out newly generated file
    std::ofstream hi(filelocation.at(setting::getsetnumber()));
    for(auto &i: prestrings) hi << i << endl;
    hi.close();


}

vector<double> minosdrift::gettimeborders(const std::string &i) {
    /// Open an individual file to get the result.

    TFile tpctimefile(i.c_str());
    assert(tpctimefile.IsOpen());

    auto tpctimehist = (TH1I*)tpctimefile.Get(i.substr(13,9).c_str());
    assert(tpctimehist);

    // WE derive the stupid histogram
    TH1I derived(*tpctimehist);
    derived.Reset();
    for(int i=55; i<derived.GetNbinsX()-1; i++){
        derived.SetBinContent(i,(-tpctimehist->GetBinContent(i+2)+
                               8*tpctimehist->GetBinContent(i+1)-
                               8*tpctimehist->GetBinContent(i-1)+
                                 tpctimehist->GetBinContent(i-2))/12.);
    }
    tpctimefile.Close();

    /// Get maximum noise bin
    double maxnoise =0;
    for(int i=60; i<120; i++) maxnoise = std::max(maxnoise, derived.GetBinContent(i));
    double t_min = derived.GetBinCenter(derived.FindFirstBinAbove(1.8*maxnoise)+7);

    if(setting::getsetname() == "78Ni"){
        for(int j=400; j<770; j++) derived.SetBinContent(j,0);
        for(int j=52; j<59; j++) derived.SetBinContent(j,0);
        if(t_min < 1250 || t_min > 1450) t_min = 1370;
    }

    double t_end = derived.GetBinCenter(derived.GetMinimumBin());

    if(setting::getsetname() == "78Ni"){
        if(t_end < 7700 || t_end > 8500) t_end = 8200;
    }

    TFile derivedfile(Form("output/MINOS/%s.root", setting::getsetname().c_str()),
                      "Update");
    derivedfile.Delete(Form("%s;1",i.substr(13,9).c_str()));
    derived.Write();

    return {t_min,t_end};
}