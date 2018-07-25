#include "histograms.hh"
#include "MakeAllTree_78Ni.hh"
#include <iostream>
#include <algorithm>
#include <thread>
#include "TThread.h"
#include "libconstant.h"

R__LOAD_LIBRARY(libanacore.so)

using namespace std;

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
    ROOT::EnableThreadSafety();

    // The Scope of this project is to analyse the 2015 SEASTAR
    printf("Welcome to the analysis program for 2015 SEASTAR DATA\n"
           "Would you like to analyse a raw file first? \n"
           " (1) yes or (0) no: ");
    int analyse_raw = 0;
    if(!(cin >> analyse_raw)) throw invalid_argument("WTF");

    // Step 1: analyse the raw data
    const vector<uint> analysedfile{57}; // 57 Index of analysed file (first, offset)
    const vector<string> input = getlist("config/minosridf.txt");

    // Define good runs
    const vector<uint> goodruns{1,3,4,5,7,8,9,10,11,13,14,19,26,29,32,33,35,37,
                                       38,40,41,42,46,49,50,51,52,54,57,58,59};

    const vector<uint> transmissionrun{54};
    const vector<uint> emptyrun{52,53};

    vector<string> output;
    for(auto run : goodruns) {
        output.push_back(runinfo::prefix +
                         input.at(analysedfile.at(0)+run).substr(34, 9)
                         + ".root");
    }

    vector<string> transmissionout;
    for(auto run: transmissionrun)
        transmissionout.push_back(runinfo::prefix + input.at(run).substr(34,9)
                                  + ".root");

    vector<string> emptyout;
    for(auto run: emptyrun)
        emptyout.push_back(runinfo::prefix + input.at(run).substr(34,9) + ".root");

    switch(analyse_raw){
        case 1:{ // Analyse SEASTAR-DATA
            cout << "Which file to analyse?"
                    " [0] empty, [1] trans, [2] file " << endl;
            uint i = 0;
            if(!(cin >> i)) throw invalid_argument("WTH no int");

            switch(i){
                case 0:{
                    cout << "Which emptyrun? [0-1] " << endl;
                    if(!(cin >> i)) throw invalid_argument("WTH");
                    cout << "(Empty) Analyzing SEASTAR:"
                         << input.at(emptyrun.at(i)) << endl;
                    generatetree(input.at(emptyrun.at(i)),
                                 emptyout.at(i));
                    break;
                }
                case 1:{
                    cout << "(Trans) Analyzing SEASTAR:"
                         << input.at(transmissionrun.at(0)) << endl;
                    generatetree(input.at(transmissionrun.at(0)), transmissionout.at(0));
                    break;
                }
                case 2:{
                    cout << "Which physicsrun? [0-30] " << endl;
                    if(!(cin >> i)) throw invalid_argument("WTH");
                    cout << "(Empty) Analyzing SEASTAR:"
                         << input.at(analysedfile.at(0)+goodruns.at(i)) << endl;
                    generatetree(input.at(analysedfile.at(0)+goodruns.at(i)),
                                 output.at(i));
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