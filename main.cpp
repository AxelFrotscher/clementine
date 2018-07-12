#include "histograms.hh"
#include "MakeAllTree_78Ni.hh"
#include <iostream>
#include <algorithm>
#include <thread>
#include "TThread.h"

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
    const vector<int> goodruns{1,3,4,5,7,8,9,10,11,13,14,19,26,29,32,33,35,37,
                                       38,40,41,42,46,49,50,51,52,54,57,58,59};

    const vector<int> emptyrun{54}; //empty 52,53

    vector<string> output;
    for(auto run : goodruns) {
        output.push_back("/d/d02-1/ag_ob/SEASTAR2_DATA/root/" +
                         input.at(analysedfile.at(0)+run).substr(34, 9)
                         + ".root");
    }

    vector<string> emptyout;
    for(auto run: emptyrun)
        emptyout.push_back("/d/d02-1/ag_ob/SEASTAR2_DATA/root/" +
                           input.at(run).substr(34,9) + ".root");

    switch(analyse_raw){
        case 1:{ // Analyse SEASTAR-DATA
            cout << "Which file to analyse? " << endl;
            uint i = 0;
            if(!(cin >> i)) throw invalid_argument("WTH");
            cout << "Analyzing SEASTAR:"
                 << input.at(analysedfile.at(0)+goodruns.at(i)) << endl;
            generatetree(input.at(analysedfile.at(0)+goodruns.at(i)),
                         output.at(i));
            break;
        }
        case 0:{
            printf("Now proceeding to make histograms\n");
            makehistograms(output);
            //makehistograms(emptyout);
            break;
        }
        default:
            __throw_invalid_argument("Option not available");
    }
    return 0;
}