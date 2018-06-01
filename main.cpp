#include "histograms.hh"
#include "MakeAllTree_78Ni.hh"
#include <iostream>
#include <thread>

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
    // The Scope of this project is to analyse the 2015 SEASTAR
    printf("Welcome to the analysis program for 2015 SEASTAR DATA\n"
           "Would you like to analyse the raw files first? \n"
           " (1) yes or (0) no: ");
    int analyse_raw = 0;
    if(!(cin >> analyse_raw)) throw invalid_argument("WTF");

    // Step 1: analyse the raw data
    const int analysedfile = 64; // Index of analysed file
    const vector<string> input = getlist("config/minosridf.txt");
    string output = "build/output/" + input.at(analysedfile).substr(16,9)+".root";

    switch(analyse_raw){
        case 1:{ // Analyse SEASTAR-DATA
            cout << "Analyzing SEASTAR:" << input.at(analysedfile) << endl;
            generatetree(input.at(analysedfile), output);
        }
        case 0:{
            printf("Now proceeding to make histograms\n");
            makehistograms(output);
            break;
        }
        default:
            __throw_invalid_argument("Option not available");
    }
    return 0;
}
