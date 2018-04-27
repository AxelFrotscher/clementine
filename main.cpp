#include "histograms.hh"
#include "MakeAllTree_78Ni.hh"
#include "DALIcalibration.hh"
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
    string output = "build/output/out.root";
    vector<string> dalioutput;
    const int analysedfile = 1; // Index of analysed file
    switch(analyse_raw){
        case 2:{
            vector<string> input = getlist("config/daliridf.txt");
            cout << "Analyzing DALI: " << input.at(analysedfile) << endl;
            //Generating output matching to input
            for(int i=0; i<input.size(); i++){
                string number = input.at(i).substr(21,21); //Start, length
                dalioutput.push_back("build/output/dali/" + number + ".root");
                generatetree(input.at(i), dalioutput.back(), true);
            }
            makehistograms(output, true);
            break;
        }
        case 1:{ // Analyse SEASTAR-DATA
            vector<string> input = getlist("config/minosridf.txt");
            cout << "Analyzing SEASTAR:" << input.at(34) << endl;
            generatetree(input.at(34), output, false);
        }
        case 0:{
            printf("Now proceeding to make histograms");
            makehistograms(output, false);
            break;
        }
        default:
            __throw_invalid_argument("Option not available");
    }
    return 0;
}
