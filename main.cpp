#include "histograms.hh"
#include "MakeAllTree_78Ni.hh"
#include <iostream>
#include <thread>

R__LOAD_LIBRARY(libanacore.so)

using namespace std;

vector<string> getlist(){
    // read in the file list for all .rdif's
    ifstream infile("config/minosrdif.txt");
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
    vector<string> input = getlist();
    string output = "build/output/out.root";

    if(analyse_raw == 1){
        generatetree(input.at(34), output);
    }
    else if (analyse_raw !=0) throw invalid_argument("Neither 0 nor 1");
    // Step 2: generate the histograms
    printf("\nNow proceeding to histogram creation ...\n");
    makehistograms(output);
    return 0;
}
