#include "histograms.hh"
#include "libconstant.h"
#include "chargestatecut.h"
#include "ICcut.h"
#include "plasticcut.h"
#include "ppaccut.h"
#include "targetcut.h"
#include "triggercut.h"
#include "higherorder.h"
#include "PID/pid.h"
#include "txtwriter.h"
#include <numeric>
#include "zdssetting.h"
#include "TGraphErrors.h"

using std::cout, std::endl, std::string, std::vector, std::atomic, std::to_string,
      std::__throw_invalid_argument;

calibpar p1;

 void makehistograms(const vector<string> input) {
    // Generating an outputfile that matches names with the input file
    // Determine run type:
    TChain chain("tree");
    for(auto &i:input) chain.Add(i.c_str());

    cout << "Beginning reconstruction of " << chain.GetEntries()
         << " Elements." << endl;
    // Store events that cannot be used
    //vector<vector<atomic<bool>>> goodevents((uint)chain.GetEntries(), vector<atomic<bool>>(2));
    //for(auto &i:goodevents) for(auto &j:i) j.exchange(true);

    vector<vector<atomic<bool>>> goodevents;
    for(uint i=0; i<(uint)chain.GetEntries(); i++){
        goodevents.emplace_back(vector<atomic<bool>>(2));
    }
    for(auto &i:goodevents) for(auto &j:i) j.exchange(true);

    setting set;
    set.setcountno((int)chain.GetEntries());
    const string settingname = set.getsetname();
    const string modename = set.getmodename();
    const bool minosbool = set.getminos();

     // lambda to generate expression
    auto gentxt = [input, settingname, modename, minosbool](auto suffix){
        if(minosbool){
            return "output/"+ settingname + "_" + modename + "_MINOS_" +
                   input.at(0).substr(34,9) + suffix;
        }
        else return "output/"+ settingname + "_" + modename + "_CS_" +
               input.at(0).substr(34,9) + suffix;
    };

    // Initialize ROOT and txt outputfile
    auto outputfile = new TFile(gentxt(".root").c_str(), "RECREATE");
    //if(!outputfile->IsOpen()) __throw_invalid_argument("Output file not valid");

    txtwriter writetotxt(gentxt(".txt")); // Writer class

     { // Make everything go out of scope to prevent memory overflow
         cout << "Making cuts..." << endl;
         triggercut(input, goodevents);
         ccsc(input, goodevents, outputfile);
         targetcut(input, goodevents, outputfile);
         ppaccut(input, goodevents, outputfile);
         plasticcut(input, goodevents, outputfile);
         iccut(input, goodevents, outputfile);
         higherorder(input, goodevents, outputfile);
     };

    // Get Z vs. A/Q
    const vector<string> reactionmodes = set.getreactions();

    // All reaction modes:
    TGraphErrors crosssection;
    crosssection.SetName(Form("cs%s", set.getsetname().c_str()));
    crosssection.SetTitle("Cross Sections Frotscher 2019");

    for(auto &i: reactionmodes){
        PID(input,goodevents,outputfile,i, crosssection);
    }

    TGraphErrors nancytcs = nancycs(set.getsetnumber());
    outputfile->cd();
    nancytcs.Write();
    crosssection.Write();

    printf("Made PID histograms in %s\n", gentxt(".root").c_str());

    if(!minosbool) writetotxt.writetofile();

    outputfile->Close();
}

TGraphErrors nancycs(const int &setnumber){
     // Create a TGraph with Nancies values
     TGraphErrors temp;
     temp.SetTitle("Cross Sections Hupin 2019");
     temp.SetName("csnancy");

     vector<vector<double>> massnumbers{
         {110,110,111,111,112,113},  //nb,mo,nb,mo,mo,tc
         {},
         {95,95,96,97},  //br,kr,kr,rb
         {100},
         {67,68,68,69,70}, // fe,fe,co,co,ni
         {72,73,74,74,75}, // ni,ni,ni,cu,cu
         {78,79,81}        // zn,zn,ga
     };

     vector<vector<double>> crosssections{
         {3.0,8,4.3,5.9,7.4,6.5},
         {},
         {2.5,6,6.7,4.7},
         {9},
         {7,9.2,5.4,8.0,12},
         {10,7.0,7.3,4.8,6.6},
         {8,5.3,4.7}
     };

     vector<vector<double>> crosssectione{
         {0.4,1,0.9,0.7,0.8,0.7},
         {},
         {0.7,3,0.7,0.6},
         {1},
         {0.3,0.6,0.5,0.5,1},
         {1,0.4,0.4,0.3,0.5},
         {1,0.4,0.3}
     };

     for(int i=0; i<massnumbers.at(setnumber).size(); i++){
         temp.SetPoint(temp.GetN(),massnumbers.at(setnumber).at(i),
                                   crosssections.at(setnumber).at(i));
         temp.SetPointError(temp.GetN()-1, 0, crosssectione.at(setnumber).at(i));
     }
     temp.SetPoint(temp.GetN(), 50,0);
     temp.SetPoint(temp.GetN(), 120,0);
     return temp;
 }