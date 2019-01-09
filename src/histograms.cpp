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

    auto gentxt = [input, settingname, modename](auto suffix){ // lambda to generate expression
        return "output/"+ settingname + "_" + modename + "_" +
               input.at(0).substr(34,9) + "hist" + to_string(input.size()) + suffix;
    };

    // Initialize ROOT and txt outputfile
    auto outputfile = new TFile(gentxt(".root").c_str(), "RECREATE");
    if(!outputfile->IsOpen()) __throw_invalid_argument("Output file not valid");
    txtwriter writetotxt(gentxt(".txt")); // Writer class

     { // Make everything go out of scope to prevent memory overflow
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

    // 111NbPPN, 111NbPP2N, 110NbPPN
    for(auto &i: reactionmodes){
        PID(input,goodevents,outputfile,i);

    }

    printf("Made PID histograms in %s\n", gentxt(".root").c_str());
    //Get ADC Spectra for DALI
    //dalicalib(alt4dtree, outputfile);
    writetotxt.writetofile();

    outputfile->Close();
}