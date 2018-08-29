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

using namespace std;

calibpar p1;
mutex consolemutex;
mutex writemutex;

void dalicalib(treereader *tree, TFile *output){
    // This Method aims to calibrate the 187 detectors of DALI
    printf("Now beginning the Calibration of the NaI crystals... \n");

    vector<string> keys{"DALINaI", "DALINaI.fADC", "DALINaI.id"};
    tree->setloopkeys(keys);

    TH2D gammadetectors("dalispectra", "Spectrum of each Gamma Detector",
                        186,0,186,4096,0,4096);
    gammadetectors.SetOption("colz");
    gammadetectors.GetXaxis()->SetTitle("Detector Number");
    gammadetectors.GetYaxis()->SetTitle("ADC Channel");

    int numdet =0;

    // Progress Bar setup
    int currevt=0; // counting variable
    uint totevents = tree->NumEntries();

    while(tree->singleloop()){
        numdet = tree->DALINaI_;
        for(int i=0; i<numdet;i++){
            if(tree->DALINaI_fADC[i]) gammadetectors.Fill(tree->DALINaI_id[i],
                tree->DALINaI_fADC[i]);
        }
        currevt++;
    }

    output->mkdir("DALI");
    output->cd("DALI");
    gammadetectors.Write();
    output->cd("");
    printf("\nFinished DALI Calibration.\n");
}

 void makehistograms(const vector<string> input) {
    // Generating an outputfile that matches names with the input file
    string output = "build/output/" + input.at(0).substr(34,9) + "hist" +
                    to_string(input.size()) + ".root";
    string txtout = "build/output/" + input.at(0).substr(34,9) + "hist" +
                    to_string(input.size()) + ".txt";

    // Initialize ROOT and txt outputfile
    auto outputfile = new TFile(output.c_str(), "RECREATE");
    if(!outputfile->IsOpen()) __throw_invalid_argument("Output file not valid");
    txtwriter writetotxt(txtout); // Writer class

    TChain chain("tree");
    for(auto &i:input) chain.Add(i.c_str());

    cout << "Beginning reconstruction of " << chain.GetEntries()
         << " Elements." << endl;
    // Store events that cannot be used
    vector<atomic<bool>> goodevents((uint)chain.GetEntries());
    for(auto &i:goodevents) i.exchange(true);

    // Determine run type:
    if(runinfo::transsize == chain.GetEntries()){
        printf("!!! Analysing an transmission run !!!\n");
    }
    if(runinfo::emptysize == chain.GetEntries()){
        printf("!!! Analysing an empty target run !!!\n");
    }

    //triggercut(input, goodevents);
    ccsc(input,goodevents,outputfile);
    targetcut(input,goodevents,outputfile);
    plasticcut(input, goodevents, outputfile);
    ppaccut(input, goodevents, outputfile);
    iccut(input,goodevents,outputfile);
    higherorder(input, goodevents, outputfile);

    // Get Z vs. A/Q
    //PID(input,goodevents,outputfile, "111NbPPN");
    //PID(input,goodevents,outputfile, "111NbPP2N");
    PID(input,goodevents,outputfile, "111NbP2P");
    //PID(input, goodevents, outputfile, "110NbPPN");
    PID(input, goodevents, outputfile, "110NbP2P");
    PID(input, goodevents, outputfile, "110MoP3P");
    PID(input, goodevents, outputfile, "111MoP3P");
    PID(input, goodevents, outputfile, "112MoP3P");
    PID(input, goodevents, outputfile, "113MoP3P");
    PID(input, goodevents, outputfile, "112TcP3P");
    PID(input, goodevents, outputfile, "113TcP3P");
    PID(input, goodevents, outputfile, "114TcP3P");

    //makepid(input, outputfile, goodevents);
    printf("Made PID histograms in %s\n", output.c_str());
    //Get ADC Spectra for DALI
    //dalicalib(alt4dtree, outputfile);
    writetotxt.writetofile();

    cout << "Run has " <<100.* accumulate(goodevents.begin(),goodevents.end(),0)
                          /goodevents.size() << " % valid Elements" << endl;

    outputfile->Close();
}