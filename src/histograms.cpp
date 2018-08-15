#include "MakeAllTree_78Ni.hh"
#include "histograms.hh"
#include "histogram_cuts.hh"
#include "libconstant.h"
#include "cutclasses/chargestatecut.h"
#include "cutclasses/ICcut.h"
#include "cutclasses/plasticcut.h"
#include "cutclasses/ppaccut.h"
#include "cutclasses/targetcut.h"
#include "cutclasses/triggercut.h"
#include "corrections/higherorder.h"
#include "PID/pid.h"

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
    Long64_t totevents = tree->NumEntries();
    const int downscale = 500; // every n-th event

    while(tree->singleloop()){
        numdet = tree->DALINaI_;
        for(int i=0; i<numdet;i++){
            if(tree->DALINaI_fADC[i]) gammadetectors.Fill(tree->DALINaI_id[i],
                tree->DALINaI_fADC[i]);
        }
        currevt++;
        if(!(currevt%downscale)) progressbar(currevt,totevents,0);
    }

    output->mkdir("DALI");
    output->cd("DALI");
    gammadetectors.Write();
    output->cd("");
    printf("\nFinished DALI Calibration.\n");
}

 void makehistograms(const vector<string> input) {

    TChain chain("tree");
    for(auto &i:input) chain.Add(i.c_str());

    // Generating an outputfile that matches names with the input file
    string output = "build/output/" + input.at(0).substr(34,9) + "hist" +
                    to_string(input.size()) + ".root";

    auto outputfile = new TFile(output.c_str(), "RECREATE");
    if(!outputfile->IsOpen()) __throw_invalid_argument("Output file not valid");

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

    triggercut(input, goodevents);
    ccsc(input,goodevents,outputfile);
    targetcut(input,goodevents,outputfile);
    plasticcut(input, goodevents, outputfile);
    ppaccut(input, goodevents, outputfile);
    iccut(input,goodevents,outputfile);
    higherorder(input, goodevents, outputfile);

    // Get Z vs. A/Q
    PID(input,goodevents,outputfile, "111NbPPN");
    PID(input,goodevents,outputfile, "111NbPP2N");
    PID(input,goodevents,outputfile, "111NbP2P");
    PID(input, goodevents, outputfile, "110NbPPN");
    PID(input, goodevents, outputfile, "110NbP2P");

    //makepid(input, outputfile, goodevents);
    printf("Made PID histograms in %s\n", output.c_str());
    //Get ADC Spectra for DALI
    //dalicalib(alt4dtree, outputfile);

    cout << "Run has " <<100.* accumulate(goodevents.begin(),goodevents.end(),0)
                          /goodevents.size() << " % valid Elements" << endl;

    outputfile->Close();
}