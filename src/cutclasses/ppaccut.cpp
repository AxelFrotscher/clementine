//
// Created by afrotscher on 7/30/18.
//

#include "cutclasses/ppaccut.h"
#include "progress.h"
#include <thread>

using namespace std;

void ppaccut::innerloop(treereader *tree,
                             std::vector<std::atomic<bool>> &goodevents,
                             std::vector<uint> range) {
    // Step 1: cloning histograms
    vector<TH1D> _effPPAC;
    vector<vector<TH2D>> _sumdiffppac;
    for(auto &i : effPPAC) _effPPAC.emplace_back(TH1D(i));

    for(auto &i: sumdiffppac){
        vector<TH2D> _sdppac;
        for (auto &j: i) _sdppac.emplace_back(TH2D(j));
        _sumdiffppac.emplace_back(_sdppac);
    }

    // Step 2: Preparing variables
    uint threadno = range.at(0)/(range.at(1)-range.at(0));
    uint i = range.at(0); // counting variable
    vector<bool> temptruth(4, true); // variable for position checking

    progressbar progress(range.at(1)-range.at(0), threadno);

    // Step 3: glorious loop
    while(i<range.at(1)){
        if(goodevents.at(i)){  // analyse good events only
            tree->getevent(i);

            // Get Efficiency relative to the last Plastic
            if (tree->BigRIPSPlastic_fTime[pl11position] >0){
                for(uint i=1; i<=numplane; i++){
                    if(tree->BigRIPSPPAC_fFiredX[i-1])
                        _effPPAC.at(0).SetBinContent(i, _effPPAC.at(0).GetBinContent(i)+1);
                    if(tree->BigRIPSPPAC_fFiredY[i-1])
                        _effPPAC.at(1).SetBinContent(i, _effPPAC.at(1).GetBinContent(i)+1);
                }
                efficiencytotal++;
            }

            // Fill Sum and differences of time signals to check
            for(uint i=0; i<numplane;i++){
                _sumdiffppac.at(0).at(0).Fill(i,tree->BigRIPSPPAC_fTSumX[i]);
                _sumdiffppac.at(0).at(1).Fill(i,tree->BigRIPSPPAC_fTDiffX[i]);
                _sumdiffppac.at(1).at(0).Fill(i,tree->BigRIPSPPAC_fTSumY[i]);
                _sumdiffppac.at(1).at(1).Fill(i,tree->BigRIPSPPAC_fTDiffY[i]);
            }

            // Make Cut conditions: each focal Plane needs to have one valid entry
            // := Sum for x and y is in range
            bool temp = true;

            for(uint k =0; k<ppacplane.size();k++){
                for(uint j=0; j<4; j++){ // Loop over all 4 PPAC's per Focal Plane
                    int cpl = ppacplane.at(k).at(j);
                    temptruth.at(j) =
                            (tree->BigRIPSPPAC_fTSumX[cpl] > ppacrange.at(k).at(j).at(0) &&
                             tree->BigRIPSPPAC_fTSumX[cpl] < ppacrange.at(k).at(j).at(1) &&
                             tree->BigRIPSPPAC_fTSumY[cpl] > ppacrange.at(k).at(j).at(2) &&
                             tree->BigRIPSPPAC_fTSumY[cpl] < ppacrange.at(k).at(j).at(3));
                    if(temptruth.at(j)) break; // We got our nice event at this F-point
                }
                // Recursively check each focal plane for at least 1 good Signal
                temp = temp*accumulate(temptruth.begin(), temptruth.end(),0);
            }
            if(!temp) goodevents.at(i).exchange(false);
        } // end of physics loop
        i++;
        progress.increaseevent();
    }

    // Step 4 reuniting the diagrams
    unitemutex.lock();
    for(uint i=0;i<effPPAC.size();i++){
        effPPAC.at(i).Add(new TH1D(_effPPAC.at(i)));
    }
    for(uint i=0; i<sumdiffppac.size();i++){
        for(uint j=0; j<sumdiffppac.at(0).size(); j++){
            sumdiffppac.at(i).at(j).Add(new TH2D(_sumdiffppac.at(i).at(j)));
        }
    }
    unitemutex.unlock();

    progress.reset();
}

void ppaccut::analyse(const std::vector<std::string> input, TFile* output){
    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto &h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    printf("now beginning with the reconstruction of the ppac's. %i threads.\n", threads);

    // Generate List of all Keys that are to be used
    vector<string> keys{"BigRIPSPlastic.fTime", "BigRIPSPPAC.fFiredX",
                        "BigRIPSPPAC.fFiredY", "BigRIPSPPAC.fTSumX",
                        "BigRIPSPPAC.fTDiffX", "BigRIPSPPAC.fTSumY",
                        "BigRIPSPPAC.fTDiffY", "BigRIPSPPAC.name"};

    for(auto &i:tree) i->setloopkeys(keys);

    //36 Values per Array (Event)
    effPPAC.emplace_back(TH1D("effPPACX", "Efficiency of PPAC X", numplane, 0,
                              numplane));
    effPPAC.emplace_back(TH1D("effPPACY", "Efficiency of PPAC Y", numplane, 0,
                              numplane));

    vector<string> arrname = { "PPACXsum","PPACXdiff","PPACYsum", "PPACYdiff"};
    vector<string> arrtitle = {
            "Sum of Signals PPACX", "Difference of Signals PPACX",
            "Sum of Signals PPACY", "Difference of Signals PPACY"};

    for(uint i = 0; i<2;i++) sumdiffppac.emplace_back(vector<TH2D>{});
    sumdiffppac.at(0).emplace_back(TH2D(arrname.at(0).c_str(),arrtitle.at(0).c_str(),
                                        numplane,0,numplane,300,000,300));
    sumdiffppac.at(0).emplace_back(TH2D(arrname.at(1).c_str(),arrtitle.at(1).c_str(),
                                        numplane,0,numplane,800,-200,200));
    sumdiffppac.at(1).emplace_back(TH2D(arrname.at(2).c_str(),arrtitle.at(2).c_str(),
                                        numplane,0,numplane,150,0,150));
    sumdiffppac.at(1).emplace_back(TH2D(arrname.at(3).c_str(),arrtitle.at(3).c_str(),
                                        numplane,0,numplane,800,-200,200));

    for(auto &hist : sumdiffppac){
        uint toggle =0;
        for(auto &h: hist){
            h.SetOption("colz");
            h.GetXaxis()->SetTitle("PPAC [ch]");
            if(!toggle) h.GetYaxis()->SetTitle("Sum [ns]");
            else h.GetYaxis()->SetTitle("diff [ns]");
            toggle++;
        }
    }

    for(auto &i: effPPAC) i.GetXaxis()->SetTitle("PPAC [ch]");
    effPPAC.at(0).GetYaxis()->SetTitle("PPACX(x)/Plastic(F11)");
    effPPAC.at(1).GetYaxis()->SetTitle("PPACY(x)/Pplastic(F11)");

    int cutpre = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&ppaccut::innerloop, this, tree.at(i), ref(goodevents),ranges));
    }

    for (auto &i: th) i.detach();

    progressbar finishcondition;
    while(finishcondition.ongoing()) finishcondition.draw();

    output->mkdir("PPAC/FiredEff");
    output->cd("PPAC/FiredEff");
    for(auto &i : effPPAC){
        i.Scale(1/(double)efficiencytotal);
        i.Write();
    }

    vector<string> xy = {"X","Y"};
    vector<string> sd = {"Sum", "Diff"};

    // Generating output structure
    for(uint j=0; j<2; j++){ // Loop over X and Y//
        string folder = "PPAC/" + sd.at(j);
        output->mkdir(folder.c_str());
        output->cd(folder.c_str());
        for(uint k=0; k<2; k++){     // Loop Sum and diff
            sumdiffppac.at(k).at(j).Write();
        }
    }

    output->cd("");

    int cutafter = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);
    printf("\nPPAC Cut out %i Events %f %%\n", cutpre-cutafter,
           100*(cutpre-cutafter)/(double)goodevents.size());

}