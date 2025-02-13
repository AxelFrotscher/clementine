//
// Created by afrotscher on 7/30/18.
//

#include "cutclasses/ppaccut.h"
#include "progress.h"
#include <thread>
#include <numeric>
#include "zdssetting.h"

using std::vector, std::string, std::thread, std::atomic;

void ppaccut::innerloop(treereader &tree, vector<int> range) {
    /// Step 1: cloning histograms
    decltype(effPPAC) _effPPAC(effPPAC);
    decltype(sumdiffppac) _sumdiffppac(sumdiffppac);

    /// Step 2: Preparing variables
    const uint threadno = range.at(0)/(range.at(1)-range.at(0));
    vector<bool> temptruth(4, true); // variable for position checking

    progressbar progress(range.at(1)-range.at(0), threadno);

    /// Step 3: glorious loop
    for(int i=range.at(0); i<range.at(1); i++){
        progress.increaseevent();
        if(!goodevents.at(i).at(0)) continue;   // analyse good events only
        tree.getevent(i);

        // Get Efficiency relative to the last Plastic
        if (tree.BigRIPSPlastic_fTime[pl11position] >0){
            for(int j=1; j<=numplane; j++){
                if(tree.BigRIPSPPAC_fFiredX[j-1])
                    _effPPAC.at(0).SetBinContent(j, _effPPAC.at(0).GetBinContent(j)+1);
                if(tree.BigRIPSPPAC_fFiredY[j-1])
                    _effPPAC.at(1).SetBinContent(j, _effPPAC.at(1).GetBinContent(j)+1);
            }
            efficiencytotal++;
        }

        // Fill Sum and differences of time signals to check
        for(int j=0; j<numplane;j++){
            _sumdiffppac.at(0).at(0).Fill(j,tree.BigRIPSPPAC_fTSumX[j]);
            _sumdiffppac.at(0).at(1).Fill(j,tree.BigRIPSPPAC_fTDiffX[j]);
            _sumdiffppac.at(1).at(0).Fill(j,tree.BigRIPSPPAC_fTSumY[j]);
            _sumdiffppac.at(1).at(1).Fill(j,tree.BigRIPSPPAC_fTDiffY[j]);
        }

        // Make Cut conditions: each focal Plane needs to have one valid entry
        // := Sum for x and y is in range
        bool temp = true;

        // Make empty and trans runs only go over points to F7
        for(unsigned long k =0; k<(ppacplane.size()); k++){
            for(unsigned int j=0; j<4; j++){ // Loop over all 4 PPAC's per Focal Plane
                int cpl = ppacplane.at(k).at(j);
                temptruth.at(j) =
                    (tree.BigRIPSPPAC_fTSumX[cpl] > ppacrange.at(k).at(j).at(0) &&
                     tree.BigRIPSPPAC_fTSumX[cpl] < ppacrange.at(k).at(j).at(1) &&
                     tree.BigRIPSPPAC_fTSumY[cpl] > ppacrange.at(k).at(j).at(2) &&
                     tree.BigRIPSPPAC_fTSumY[cpl] < ppacrange.at(k).at(j).at(3));
                if(temptruth.at(j)) break; // We got our nice event at this F-point
            }
            // Recursively check each focal plane for at least 1 good Signal
            temp *= accumulate(temptruth.begin(), temptruth.end(),0);

            if(k==2 && !temp){ // if fault before F8 then set both to false
                for(auto &l:goodevents.at(i)) l.exchange(false);
                break;
            }
        }
        if(!temp) goodevents.at(i).at(1).exchange(false);
    }

    // Step 4 reuniting the diagrams
    unitemutex.lock();
    for(uint i=0;i<effPPAC.size();i++){
        effPPAC.at(i).Add(&_effPPAC.at(i));
    }
    for(uint i=0; i<sumdiffppac.size();i++){
        for(uint j=0; j<sumdiffppac.at(0).size(); j++){
            sumdiffppac.at(i).at(j).Add(&_sumdiffppac.at(i).at(j));
        }
    }
    unitemutex.unlock();

    progressbar::reset();
}

void ppaccut::analyse(const std::vector<std::string> &input, TFile* output){

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Reconstruction of PPAC's. %i threads.\n", threads);

    // Generate List of all Keys that are to be used
    vector<string> keys{"BigRIPSPlastic.fTime", "BigRIPSPPAC.fFiredX",
                        "BigRIPSPPAC.fFiredY", "BigRIPSPPAC.fTSumX",
                        "BigRIPSPPAC.fTDiffX", "BigRIPSPPAC.fTSumY",
                        "BigRIPSPPAC.fTDiffY", "BigRIPSPPAC.name"};

    for(auto &i:tree) i.setloopkeys(keys);

    //36 Values per Array (Event)
    effPPAC = {{"effPPACX","Efficiency of PPAC X",numplane,0,(double)numplane},
               {"effPPACY","Efficiency of PPAC Y",numplane,0,(double)numplane}};

    std::array<string, 4> arrname = { "PPACXsum","PPACXdiff","PPACYsum",
                                 "PPACYdiff"};
    std::array<string, 4> arrtitle = {
            "Sum of Signals PPACX", "Difference of Signals PPACX",
            "Sum of Signals PPACY", "Difference of Signals PPACY"};

    sumdiffppac = {{{arrname.at(0).c_str(),arrtitle.at(0).c_str(),
                        numplane,0,(double)numplane,300,000,300},
                    {arrname.at(1).c_str(),arrtitle.at(1).c_str(),
                        numplane,0,(double)numplane,800,-200,200}},
                   {{arrname.at(2).c_str(),arrtitle.at(2).c_str(),
                        numplane,0,(double)numplane,150,0,150},
                    {arrname.at(3).c_str(),arrtitle.at(3).c_str(),
                        numplane,0,(double)numplane,800,-200,200}}};

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

    std::array<int, 2> cutpre = {0,0};
    for(auto &i: goodevents){
        cutpre.at(0) += i.at(0);
        cutpre.at(1) += i.at(1);
    }

    progressbar finishcondition;
    vector<thread> th;
    for(int i=0; i<threads; i++){
        vector<int> ranges = {(int)(i*goodevents.size()/threads),
                               (int)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&ppaccut::innerloop, this, std::ref(tree.at(i)),
                               ranges));
    }

    for (auto &i: th) i.detach();

    while(progressbar::ongoing()) finishcondition.draw();

    output->mkdir("PPAC/FiredEff");
    output->cd("PPAC/FiredEff");
    for(auto &i : effPPAC){
        i.Scale(1/(double)efficiencytotal);
        i.Write();
    }

    const std::array<string, 2> xy = {"X","Y"};
    const std::array<string, 2> sd = {"Sum", "Diff"};

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

    std::array<int, 2> cutafter = {0,0};
    for(auto &i:goodevents){
        cutafter.at(0) += i.at(0);
        cutafter.at(1) += i.at(1);
    }

    printf("\nPPAC Cut out %i F1-7 Events (%.2f %%) %i F1-11 Events (%.2f %%) \n",
           cutpre.at(0)-cutafter.at(0), 100*(cutpre.at(0)-cutafter.at(0))/
                                        (double)goodevents.size(),
           cutpre.at(1)-cutafter.at(1), 100*(cutpre.at(1)-cutafter.at(1))/
                                        (double)goodevents.size());
}