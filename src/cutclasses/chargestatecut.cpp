//
// Created by afrotscher on 8/1/18.
//

#include "cutclasses/chargestatecut.h"
#include "libconstant.h"
#include "zdssetting.h"
#include <thread>
#include <numeric>

using namespace std;

void ccsc::innerloop(treereader *tree, std::vector<std::atomic<bool>> &goodevents,
                       std::vector<uint> range) {
    // precious tight inner loop
    // Cloning histograms
    vector<TH2D> _cschist;
    for(auto &i: cschist) _cschist.emplace_back(TH2D(i));

    // Step 2: Preparing variables
    uint threadno = range.at(0)/(range.at(1)-range.at(0));
    double brhoratio = 0;

    progressbar progress(range.at(1)-range.at(0), threadno);

    for(int i=range.at(0); i<range.at(1); i++){
         if(goodevents.at(i)){
            tree->getevent(i);

            brhoratio = tree->BigRIPSBeam_brho[3]/tree->BigRIPSBeam_brho[2];
            _cschist.at(0).Fill(brhoratio, tree->BigRIPSBeam_brho[2]);

            if(mycut.at(threadno)->IsInside(brhoratio,tree->BigRIPSBeam_brho[2])){
                _cschist.at(1).Fill(brhoratio,tree->BigRIPSBeam_brho[2]);
            }
            else goodevents.at(i).exchange(false);
        }

        progress.increaseevent();
    }

    // Step 3: rejoin the data
    unitemutex.lock();
    for(uint i=0; i<_cschist.size();i++)
        cschist.at(i).Add(new TH2D(_cschist.at(i)));
    unitemutex.unlock();
    progress.reset();
}

void ccsc::analyse(const std::vector<std::string> input, TFile* output){
    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto &h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    printf("Beginning with the CCSC. %i threads.\n", threads);

    // Generate output histogram
    cschist.emplace_back(
            TH2D("csc", "Charged state change", 500,0.9,1.1,1000,3,7));
    cschist.emplace_back(
            TH2D("csccut", "Charged state change cut", 500,0.95,1.05,1000,3,7));

    for(auto &histo: cschist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("B#rho F8-9");
        histo.GetXaxis()->SetTitle("B#rho F8-9 / B#rho F9-11");
    }

    // Get relevant keys
    vector<string> keys{"BigRIPSBeam.brho"};
    for(auto &i: tree) i->setloopkeys(keys);

    // Get Cut from right setting | mode
    setting set;
    if(set.isemptyortrans()){ // Cut after F7 not sensible for empty/trans runs
        cout << "Empty or trans run. Not doing CCSC cut." << endl; return;
    }
    for(auto &i: tree) mycut.push_back(set.getbrhocut());

    int cutpre = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&ccsc::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for (auto &t: th) t.detach();

    // Setup synchronization class
    progressbar finishcondition;
    while(finishcondition.ongoing()) finishcondition.draw();

    int cutafter = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    printf("\nCCSC Cut out %i Events %f %%\n", cutpre-cutafter,
           100.*(cutpre-cutafter)/(double)goodevents.size());

    output->mkdir("CSC");
    output->cd("CSC");
    for(auto hist: cschist) hist.Write();
    output->cd("");
}