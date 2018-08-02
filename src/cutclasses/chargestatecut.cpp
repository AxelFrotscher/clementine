//
// Created by afrotscher on 8/1/18.
//

#include "cutclasses/chargestatecut.h"
#include "libconstant.h"

using namespace std;

thread ccsc::innerloop(treereader *tree, std::vector<std::atomic<bool>> &goodevents,
                       std::vector<int> range) {
    // precious tight inner loop
    // Cloning histograms
    vector<TH2D> _cschist;
    for(auto &i: cschist) _cschist.emplace_back(TH2D(i));

    // Step 2: Preparing variables
    const int downscale = (int)((range.at(1)-range.at(0))/100.);
    int threadno = range.at(0)/(range.at(1)-range.at(0));
    int i = range.at(0); // counting variable
    double brhoratio = 0;

    while(i<range.at(1)){
        if(goodevents.at(i)){
            tree->getevent(i);
            brhoratio = tree->BigRIPSBeam_brho[3]/tree->BigRIPSBeam_brho[2];
            _cschist.at(0).Fill(brhoratio, tree->BigRIPSBeam_brho[2]);

            if(mycut.at(threadno)->IsInside(brhoratio,tree->BigRIPSBeam_brho[2])){
                _cschist.at(1).Fill(brhoratio,tree->BigRIPSBeam_brho[2]);
            }
            else goodevents.at(i).exchange(false);
        }
        i++;

        if(!((i-range.at(0))%downscale)){
            consolemutex.lock();
            progressbar(i-range.at(0), range.at(1)-range.at(0), threadno);
            consolemutex.unlock();
        }
    }

    // Step 3: rejoin the data
    unitemutex.lock();
    for(uint i=0; i<_cschist.size();i++) cschist.at(i).Add(new TH2D(_cschist.at(i)));
    unitemutex.unlock();
}

void ccsc::analyse(std::vector<treereader*> tree, TFile* output){
    threads = (int)tree.size();
    printf("now beginning with the CCSC. %i threads.\n", threads);

    // Generate output histogram
    cschist.emplace_back(
            TH2D("csc", "Charged state change", 1000,0.8,1.2,1000,3,7));
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

    // Get Cut
    TFile cutfile("config/cut.root");
    if(!(cutfile.IsOpen())) __throw_invalid_argument("Could not open cut file "
                                                     "at config/rhocut.root");

    if(runinfo::transsize == goodevents.size())
        for(auto &i: tree) mycut.push_back((TCutG*)cutfile.Get("brhoempty"));
    else if(runinfo::emptysize == goodevents.size()){
        for(auto &i: tree) mycut.push_back((TCutG*)cutfile.Get("emptybrhocut"));
    }
    else for(auto &i: tree) mycut.push_back((TCutG*)cutfile.Get("brhocut"));
    if(!mycut.at(0)) __throw_invalid_argument("Could not load cut from file!\n");

    int cutpre = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<int> ranges = {i*goodevents.size()/threads,
                              (i+1)*goodevents.size()/threads-1};
        th.emplace_back(thread(&ccsc::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for (auto &i: th) i.join();

    int cutafter = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    printf("\nCCSC Cut out %i Events %f %%\n", cutpre-cutafter,
           100*(cutpre-cutafter)/(double)goodevents.size());

    output->mkdir("CSC");
    output->cd("CSC");
    for(auto hist: cschist) hist.Write();
    output->cd("");

}