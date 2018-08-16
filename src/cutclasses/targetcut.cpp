//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/targetcut.h"
#include "histogram_cuts.hh"
#include "progress.h"
#include <thread>

using namespace std;

void targetcut::innerloop(treereader *tree, std::vector<std::atomic<bool>>
                            &goodevents, std::vector<uint> range) {
    //Step 1: Cloning histograms
    vector<TH2D> _tarhist;
    for(auto &i:tarhist) _tarhist.emplace_back(TH2D(i));

    //Step 2: preparing variables
    const auto downscale = (int)((range.at(1)-range.at(0))/100.);
    uint threadno = range.at(0)/(range.at(1)-range.at(0));

    progressbar progress(range.at(1)-range.at(0), threadno);

    int i = range.at(0); // counting variable

    vector<vector<double>> slopex(2,vector<double>(0)); //0:x,y 1:val, z
    vector<vector<double>> slopey(2,vector<double>(0)); //0:x,y 1:val, z

    // Get temporary mean value:
    double meanx = 0, meany = 0, meanxz=0, meanyz =0, fXTar=0, fYTar=0;

    while(i < range.at(1)){
        if(goodevents.at(i)){
            // Determine cut only on good events
            tree->getevent(i);
            // 1. Fill valid events
            for(uint i =ppacoffset; i<(ppacoffset+f8ppac); i++) {
                if (abs(tree->BigRIPSPPAC_fX[i] - magnum) > 1){
                    slopex.at(0).push_back(tree->BigRIPSPPAC_fX[i]);
                    slopex.at(1).push_back(f8z.at(i-ppacoffset).at(0));
                }
                if (abs(tree->BigRIPSPPAC_fY[i] - magnum) > 1){
                    slopey.at(0).push_back(tree->BigRIPSPPAC_fY[i]);
                    slopey.at(1).push_back(f8z.at(i-ppacoffset).at(1));
                }
            }
            if((slopex.at(0).size() < 2) || (slopey.at(0).size() <2)){
                goodevents.at(i).exchange(false);
            }
            else{
                meanx = accumulate(slopex.at(0).begin(), slopex.at(0).end(), 0.0)/
                        slopex.at(0).size();
                meany = accumulate(slopey.at(0).begin(), slopey.at(0).end(), 0.0)/
                        slopey.at(0).size();
                meanxz= accumulate(slopex.at(1).begin(), slopex.at(1).end(), 0.0)/
                        slopex.at(1).size();
                meanyz= accumulate(slopey.at(1).begin(), slopey.at(1).end(), 0.0)/
                        slopey.at(1).size();
                fXTar= meanx - (meanxz-f8offset)*slope(slopex.at(1),slopex.at(0));
                fYTar= meany - (meanyz-f8offset)*slope(slopey.at(1),slopey.at(0));

                _tarhist.at(0).Fill(fXTar,fYTar);

                if(pow(fXTar*fXTar+fYTar*fYTar,0.5) < targetradius){ // target cut
                    _tarhist.at(1).Fill(fXTar,fYTar);
                }
                else goodevents.at(i).exchange(false);

            }

            slopex.at(0).clear();
            slopey.at(0).clear();
            slopex.at(1).clear();
            slopey.at(1).clear();
        }

        i++;
        progress.increaseevent();
    }

    // Step 3: rejoining histograms
    unitemutex.lock();
    for(uint i=0; i<tarhist.size();i++){
        tarhist.at(i).Add(new TH2D(_tarhist.at(i)));
    }
    unitemutex.unlock();
    progress.reset();
}

void targetcut::analyse(const std::vector<std::string> input, TFile *output) {
    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto &h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    printf("Targetcut with %ix power!\n", threads);

    // Get relevant keys
    vector<string> keys{"BigRIPSPPAC.fX", "BigRIPSPPAC.fY"};
    for(auto &i: tree) i->setloopkeys(keys);

    tarhist.emplace_back(TH2D("target", "Beam Profile F8", 1000,-50,50,1000,-50,50));
    tarhist.emplace_back(TH2D("targetcut", "Cut Beam Profile F8", 1000,-50,50,1000,-50,50));

    for(auto &histo: tarhist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("x [mm]");
        histo.GetXaxis()->SetTitle("y [mm]");
    }

    int cutpre = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                              (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&targetcut::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for (auto &i: th) i.detach();

    progressbar finishcondition;
    while(finishcondition.ongoing()) finishcondition.draw();

    int cutafter = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    printf("\nTarget Cut out %i Events %f %%\n", cutpre-cutafter,
           100*(cutpre-cutafter)/(double)goodevents.size());

    output->mkdir("Target");
    output->cd("Target");
    for(auto hist: tarhist) hist.Write();
    output->cd("");

}