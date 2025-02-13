//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/targetcut.h"
#include "histogram_cuts.hh"
#include "progress.h"
#include <thread>
#include <numeric>

using std::vector, std::cout, std::endl, std::atomic, std::string, std::thread;

void targetcut::innerloop(treereader &tree, vector<int> range) {
    /// Step 1: Cloning histograms
    decltype(tarhist) _tarhist(tarhist);

    /// Step 2: preparing variables
    const uint threadno = range.at(0)/(range.at(1)-range.at(0));
    
    progressbar progress(range.at(1)-range.at(0), threadno);

    vector<vector<double>> slopex(2,vector<double>(0)); //0:x,y 1:val, z
    vector<vector<double>> slopey(2,vector<double>(0)); //0:x,y 1:val, z

    // Get temporary mean value:
    double meanx = 0, meany = 0, meanxz=0, meanyz =0, fXTar=0, fYTar=0;

    // Do test to have right PPAC's
    tree.getevent(range.at(0));
    for(int i=ppacoffset; i<(ppacoffset+f8ppac); i++){
        if(tree.BigRIPSPPAC_fpl[i] != 8){
            cout << "PPAC FPl from Index " << i << " not 8 but "
                 << tree.BigRIPSPPAC_fpl[i] << endl;
            std::__throw_invalid_argument("Wrong PPAC F8 indizes");
        }
    }

    for(int i =range.at(0); i < range.at(1); i++){
        if(goodevents.at(i).at(0)){
            // Determine cut only on good events
            tree.getevent(i);
            // 1. Fill valid events
            for(int j =ppacoffset; j<(ppacoffset+f8ppac); j++) {
                if (std::abs(tree.BigRIPSPPAC_fX[j] - magnum) > 1){
                    slopex.at(0).push_back(tree.BigRIPSPPAC_fX[j]);
                    slopex.at(1).push_back(f8z.at(j-ppacoffset).at(0));
                }
                if (std::abs(tree.BigRIPSPPAC_fY[j] - magnum) > 1){
                    slopey.at(0).push_back(tree.BigRIPSPPAC_fY[j]);
                    slopey.at(1).push_back(f8z.at(j-ppacoffset).at(1));
                }
            }
            if((slopex.at(0).size() < 2) || (slopey.at(0).size() <2)){
                for(auto &j:goodevents.at(i)) j.exchange(false);
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
                fXTar= meanx + (f8offset - meanxz)*slope(slopex.at(1),slopex.at(0));
                fYTar= meany + (f8offset - meanyz)*slope(slopey.at(1),slopey.at(0));

                _tarhist.at(0).Fill(fXTar,fYTar);

                if(pow(fXTar*fXTar+fYTar*fYTar,0.5) < targetradius){ // target cut
                    _tarhist.at(1).Fill(fXTar,fYTar);
                }
                else for(auto &j:goodevents.at(i)) j.exchange(false);

            }

            slopex.at(0).clear();
            slopey.at(0).clear();
            slopex.at(1).clear();
            slopey.at(1).clear();
        }
        progress.increaseevent();
    }

    // Step 3: rejoining histograms
    unitemutex.lock();
    for(ulong i=0; i<tarhist.size();i++)  tarhist.at(i).Add(&_tarhist.at(i));
    unitemutex.unlock();
    progressbar::reset();
}

void targetcut::analyse(const vector<string> &input, TFile *output) {

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Targetcut with %ix power!\n", threads);

    // Get relevant keys
    vector<string> keys{"BigRIPSPPAC.fX", "BigRIPSPPAC.fY", "BigRIPSPPAC.fpl"};
    for(auto &i: tree) i.setloopkeys(keys);

    tarhist = {{"target", "Beam Profile F8", 1000,-50,50,1000,-50,50},
               {"targetcut", "Cut Beam Profile F8", 1000,-50,50,1000,-50,50}};

    for(auto &histo: tarhist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("y [mm]");
        histo.GetXaxis()->SetTitle("x [mm]");
    }

    int cutpre =  0;
    for(auto &i:goodevents) cutpre += i.at(0);

    progressbar finishcondition;
    vector<thread> th;
    for(int i=0; i<threads; i++){
        vector<int> ranges = {(int)(i*goodevents.size()/threads),
                              (int)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&targetcut::innerloop, this, std::ref(tree.at(i)),
                               ranges));
    }

    for (auto &i: th) i.detach();

    while(progressbar::ongoing()) finishcondition.draw();

    int cutafter = 0;
    for(auto &i:goodevents) cutafter += i.at(0);

    printf("\nTarget Cut out %i Events %.3f %%\n", cutpre-cutafter,
           100*(cutpre-cutafter)/(double)goodevents.size());

    output->mkdir("Target");
    output->cd("Target");
    for(auto hist: tarhist) hist.Write();
    output->cd("");
}