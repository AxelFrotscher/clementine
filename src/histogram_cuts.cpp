//
// Created by axel on 05.06.18.
//

#include "histogram_cuts.hh"
#include "MakeAllTree_78Ni.hh"
#include "libconstant.h"
#include "histograms.hh"
#include "cutclasses/triggercut.h"

using namespace std;

const bool closeness(const vector<double> &d, double sigma){
    // Standard deviation for vector d
    double sum = accumulate(d.begin(),d.end(),0.0);
    double mean = sum /d.size();
    vector<double> diff(d.size());

    transform(d.begin(),d.end(), diff.begin(),
              [mean](double x) { return x - mean; });

    double sq_sum= inner_product(diff.begin(),diff.end(),diff.begin(),0);
    return sigma * mean > sqrt(sq_sum / d.size());
}

double linfit(double *x, double *par){
    // Linear fit function
    return par[0] + par[1]*x[0];
}

/*void targetcut(treereader *tree, TFile *output,
               vector<atomic<bool>> &goodevents){
    // This method aims to reconstruct the target impinging position at F8
    printf("Now performing the Target cut... \n");

    vector<TH2D> tarhist{
        TH2D("target", "Beam Profile F8", 1000,-50,50,1000,-50,50),
        TH2D("targetcut", "Cut Beam Profile F8", 1000,-50,50,1000,-50,50)};

    for(auto &histo: tarhist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("x [mm]");
        histo.GetXaxis()->SetTitle("y [mm]");
    }

    // Get relevant keys
    vector<string> keys{"BigRIPSPPAC.fX", "BigRIPSPPAC.fY"};
    tree->setloopkeys(keys);

    const int ppacoffset = 18; // [18] is F8PPAC-1A
    const int f8ppac = 4;
    const int  magnum = -9999;
    const int f8offset = -2; // mm

    // F8 PPAC positions relative to F8 0:PPAcno 1: z(x),z(y)
    const vector<vector<double>> f8z{{-1310.7,-1302.1},{-1273.3,-1281.9},
                                     {-810.7,-802.1},{-773.3,-781.9}};

    vector<vector<double>> slopex(2,vector<double>(0)); //0:x,y 1:val, z
    vector<vector<double>> slopey(2,vector<double>(0)); //0:x,y 1:val, z

    // Progress Bar setup
    uint eventno=0; // counting variable
    uint totevents = goodevents.size();
    const int downscale = totevents/100; // every percent
    // Get Number of Good events before
    int  cutcount =0;
    // Get temporary mean value:
    double meanx = 0, meany = 0, meanxz=0, meanyz =0, fXTar=0, fYTar=0;

    while(eventno < totevents){
        if(goodevents.at(eventno)){
            // Determine cut only on good events
            tree->getevent(eventno);
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
                // cutting
                cutcount++;
                if(goodevents.at(eventno)) goodevents.at(eventno).exchange(false);
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
                tarhist.at(0).Fill(fXTar,fYTar);
                if(pow(fXTar*fXTar+fYTar*fYTar,0.5) < 20){ // target cut
                    tarhist.at(1).Fill(fXTar,fYTar);
                }
                else{
                    cutcount++;
                    goodevents.at(eventno).exchange(false);
                }
            }

            slopex.at(0).clear();
            slopey.at(0).clear();
            slopex.at(1).clear();
            slopey.at(1).clear();
        }

        eventno++;
        if(!(eventno%downscale)){
            consolemutex.lock();
            progressbar(eventno, totevents,5);
            consolemutex.unlock();
        }
    }

    writemutex.lock();
    printf("\nTarget Cut out %i Events %f %%\n", cutcount,
           100*cutcount/(double)totevents);

    output->mkdir("Target");
    output->cd("Target");
    for(auto hist: tarhist) hist.Write();
    output->cd("");
    writemutex.unlock();
}*/