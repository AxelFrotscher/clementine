//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/ICcut.h"

using namespace std;

thread iccut::innerloop(treereader *tree, std::vector<std::atomic<bool>>
                        &goodevents, std::vector<int> range) {
    //Step 1: Cloning histograms
    vector<vector<TH2D>> _comparediag;
    for(auto &i: comparediag){
        vector<TH2D> temp1d;
        for(auto &j :i) temp1d.emplace_back(TH2D(j));
        _comparediag.emplace_back(temp1d);
    }

    //Step 2: Preparing Variables
    const int downscale = (int)((range.at(1)-range.at(0))/100.);
    int threadno = range.at(0)/(range.at(1)-range.at(0));
    int i = range.at(0); // counting variable


    while(i < range.at(1)){
        if(goodevents.at(i)){
            // Determine cut only on good events
            tree->getevent(i);

            if((tree->BigRIPSIC_nhitchannel[0]*
                tree->BigRIPSIC_nhitchannel[1]) <16){//Require 4 hits per IC
                goodevents.at(i).exchange(false);
            }
            if((tree->BigRIPSIC_fADC[0][0] <0) ||
               (tree->BigRIPSIC_fADC[1][0] <0)) continue;
            for(int i=0; i<(numchannel-1);i++){
                if(tree->BigRIPSIC_fADC[0][i] >0)
                    _comparediag.at(i).at(0).Fill(
                            tree->BigRIPSIC_fADC[0][0],
                            tree->BigRIPSIC_fADC[0][i+1]- tree->BigRIPSIC_fADC[0][0]);
                if(tree->BigRIPSIC_fADC[1][i] >0)
                    _comparediag.at(i).at(1).Fill(
                            tree->BigRIPSIC_fADC[1][0],
                            tree->BigRIPSIC_fADC[1][i+1]-tree->BigRIPSIC_fADC[1][0]);
            }
        }
        i++;
        if(!((i-range.at(0))%downscale)){
            consolemutex.lock();
            progressbar(i-range.at(0),range.at(1)-range.at(0),threadno);
            consolemutex.unlock();
        }
    }

    //Step 3: rejoin data histograms
    unitemutex.lock();
    for(uint i=0; i<comparediag.size(); i++){
        for(uint j=0;j<comparediag.at(0).size();j++){
            comparediag.at(i).at(j).Add(new TH2D(_comparediag.at(i).at(j)));
        }
    }
    unitemutex.unlock();
}

void iccut::analyse(std::vector<treereader *> tree, TFile *output) {
    threads = (int)tree.size();
    printf("Now beginning analysis of the IC's (Fpl 7 & 11)\n");

    vector <string> readoutkeys{"BigRIPSIC.nhitchannel", "BigRIPSIC.fADC[32]"};
    for(auto &i: tree) i->setloopkeys(readoutkeys);


    for(int i=1; i<numchannel; i++){
        vector<string> arrname = {"ICratio" + to_string(i) + "to0fpl7",
                                  "ICratio" + to_string(i) + "to0fpl11"};
        string arrtitle = "Peak ADC Ratio " + to_string(i) + " to 0";
        comparediag.push_back({TH2D(arrname.at(0).c_str(),arrtitle.c_str(),
                                    2048,0,16384,750,-1500,1500),
                               TH2D(arrname.at(1).c_str(),arrtitle.c_str(),
                                    2048,0,16384,750,-1500,1500)});
        for(auto &elem: comparediag.back()){
            elem.SetOption("colz");
            elem.GetXaxis()->SetTitle("ADC0 [ch]");
            elem.GetYaxis()->SetTitle(("ADC"+to_string(i)+ " - ADC0 [ch]").c_str());
        }
    }
    printf("Successfully generated Histograms for the IC-cut...\n");

    int cutpre = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<int> ranges = {i*goodevents.size()/threads,
                              (i+1)*goodevents.size()/threads-1};
        th.emplace_back(thread(&iccut::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for (auto &i: th) i.join();

    int cutafter = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    printf("\nIC Cut out %i Events %f %%\n", cutpre-cutafter,
           100*(cutpre-cutafter)/(double)goodevents.size());

    output->mkdir("IC/IC7");
    output->mkdir("IC/IC11");
    output->cd("IC/IC7");
    for(auto &elem: comparediag) elem.at(0).Write();
    output->cd("IC/IC11");
    for(auto &elem: comparediag) elem.at(1).Write();
    output->cd("");

}