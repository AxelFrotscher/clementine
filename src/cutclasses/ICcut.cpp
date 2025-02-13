//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/ICcut.h"
#include <thread>
#include "progress.h"
#include "zdssetting.h"

using std::vector, std::atomic, std::string, std::thread, std::to_string;

void iccut::innerloop(treereader &tree, vector<int> range) {
    ///Step 1: Cloning histograms
    decltype(comparediag) _comparediag(comparediag);

    ///Step 2: Preparing Variables
    const uint threadno = range.at(0)/(range.at(1)-range.at(0));

    progressbar progress(range.at(1)-range.at(0),threadno);
    //setting set;
    //const bool trans = setting::isemptyortrans();
    vector<double> limits = setting::geticlimits();
    assert(limits.size() == 4);

    for(int i=range.at(0); i < range.at(1); i++){
        if(goodevents.at(i).at(0)){
            // Determine cut only on good events
            tree.getevent(i);

            vector<bool> inrange = {false, false};
            for(int j=0; j<numchannel; j++){
                if(tree.BigRIPSIC_fADC[0][j] > limits.at(0) &&
                   tree.BigRIPSIC_fADC[0][j] < limits.at(1)){
                    inrange.at(0) = true;
                }
                if(tree.BigRIPSIC_fADC[1][j] > limits.at(2) &&
                   tree.BigRIPSIC_fADC[1][j] < limits.at(3)){
                    inrange.at(1) = true;
                }
            }

            if(!inrange.at(0)) for(auto &j:goodevents.at(i)) j.exchange(false);
            if(!inrange.at(1)) goodevents.at(i).at(1).exchange(false);

            if((tree.BigRIPSIC_fADC[0][0] <0) ||
               (tree.BigRIPSIC_fADC[1][0] <0)) continue;
            for(int j=0; j<(numchannel-1); j++){
                if(tree.BigRIPSIC_fADC[0][j] >0)
                    _comparediag.at(j).at(0).Fill(
                            tree.BigRIPSIC_fADC[0][0],
                            tree.BigRIPSIC_fADC[0][j+1]- tree.BigRIPSIC_fADC[0][0]);
                if(tree.BigRIPSIC_fADC[1][j] >0)
                    _comparediag.at(j).at(1).Fill(
                            tree.BigRIPSIC_fADC[1][0],
                            tree.BigRIPSIC_fADC[1][j+1]-tree.BigRIPSIC_fADC[1][0]);
            }
        }
        progress.increaseevent();
    }

    //Step 3: rejoin data histograms
    unitemutex.lock();
    for(unsigned long i=0; i<comparediag.size(); i++){
        for(unsigned long j=0; j<comparediag.at(0).size();j++){
            comparediag.at(i).at(j).Add(&_comparediag.at(i).at(j));
        }
    }
    unitemutex.unlock();

    progressbar::reset();
}

void iccut::analyse(const std::vector<std::string> &input, TFile *output) {

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Now beginning analysis of the IC's (Fpl 7 & 11)\n");

    vector <string> readoutkeys{"BigRIPSIC.nhitchannel", "BigRIPSIC.fADC[32]"};
    for(auto &i: tree) i.setloopkeys(readoutkeys);

    for(int i=1; i<numchannel; i++){
        comparediag.emplace_back();
        vector<string> arrname = {"ICratio" + to_string(i) + "to0fpl7",
                                  "ICratio" + to_string(i) + "to0fpl11"};
        string arrtitle = "Peak ADC Ratio " + to_string(i) + " to 0";

        comparediag.back() =
         {{arrname.at(0).c_str(),arrtitle.c_str(),1024,0,16384,375,-1500,1500},
          {arrname.at(1).c_str(),arrtitle.c_str(),1024,0,16384,375,-1500,1500}};

        for(auto &elem: comparediag.back()){
            elem.SetOption("colz");
            elem.GetXaxis()->SetTitle("ADC0 [ch]");
            elem.GetYaxis()->SetTitle(("ADC"+to_string(i)+ " - ADC0 [ch]").c_str());
        }
    }
    printf("Successfully generated Histograms for the IC-cut...\n");

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
        th.emplace_back(thread(&iccut::innerloop, this, std::ref(tree.at(i)),
                               ranges));
    }

    for (auto &i: th) i.detach();

    while(progressbar::ongoing()) finishcondition.draw();

    vector<int> cutafter = {0,0};
    for(auto &i:goodevents){
        cutafter.at(0) += i.at(0);
        cutafter.at(1) += i.at(1);
    }

    printf("\nIC Cut out %i F1-7 Events (%.2f %%) %i F1-11 Events (%.2f %%) \n",
           cutpre.at(0)-cutafter.at(0), 100*(cutpre.at(0)-cutafter.at(0))/(double)goodevents.size(),
           cutpre.at(1)-cutafter.at(1), 100*(cutpre.at(1)-cutafter.at(1))/(double)goodevents.size());

    output->mkdir("IC/IC7");
    output->mkdir("IC/IC11");
    output->cd("IC/IC7");
    for(auto &elem: comparediag) elem.at(0).Write();
    output->cd("IC/IC11");
    for(auto &elem: comparediag) elem.at(1).Write();
    output->cd("");
}