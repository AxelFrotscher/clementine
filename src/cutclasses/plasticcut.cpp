//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/plasticcut.h"
#include "progress.h"
#include "TMath.h"
#include <thread>
#include <numeric>
#include "zdssetting.h"

using namespace std;

void plasticcut::innerloop(treereader *tree, std::vector<std::vector<std::atomic<bool>>>
                             &goodevents, std::vector<uint> range) {
    // Step 1: cloning histograms
    vector<TH1D> _qcorr;
    vector<TH2D> _qcorr2D;
    vector<TH2D> _tqcorr2D;

    for(auto &i : qcorr)    _qcorr.emplace_back(TH1D(i));
    for(auto &i : qcorr2D)  _qcorr2D.emplace_back(TH2D(i));
    for(auto &i : tqcorr2D) _tqcorr2D.emplace_back(TH2D(i));

    // Step 2: Preparing Variables
    uint threadno = range.at(0)/(range.at(1)-range.at(0));

    progressbar progress(range.at(1)-range.at(0), threadno);

    for(int i=range.at(0);i < range.at(1);i++){
        if(goodevents.at(i).at(0)){
            tree->getevent(i);
            if(sqrt(tree->BigRIPSPlastic_fQLRaw[F7pos]*
                    tree->BigRIPSPlastic_fQRRaw[F7pos])> threshhold) // F7-cut
                for(uint j=0; j<numplastic; j++){              // plastic loop
                    if((tree->BigRIPSPlastic_fQLRaw[j] >0 )&&
                       (tree->BigRIPSPlastic_fQRRaw[j] >0)){   // 0 charge veto
                        qcorr2D.at(j).Fill(tree->BigRIPSPlastic_fQLRaw[j],
                                           tree->BigRIPSPlastic_fQRRaw[j]);
                        qcorr.at(j).Fill(sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                                              tree->BigRIPSPlastic_fQRRaw[j]));
                        if((tree->BigRIPSPlastic_fTLRaw[j] >0) &&
                           (tree->BigRIPSPlastic_fTRRaw[j] >0) ){ // 0 time veto
                            tqcorr2D.at(j).Fill(tree->BigRIPSPlastic_fTLRaw[j]-
                                                tree->BigRIPSPlastic_fTRRaw[j],
                                      TMath::Log(tree->BigRIPSPlastic_fQLRaw[j]/
                                       (double)tree->BigRIPSPlastic_fQRRaw[j]));
                        }
                    }
                }
            bool temp = true;
            // Applying cut to data
            for(int j=0;j<numplastic;j++){
                temp = temp * (sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                       tree->BigRIPSPlastic_fQRRaw[j]) > acceptance_range[j][0])*
                       (sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                       tree->BigRIPSPlastic_fQRRaw[j]) < acceptance_range[j][1]);
                if(!temp && j==2){
                    for(auto &j:goodevents.at(i)) j.exchange(false);
                    break;
                }
            }
            if(!temp) goodevents.at(i).at(1).exchange(false);
        }
        progress.increaseevent();
    }

    // Step 3 reuniting the diagrams
    unitemutex.lock();
    for(uint i=0;i<_qcorr.size();i++) qcorr.at(i).Add(new TH1D(_qcorr.at(i)));
    for(uint i=0;i<_qcorr2D.size();i++) qcorr2D.at(i).Add(new TH2D(_qcorr2D.at(i)));
    for(uint i=0;i<_tqcorr2D.size();i++) tqcorr2D.at(i).Add(new TH2D(_tqcorr2D.at(i)));
    unitemutex.unlock();

    progress.reset();
}

void plasticcut::analyse(const std::vector<std::string> input, TFile *output) {
    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto &h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    printf("Performing Plastic Cut with %i threads.\n", threads);

    // Setting appropriate keys
    vector<string> keys{"BigRIPSPlastic.fQLRaw", "BigRIPSPlastic.fQRRaw",
                        "BigRIPSPlastic.fTLRaw", "BigRIPSPlastic.fTRRaw",
                        "BigRIPSPlastic.fpl"};
    for(auto &i: tree) i->setloopkeys(keys);

    // Generating the histograms:
    tree.at(0)->getevent(0);
    for(uint i=0; i<numplastic; i++){
        arrayname.push_back(vector<string>{
            "f7pltrigQ."  + to_string(tree.at(0)->BigRIPSPlastic_fpl[i]),
            "PlasticQ2D." + to_string(tree.at(0)->BigRIPSPlastic_fpl[i]),
            "tQcorr."     + to_string(tree.at(0)->BigRIPSPlastic_fpl[i])
        });

        arraytitle.push_back(vector<string>{
            "Charge deposition by Plastic F"     +
                                   to_string(tree.at(0)->BigRIPSPlastic_fpl[i]),
            "Charge distribution by Plastic F"   +
                                   to_string(tree.at(0)->BigRIPSPlastic_fpl[i]),
            "Charge Ratio/time difference Pl. F" +
                                    to_string(tree.at(0)->BigRIPSPlastic_fpl[i])
        });

        qcorr2D.emplace_back(
            TH2D(arrayname.at(i).at(1).c_str(), arraytitle.at(i).at(1).c_str(),
                 500,0,1500, 500,0,1500));
        qcorr.emplace_back(
            TH1D(arrayname.at(i).at(0).c_str(), arraytitle.at(i).at(0).c_str(),
                 750,0,1500));
        tqcorr2D.emplace_back(
            TH2D(arrayname.at(i).at(2).c_str(), arraytitle.at(i).at(2).c_str(),
                 300,-150,150, 400,-2,2));
    }

    for(auto &i : qcorr){
        i.GetXaxis()->SetTitle("#sqrt{Q_{l}Q_{r}} [ch]");
        i.GetYaxis()->SetTitle("N");
    }

    for(auto &i : qcorr2D){
        i.GetXaxis()->SetTitle("Q_{l} [ch]");
        i.GetYaxis()->SetTitle("Q_{r} [ch]");
        i.SetOption("colz");
    }

    for(auto &i : tqcorr2D){
        i.GetXaxis()->SetTitle("t_{l}-t_{r} [ns]");
        i.GetYaxis()->SetTitle("ln(Q_{l}/Q_{r})");
        i.SetOption("colz");
    }

    // Start the looop
    vector<int> cutpre = {0,0};
    for(auto &i: goodevents){
        cutpre.at(0) += i.at(0);
        cutpre.at(1) += i.at(1);
    }

    progressbar finishcondition;
    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&plasticcut::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for (auto &i: th) i.detach();

    while(finishcondition.ongoing()) finishcondition.draw();

    vector<string> folders{"Plastics/2D", "Plastics/Q1Q2", "Plastics/TQCorr"};
    for (auto &str:folders) output->mkdir(str.c_str());

    output->cd("Plastics/2D");
    for(auto histo: qcorr2D) histo.Write();
    output->cd("Plastics/Q1Q2");
    for(auto histo: qcorr) histo.Write();
    output->cd("Plastics/TQCorr");
    for(auto histo: tqcorr2D) histo.Write();
    output->cd("");

    vector<int> cutafter ={0,0};
    for(auto &i: goodevents){
        cutafter.at(0) += i.at(0);
        cutafter.at(1) += i.at(1);
    }

    printf("\nPlastic Cut out %i F1-7 Events (%.3f %%) %i F1-11 Events (%.3f %%) \n",
           cutpre.at(0)-cutafter.at(0), 100*(cutpre.at(0)-cutafter.at(0))/(double)goodevents.size(),
           cutpre.at(1)-cutafter.at(1), 100*(cutpre.at(1)-cutafter.at(1))/(double)goodevents.size());
}
