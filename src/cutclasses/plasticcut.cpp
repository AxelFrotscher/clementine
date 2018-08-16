//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/plasticcut.h"
#include "progress.h"
#include "TMath.h"
#include <thread>

using namespace std;

void plasticcut::innerloop(treereader *tree, std::vector<std::atomic<bool>>
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
    uint i = range.at(0); // counting variable

    progressbar progress(range.at(1)-range.at(0), threadno);

    while(i < range.at(1)){
        if(goodevents.at(i)){
            tree->getevent(i);
            if(sqrt(tree->BigRIPSPlastic_fQLRaw[F7pos]*
                    tree->BigRIPSPlastic_fQRRaw[F7pos])> threshhold) // F7-cut
                for(uint i=0; i<numplastic; i++){              // plastic loop
                    if((tree->BigRIPSPlastic_fQLRaw[i] >0 )&&
                       (tree->BigRIPSPlastic_fQRRaw[i] >0)){   // 0 charge veto
                        qcorr2D.at(i).Fill(tree->BigRIPSPlastic_fQLRaw[i],
                                           tree->BigRIPSPlastic_fQRRaw[i]);
                        qcorr.at(i).Fill(sqrt(tree->BigRIPSPlastic_fQLRaw[i]*
                                              tree->BigRIPSPlastic_fQRRaw[i]));
                        if((tree->BigRIPSPlastic_fTLRaw[i] >0) &&
                           (tree->BigRIPSPlastic_fTRRaw[i] >0) ){ // 0 time veto
                            tqcorr2D.at(i).Fill(tree->BigRIPSPlastic_fTLRaw[i]-
                                                tree->BigRIPSPlastic_fTRRaw[i],
                                      TMath::Log(tree->BigRIPSPlastic_fQLRaw[i]/
                                       (double)tree->BigRIPSPlastic_fQRRaw[i]));
                        }
                    }
                }
            bool temp = true;
            for(int j=0;j<numplastic;j++){
                temp = temp * (sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                       tree->BigRIPSPlastic_fQRRaw[j]) > acceptance_range[j][0])*
                       (sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                       tree->BigRIPSPlastic_fQRRaw[j]) < acceptance_range[j][1]);
            }
            if(!temp) goodevents.at(i).exchange(false);

        }
        i++;
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

    printf("Performing  Plastic Cut with %i threads.\n", threads);

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
    int cutpre = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&plasticcut::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for (auto &i: th) i.detach();

    progressbar finishcondition;
    while(finishcondition.ongoing()) finishcondition.draw();

    int cutafter = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);

    printf("\nPlastic Cut out %i Events %f %%\n", cutpre-cutafter,
           100*(cutpre-cutafter)/(double)goodevents.size());

    vector<string> folders{"Plastics/2D", "Plastics/Q1Q2", "Plastics/TQCorr"};
    for (auto &str:folders) output->mkdir(str.c_str());

    output->cd("Plastics/2D");
    for(auto histo: qcorr2D) histo.Write();
    output->cd("Plastics/Q1Q2");
    for(auto histo: qcorr) histo.Write();
    output->cd("Plastics/TQCorr");
    for(auto histo: tqcorr2D) histo.Write();
    output->cd("");
    printf("Finished writing plastic histograms!\n");

}
