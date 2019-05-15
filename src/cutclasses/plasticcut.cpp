//
// Created by afrotscher on 8/2/18.
//

#include "cutclasses/plasticcut.h"
#include "progress.h"
#include "TMath.h"
#include <thread>
#include <numeric>
#include "zdssetting.h"

using std::vector, std::atomic, std::string, std::to_string, std::thread;

void plasticcut::innerloop(treereader &tree, vector<uint> range) {
    // Step 1: cloning histograms
    decltype(qcorr)    _qcorr;
    decltype(qcorr2D)  _qcorr2D;
    decltype(tqcorr2D) _tqcorr2D;
    decltype(qxF11)    _qxF11;

    for(auto &i : qcorr)    _qcorr.emplace_back(i);
    for(auto &i : qcorr2D)  _qcorr2D.emplace_back(i);
    for(auto &i : tqcorr2D) _tqcorr2D.emplace_back(i);
    for(auto &i : qxF11)    _qxF11.emplace_back(i);

    // Step 2: Preparing Variables
    uint threadno = range.at(0)/(range.at(1)-range.at(0));

    progressbar progress(range.at(1)-range.at(0), threadno);

    bool extraF11analyse = false;
    if(mycut.at(threadno).size() == 7) extraF11analyse = true;

    for(int i=range.at(0);i < range.at(1);i++){
        if(goodevents.at(i).at(0)){
            tree.getevent(i);

            //F7 banana-cut, "uncut histogram's
            if(mycut.at(threadno).at(1)->IsInside(
                    tree.BigRIPSPlastic_fTLRaw[1]-tree.BigRIPSPlastic_fTRRaw[1],
                    TMath::Log(tree.BigRIPSPlastic_fQLRaw[1]/
                                       (double)tree.BigRIPSPlastic_fQRRaw[1]))){
                for(uint j=0; j<numplastic; j++){              // plastic loop
                    if((tree.BigRIPSPlastic_fQLRaw[j] >0 )&&
                       (tree.BigRIPSPlastic_fQRRaw[j] >0)){   // 0 charge veto
                        _qcorr2D.at(2*j).Fill(tree.BigRIPSPlastic_fQLRaw[j],
                                           tree.BigRIPSPlastic_fQRRaw[j]);
                        _qcorr.at(2*j).Fill(sqrt(tree.BigRIPSPlastic_fQLRaw[j]*
                                              tree.BigRIPSPlastic_fQRRaw[j]));
                        if((tree.BigRIPSPlastic_fTLRaw[j] >0) &&
                           (tree.BigRIPSPlastic_fTRRaw[j] >0) ){ // 0 time veto
                            _tqcorr2D.at(2*j).Fill(tree.BigRIPSPlastic_fTLRaw[j]-
                                                  tree.BigRIPSPlastic_fTRRaw[j],
                                       TMath::Log(tree.BigRIPSPlastic_fQLRaw[j]/
                                        (double)tree.BigRIPSPlastic_fQRRaw[j]));
                        }
                    }
                }
            }

            //Cut out events properly:
            if(!(mycut.at(threadno).at(0)->IsInside(
                    tree.BigRIPSPlastic_fTLRaw[0]-tree.BigRIPSPlastic_fTRRaw[0],
                    TMath::Log(tree.BigRIPSPlastic_fQLRaw[0]/
                               (double)tree.BigRIPSPlastic_fQRRaw[0])) &&
                 mycut.at(threadno).at(1)->IsInside(
                    tree.BigRIPSPlastic_fTLRaw[1]-tree.BigRIPSPlastic_fTRRaw[1],
                    TMath::Log(tree.BigRIPSPlastic_fQLRaw[1]/
                               (double)tree.BigRIPSPlastic_fQRRaw[1])))){

                // If not in F3 and F7, cut out event completely
                for(auto &k:goodevents.at(i)) k.exchange(false);
            }
            else{ // Test F7 and F11, no because we want cross sections
                for(uint j=0; j<numplastic; j++){              // plastic loop
                    if(!(tree.BigRIPSPlastic_fTLRaw[j] >0 && // 0 time veto
                         tree.BigRIPSPlastic_fTRRaw[j] >0) ) continue;
                    _qcorr2D.at(2*j+1).Fill(tree.BigRIPSPlastic_fQLRaw[j],
                                            tree.BigRIPSPlastic_fQRRaw[j]);
                    _qcorr.at(2*j+1).Fill(sqrt(tree.BigRIPSPlastic_fQLRaw[j]*
                                               tree.BigRIPSPlastic_fQRRaw[j]));
                    _tqcorr2D.at(2*j+1).Fill(tree.BigRIPSPlastic_fTLRaw[j]-
                                             tree.BigRIPSPlastic_fTRRaw[j],
                                             TMath::Log(tree.BigRIPSPlastic_fQLRaw[j]/
                                             (double)tree.BigRIPSPlastic_fQRRaw[j]));

                }
                // Extra special Obertelli plastic F11 (index=3) analysis
                if(!extraF11analyse) continue;
                for(int j=4; j<7; j++){
                    if(mycut.at(threadno).at(j)->IsInside(tree.BigRIPSPlastic_fQLRaw[3],
                                                      tree.BigRIPSPlastic_fQRRaw[3]))
                    _qxF11.at(j-4).Fill(tree.F11X);
                }
            }
        }
        progress.increaseevent();
    }

    // Step 3 reuniting the diagrams
    unitemutex.lock();
    for(uint i=0; i<_qcorr.size(); i++) qcorr.at(i).Add(&_qcorr.at(i));
    for(uint i=0; i<_qxF11.size(); i++) qxF11.at(i).Add(&_qxF11.at(i));
    for(uint i=0; i<_qcorr2D.size(); i++) qcorr2D.at(i).Add(&_qcorr2D.at(i));
    for(uint i=0; i<_tqcorr2D.size(); i++) tqcorr2D.at(i).Add(&_tqcorr2D.at(i));
    unitemutex.unlock();

    progressbar::reset();
}

void plasticcut::analyse(const vector<string> &input, TFile *output) {

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Performing Plastic Cut with %i threads.\n", threads);

    // Setting appropriate keys
    vector<string> keys{"BigRIPSPlastic.fQLRaw", "BigRIPSPlastic.fQRRaw",
                        "BigRIPSPlastic.fTLRaw", "BigRIPSPlastic.fTRRaw",
                        "BigRIPSPlastic.fpl", "F11X"};
    for(auto &i: tree) i.setloopkeys(keys);

    // Generating the histograms:
    tree.at(0).getevent(0);
    for(uint i=0; i<numplastic; i++){
        arrayname.push_back(vector<string>{
            "f7pltrigQ."  + to_string(tree.at(0).BigRIPSPlastic_fpl[i]),
            "PlasticQ2D." + to_string(tree.at(0).BigRIPSPlastic_fpl[i]),
            "tQcorr."     + to_string(tree.at(0).BigRIPSPlastic_fpl[i])
        });

        arraytitle.push_back(vector<string>{
            "Charge deposition by Plastic F"     +
                                   to_string(tree.at(0).BigRIPSPlastic_fpl[i]),
            "Charge distribution by Plastic F"   +
                                   to_string(tree.at(0).BigRIPSPlastic_fpl[i]),
            "Charge Ratio/time difference Pl. F" +
                                    to_string(tree.at(0).BigRIPSPlastic_fpl[i])
        });

        qcorr2D.emplace_back(arrayname.at(i).at(1).c_str(),
                        arraytitle.at(i).at(1).c_str(), 700,0,4500, 700,0,4500);
        qcorr2D.emplace_back((arrayname.at(i).at(1) + "cut").c_str(),
                             arraytitle.at(i).at(1).c_str(), 700,0,4500, 700,0,4500);
        qcorr.emplace_back(arrayname.at(i).at(0).c_str(),
                           arraytitle.at(i).at(0).c_str(), 2500,0,5000);
        qcorr.emplace_back( (arrayname.at(i).at(0)+ "cut").c_str(),
                             arraytitle.at(i).at(0).c_str(), 2500,0,5000);
        tqcorr2D.emplace_back(arrayname.at(i).at(2).c_str(),
                        arraytitle.at(i).at(2).c_str(), 450,-250,200, 500,-2.5,2.5);
        tqcorr2D.emplace_back((arrayname.at(i).at(2)+"cut").c_str(),
                              arraytitle.at(i).at(2).c_str(), 450,-250,200, 500,-2.5,2.5);
    }
    qxF11.emplace_back("F11ped", "F11 pedestal plastic position", 80,-40,40);
    qxF11.emplace_back("F11reg", "F11 regular plastic position", 80,-40,40);
    qxF11.emplace_back("F11wei", "F11 weird plastic position", 80,-40,40);


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

    //Get plasticcuts
    for(int i=0; i<threads; i++){
        mycut.push_back(setting::getplasticcut());
    }

    progressbar finishcondition;
    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&plasticcut::innerloop, this, std::ref(tree.at(i)),
                               ranges));
    }

    for (auto &i: th) i.detach();

    while(progressbar::ongoing()) finishcondition.draw();

    vector<string> folders{"Plastics/2D", "Plastics/Q1Q2", "Plastics/TQCorr",
                           "Plastics/F11weird"};
    for (auto &str:folders) output->mkdir(str.c_str());

    output->cd("Plastics/2D");
    for(auto &histo: qcorr2D) histo.Write();
    output->cd("Plastics/Q1Q2");
    for(auto &histo: qcorr) histo.Write();
    output->cd("Plastics/TQCorr");
    for(auto &histo: tqcorr2D) histo.Write();
    output->cd("Plastics/F11weird");
    for (auto &histo: qxF11) histo.Write();
    output->cd("");

    vector<int> cutafter ={0,0};
    for(auto &i: goodevents){
        cutafter.at(0) += i.at(0);
        cutafter.at(1) += i.at(1);
    }

    printf("\nPlastic Cut out %i F1-7 Events (%.3f %%) %i F1-11 Events (%.3f %%) \n",
           cutpre.at(0)-cutafter.at(0), 100*(cutpre.at(0)-cutafter.at(0))/(double)goodevents.size(),
           cutpre.at(1)-cutafter.at(1), 100*(cutpre.at(1)-cutafter.at(1))/(double)goodevents.size());

    for(auto &I: mycut) for(auto &j:I) delete j;
}
