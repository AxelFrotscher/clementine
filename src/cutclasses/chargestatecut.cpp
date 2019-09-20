//
// Created by afrotscher on 8/1/18.
//

#include "cutclasses/chargestatecut.h"
#include "libconstant.h"
#include "zdssetting.h"
#include <thread>

using std::vector, std::string, std::thread, std::atomic, std::endl;

void ccsc::innerloop(treereader &tree, std::vector<int> range) {
    // precious tight inner loop
    // Cloning histograms
    decltype(cschist) _cschist(cschist);

    // Step 2: Preparing variables
    const uint threadno = range.at(0)/(range.at(1)-range.at(0));
    
    double brhoratio = 0;

    progressbar progress(range.at(1)-range.at(0), threadno);

    for(int i=range.at(0); i<range.at(1); i++){
         if(goodevents.at(i).at(1)){
            tree.getevent(i);

            brhoratio = tree.BigRIPSBeam_brho[3]/tree.BigRIPSBeam_brho[2];
            _cschist.at(0).Fill(brhoratio, tree.BigRIPSBeam_brho[2]);

            if(mycut.at(threadno) &&
               mycut.at(threadno)->IsInside(brhoratio,tree.BigRIPSBeam_brho[2])){
                _cschist.at(1).Fill(brhoratio,tree.BigRIPSBeam_brho[2]);
            }
            else goodevents.at(i).at(1).exchange(false);
        }

        progress.increaseevent();
    }

    // Step 3: rejoin the data
    unitemutex.lock();
    for(unsigned long i=0; i<_cschist.size();i++) cschist.at(i).Add(&_cschist.at(i));
    unitemutex.unlock();
    progressbar::reset();
}

void ccsc::analyse(const std::vector<std::string> &input, TFile* output){

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Beginning with the CCSC. %i threads.\n", threads);

    // Generate output histogram
    cschist = {{"csc", "Charged state change", 500,0.9,1.1,1000,3,7},
               {"csccut", "Charged state change cut", 500,0.95,1.05,1000,3,7}};

    for(auto &histo: cschist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("B#rho F8-9");
        histo.GetXaxis()->SetTitle("B#rho F8-9 / B#rho F9-11");
    }

    // Get relevant keys
    vector<string> keys{"BigRIPSBeam.brho"};
    for(auto &i: tree) i.setloopkeys(keys);

    // Get Cut from right setting | mode
    // setting set;
    //if(set.isemptyortrans()){ // Cut after F7 not sensible for empty/trans runs
    //    cout << "Empty or trans run. Not doing CCSC cut." << endl; return;
    //}
    for([[gnu::unused]] auto &i: tree) mycut.push_back(setting::getbrhocut());

    int cutpre = 0;
    for(auto &i:goodevents) cutpre += i.at(1);

    progressbar finishcondition;
    vector<thread> th;
    for(int i=0; i<threads; i++){
        vector<int> ranges = {(int)(i*goodevents.size()/threads),
                               (int)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&ccsc::innerloop, this, std::ref(tree.at(i)),
                               ranges));
    }

    for (auto &t: th) t.detach();

    // Setup synchronization class
    while(progressbar::ongoing()) finishcondition.draw();

    int cutafter = 0;
    for(auto &i:goodevents) cutafter += i.at(1);

    printf("\nCCSC Cut out (F1-F11) %i Events %.3f %%\n", cutpre-cutafter,
           100.*(cutpre-cutafter)/(double)goodevents.size());

    output->mkdir("CSC");
    output->cd("CSC");
    for(auto hist: cschist) hist.Write();
    output->cd("");

    for(auto &I: mycut) delete I;
}