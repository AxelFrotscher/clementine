//
// Created by afrotscher on 7/27/18.
//

#include "cutclasses/triggercut.h"
#include "time.h"
#include "progress.h"
#include "zdssetting.h"
#include "txtwriter.h"
#include <thread>
#include <numeric>
#include <deque>

using namespace std;

void triggercut::innerloop(treereader &tree,
                           std::vector<std::vector <std::atomic<bool>>> &goodevents,
                           const std::vector<uint> range) {
    uint i = range.at(0);
    uint threadno = range.at(0)/(range.at(1)-range.at(0));

    // Construct progressbar object
    progressbar progress(range.at(1)-range.at(0),threadno);

    while(i<range.at(1)){
        if(goodevents.at(i).at(0)){
            tree.getevent(i);
            if(tree.EventInfo_fBit[0] == badtrg || tree.EventInfo_fBit[0] == 11){
                for(auto &j:goodevents.at(i))  j.exchange(false);
            }
        }
        i++;
        progress.increaseevent();
    }

    progress.reset();
}

void triggercut::analyse(const vector<string> input){
    setting set;
    txtwriter txt;
    
    if(!set.isemptyortrans() && set.getminos()){
        cout << "Physics Run. Omitting trigger cut to gain statistics" << endl;
        txt.addline("No trigger cut on DALI Trigger applied.");
        return;
    }

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Now performing the Trigger cut with %i threads.\n", threads);

    // This method aims at analysing the trigger bits
    vector<string> keys{"EventInfo.fBit"};
    for(auto &i:tree) i.setloopkeys(keys);

    progressbar finishcondition;
    vector<thread> th;
    for(uint i=0; i<threads;i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                              (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&triggercut::innerloop, this,
                               ref(tree.at(i)),ref(goodevents), ranges));
    }

    for(auto &t : th) t.detach();

    // Setup synchronization class
    while(finishcondition.ongoing()) finishcondition.draw();

    int cutout = 0;
    for(auto &i:goodevents) cutout += i.at(1);

    printf("\nTrigger Cut out %lu Events %.2f %%\n", goodevents.size()-cutout,
           100*(1-cutout/(double)goodevents.size()));

    txt.addline(Form("Trigger Cut for DALI %i events (%.2f %%).", cutout,
                     100*(1-cutout/(double)goodevents.size())));
}