//
// Created by afrotscher on 7/27/18.
//

#include "cutclasses/triggercut.h"

using namespace std;

void triggercut::innerloop(treereader *tree, std::vector <std::atomic<bool>> &goodevents,
                                  const std::vector<uint> range) {
    int i = range.at(0);
    int threadno = range.at(0)/(range.at(1)-range.at(0));
    const int downscale = (range.at(1)-range.at(0))/100; // every percent

    while(i<range.at(1)){
        if(goodevents.at(i)){
            tree->getevent(i);
            if(tree->EventInfo_fBit[0] == badtrg){
                goodevents.at(i).exchange(false);
            }
        }
        i++;
        if(!((i-range.at(0))%downscale)){
            consolemutex.lock();
            progressbar(i-range.at(0), range.at(1)-range.at(0),threadno);
            consolemutex.unlock();
        }
    }
}

void triggercut::analyse(const vector<string> input){
    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    printf("Now performing the Trigger cut with %i threads.", threads);

    // This method aims at analysing the trigger bits
    vector<string> keys{"EventInfo.fBit"};
    for(auto &i:tree) i->setloopkeys(keys);

    vector<thread> th;
    for(uint i=0; i<threads;i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                              (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&triggercut::innerloop, this,
                               tree.at(i),ref(goodevents), ranges));
    }

    for(auto &t : th) t.join();

    int cutout = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);
    printf("\nTrigger Cut out %lu Events %f %%\n", goodevents.size()-cutout,
           100*(1-cutout/(double)goodevents.size()));
}