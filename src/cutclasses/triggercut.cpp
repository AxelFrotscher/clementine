//
// Created by afrotscher on 7/27/18.
//

#include "cutclasses/triggercut.h"

using namespace std;

thread triggercut::innerloop(treereader *tree, std::vector <std::atomic<bool>> &goodevents,
                                  const std::vector<int> range) {
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

void triggercut::analyse(vector<treereader*> tree) {
    threads = (int)tree.size();
    printf("Now performing the Trigger cut with %i threads.", threads);
    // This method aims at analysing the trigger bits
    vector<string> keys{"EventInfo.fBit"};
    for(auto &i:tree) i->setloopkeys(keys);

    vector<thread> th;
    for(uint i=0; i<threads;i++){
        vector<int> ranges = {i*goodevents.size()/threads,
                                 (i+1)*goodevents.size()/threads-1};
        th.emplace_back(thread(&triggercut::innerloop, this,
                               tree.at(i),ref(goodevents), ranges));
    }

    for(auto &t : th) t.join();

    int cutout = (int)accumulate(goodevents.begin(), goodevents.end(), 0.0);
    printf("\nTrigger Cut out %i Events %f %%\n", cutout,
           100*(1-cutout/(double)goodevents.size()));
}