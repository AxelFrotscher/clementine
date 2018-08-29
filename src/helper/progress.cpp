//
// Created by afrotscher on 8/16/18.
//

#include "progress.h"
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <iostream>

using namespace std;

// static is like 'extern' separate declaration needed
vector<int> progressbar::currevt;
vector<uint> progressbar::totevent;
atomic<int> progressbar::finishercount;
mutex progressbar::currevtmutex;
atomic<int> progressbar::currentthreads;
atomic<bool> progressbar::ongoinganalysis;


void progressbar::draw(){
    // Worker Instances must not draw anything
    if(worker) return;

    vector<int> pos; // position for each bar

    for(uint i=0; i<currevt.size(); i++){ // Loop over each bar
        pos.push_back((int)(i*(barwidth+constadd)+
                            barwidth*(float)currevt.at(i)/totevent.at(i)));
        cout << "[";
        for(int j=i*(barwidth+constadd); j<(i*(barwidth+constadd)+barwidth); j++){
            // determine output for each bar
            if(j<pos.at(i))        cout << "=";
            else if (j==pos.at(i)) cout << ">";
            else                   cout << " ";
        }
        cout << "] " << max(int(100.*currevt.at(i)/totevent.at(i)),0) << "% ";
    }
    cout << "\r";
    cout.flush();

    this_thread::__sleep_for(chrono::seconds((long)updateinterval), chrono::nanoseconds(0));
}

void progressbar::reset() {
    // Check for last thread and reset if needed
    if(++finishercount == currentthreads){
        progressbar::currevt  = vector<int>(displaythreads, 0);
        progressbar::totevent = vector<uint>(displaythreads, 0);

        // After last reset the analysis is finished completely
        ongoinganalysis.exchange(false);
    }
}

void progressbar::initstatics() {
    currevtmutex.lock();
    if(currevt.size() != displaythreads){
        currevt  = std::vector<int>(displaythreads, 0);
        totevent = std::vector<uint>(displaythreads, 0);
        finishercount.exchange(0);
        ongoinganalysis.exchange(true);
    }

    // Reset threadcount for subsequent calls
    if(finishercount == currentthreads){
        finishercount.exchange(0);
        currentthreads.exchange(0);
        ongoinganalysis.exchange(true);
    }
    currevtmutex.unlock();
}