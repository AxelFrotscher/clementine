//
// Created by afrotscher on 8/16/18.
//

#include "progress.h"
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <iostream>
#include <sstream>
#include <cmath>
#include <TString.h>

using namespace std;

// static is like 'extern' separate declaration needed
vector<int> progressbar::currevt;
vector<uint> progressbar::totevent;
atomic<int> progressbar::finishercount;
mutex progressbar::currevtmutex;
atomic<int> progressbar::currentthreads;
atomic<bool> progressbar::ongoinganalysis;
vector<int> progressbar::lastevent;
long int progressbar::lasttime;
int progressbar::barwidth = 100;

void progressbar::draw(){
    // Worker Instances must not draw anything
    if(worker && (currentthreads>1)) return;

    vector<int> pos; // position for each bar

    for(uint i=0; i<currevt.size(); i++){ // Loop over each bar
        pos.push_back((int)(barwidth*(float)currevt.at(i)/totevent.at(i)));
        if(currevt.at(i) == 0) continue; // Do not display empty threads

        stringstream stream;
        stream << "[";
        for(int j=0; j<barwidth; j++){
            // determine output for each bar
            if(j<pos.at(i))        stream << "=";
            else if (j==pos.at(i)) stream << ">";
            else                   stream << " ";
        }
        stream << "] ";

        string pstr = stream.str(); // Put out string for insertion of speed

        // Prepare time difference and speed
        long int thistime = chrono::duration_cast<chrono::milliseconds>(
                        chrono::steady_clock::now().time_since_epoch()).count();
        long int speed = 0;
        if(thistime-lasttime > 0) speed = 1000*(currevt.at(i)-lastevent.at(i))/
                                          (thistime-lasttime);

        // Replace with speed so that it doesnt block
        string strsp = " " + to_string(speed) + "/s " +
                       Form("%.1f", max((int(1000.*currevt.at(i)/
                       totevent.at(i)))/10.,0.)) + "% ";
        if(barwidth - pos.at(i) > strsp.size()) // Display speed on right side
            pstr.replace(barwidth-strsp.size()+1, strsp.size(),
                         strsp);
        else if(pos.at(i) > (strsp.size()))
            pstr.replace(2, strsp.size(), strsp);

        cout << pstr;
        lastevent.at(i) = currevt.at(i); // Replace value for next draw() call
    }
    cout << "\r";
    cout.flush();

    // Update the "old" values
    lasttime = chrono::duration_cast<chrono::milliseconds>(
            chrono::steady_clock::now().time_since_epoch()).count();

    if(currentthreads>1){
        barwidth = 31;
        this_thread::sleep_for(1s);
    }
}

void progressbar::reset() {
    // Check for last thread and reset if needed
    if(++finishercount == currentthreads){
        progressbar::currevt  = vector<int>(displaythreads, 0);
        progressbar::totevent = vector<uint>(displaythreads, 0);

        // Last initialization of last step time and events
        lasttime = chrono::duration_cast<chrono::milliseconds>(
                chrono::steady_clock::now().time_since_epoch()).count();
        lastevent = vector<int>(displaythreads, 0);

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

        // First initialization of last step time and events
        lasttime = chrono::duration_cast<chrono::milliseconds>(
                    chrono::steady_clock::now().time_since_epoch()).count();
        lastevent = vector<int>(displaythreads, 0);
    }

    // Reset threadcount for subsequent calls
    if(finishercount == currentthreads){
        finishercount.exchange(0);
        currentthreads.exchange(0);
        ongoinganalysis.exchange(true);
    }
    currevtmutex.unlock();
}