//
// Created by afrotscher on 8/16/18.
//

#include "progress.h"
#include <cstdlib>
#include <chrono>
#include <thread>
#include <iostream>
#include <sstream>
#include <cmath>
#include <TString.h>

using std::stringstream, std::string, std::max, std::to_string, std::cout,
      std::vector, std::chrono::duration_cast, std::chrono::steady_clock,
      std::chrono::milliseconds, std::this_thread::sleep_for,
      std::chrono_literals::operator""s;

void progressbar::draw(){
    /// Worker Instances must not draw anything
    if(worker && (currentthreads>1)) return;

    vector<unsigned long> pos; // position for each bar

    for(ulong i=0; i<currevt.size(); i++){ // Loop over each bar
        pos.push_back((int)(barwidth*(float)currevt.at(i)/totevent.at(i)));
        if(currevt.at(i) == 0) continue; // Do not display empty threads

        stringstream stream;
        stream << "[";
        for(unsigned long j=0; j<barwidth; j++){
            // determine output for each bar
            if(j<pos.at(i))        stream << "=";
            else if (j==pos.at(i)) stream << ">";
            else                   stream << " ";
        }
        stream << "] ";

        string pstr = stream.str(); // Put out string for insertion of speed

        // Prepare time difference and speed
        long int thistime = duration_cast<milliseconds>(
                                steady_clock::now().time_since_epoch()).count();
        long int speed = 0;
        if(thistime-lasttime > 0) speed = 1000*(currevt.at(i)-lastevent.at(i))/
                                          (thistime-lasttime);

        // Replace with speed so that it doesnt block
        string strsp = " " + to_string(speed) + "/s " +
                       Form("%.1f", max((int(1000.*currevt.at(i)/
                       totevent.at(i)))/10.,0.)) + "% ";

        if(pos.at(i) > (strsp.size()))
            pstr.replace(2, strsp.size(), strsp);
        else if(barwidth - pos.at(i) > strsp.size()) // Display speed on right side
            pstr.replace(barwidth-strsp.size()+1, strsp.size(), strsp);

        cout << pstr;
        lastevent.at(i) = currevt.at(i); // Replace value for next draw() call
    }
    cout.flush();
    cout << "\r"; // for remote execution AFTER the flushing...

    // Update the "old" values
    lasttime = duration_cast<milliseconds>(
                                steady_clock::now().time_since_epoch()).count();

    if(currentthreads>1 || !worker){
        barwidth = 31;
        sleep_for(1s);
    }
}

void progressbar::reset() {
    /// Check for last thread and reset if needed
    if(++finishercount == currentthreads){
        progressbar::currevt  = vector<int>(displaythreads, 0);
        progressbar::totevent = vector<uint>(displaythreads, 0);

        // Last initialization of last step time and events
        lasttime = duration_cast<milliseconds>(
                                steady_clock::now().time_since_epoch()).count();
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
        lasttime = duration_cast<milliseconds>(
                                steady_clock::now().time_since_epoch()).count();
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