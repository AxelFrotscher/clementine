//
// Created by afrotscher on 8/16/18.
//

#pragma once

#include <vector>
#include <mutex>
#include <atomic>

class progressbar {
public:
    // Regular constructor for updating the database
    progressbar(uint totevent_, uint threadpos_){
        initstatics();
        currentthreads++;
        if(threadpos_ < displaythreads){
            threadpos = threadpos_;
            totevent.at(threadpos_)  = totevent_;
        }
        else worker = false;
    };
    // Second constructor, used by the display thread
    explicit progressbar(){
        worker = false;
        initstatics();
    };
    void draw();
    void increaseevent(){if(worker) currevt[threadpos]++;};
    void reset();
    const bool ongoing(){return ongoinganalysis;};

private:
    void initstatics();

    static std::vector <int> currevt;
    static const int constadd = 10;
    static const int barwidth = 30;
    static const int displaythreads = 7;
    static constexpr double updateinterval = 1; // s
    static std::atomic<int> currentthreads;  // current working (displayed) threads
    static std::mutex currevtmutex;
    bool worker = true ; // Defines a thread to be displayed
    static std::vector<uint> totevent;  // number of object to be calculated
    uint threadpos = 0;  // thread position

    static std::atomic<int> finishercount;  // number of finished threads
    static std::atomic<bool> ongoinganalysis; // switch for status of analysis
};