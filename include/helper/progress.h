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
    static void reset();
    static const bool ongoing(){return ongoinganalysis;};

private:
    static void initstatics();

    inline static std::vector <int> currevt; // current progress on the matter

    inline static std::vector <int> lastevent; // progress on last draw() call
    inline static long int lasttime; // in milliseconds

    inline static int barwidth = 100;
    static const int displaythreads = 8;
    inline static std::atomic<int> currentthreads;  // current working (displayed) threads
    inline static std::mutex currevtmutex;
    bool worker = true ; // Defines a thread to be displayed
    inline static std::vector<uint> totevent;  // number of object to be calculated
    uint threadpos = 0;  // thread position

    inline static std::atomic<int> finishercount;  // number of finished threads
    inline static std::atomic<bool> ongoinganalysis; // switch for status of analysis
};