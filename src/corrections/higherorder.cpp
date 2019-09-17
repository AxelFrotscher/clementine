//
// Created by afrotscher on 8/2/18.
//

#include "corrections/higherorder.h"
#include "libconstant.h"
#include "progress.h"
#include "TProfile.h"
#include "TF1.h"
#include <thread>
#include "zdssetting.h"

using std::vector, std::atomic, std::string, std::to_string, std::thread;

void higherorder::innerloop(treereader &tree, const vector<int> &range) {
    // Step 1: Cloning histograms
    decltype(culpritdiag) _culpritdiag(culpritdiag);

    // Step 2: Preparing Variables
    const uint threadno = range.at(0)/(range.at(1)-range.at(0));

    progressbar progress(range.at(1)-range.at(0), threadno);

    vector<vector<double>> fillvals(2,vector<double>(9,0)); // Fill dependent variable and

    for(int i=range.at(0); i<range.at(1); i++){
        if(goodevents.at(i).at(0)) { // We absolutely need CCSC cuts for HOC
            tree.getevent(i);
            fillvals.at(0).at(0) = tree.F3X;
            fillvals.at(0).at(1) = tree.F3A;
            fillvals.at(0).at(2) = tree.F5X;
            fillvals.at(0).at(3) = tree.F5A;
            fillvals.at(1).at(0) = tree.F9X;
            fillvals.at(1).at(1) = tree.F9A;
            fillvals.at(1).at(2) = tree.F11X;
            fillvals.at(1).at(3) = tree.F11A;

            // Fill corrected Values
            fillvals.at(0).at(6) = tree.BigRIPSBeam_aoq[beam.at(0)] + p1.F7absF5X0 -
                                   (p1.F7absF5X + tree.F5X * p1.F7linF5X);
            fillvals.at(0).at(7) = fillvals.at(0).at(6) - tree.F5A * p1.F7linF5A;
            fillvals.at(0).at(8) = fillvals.at(0).at(7) - tree.F3X * p1.F7linF3X;

            fillvals.at(1).at(6) = tree.BigRIPSBeam_aoq[beam.at(1)] + p1.F11absF9X0 -
                                   (p1.F11absF9X + tree.F9X * p1.F11linF9X);
            fillvals.at(1).at(7) = fillvals.at(1).at(6) - tree.F9A * p1.F11linF9A;
            fillvals.at(1).at(8) = fillvals.at(1).at(7) - tree.F11A * p1.F11linF11A;

            // Fill pre and post events (post only when gevts[i][1] ok is)
            for (ulong ii = 0; ii < (beam.size()-1+goodevents.at(i).at(1)); ii++){
                if ((pow(1./cutval.at(ii).at(2)*(tree.BigRIPSBeam_aoq[beam.at(ii)] -
                                                 cutval.at(ii).at(0)), 2) +
                     pow(1./cutval.at(ii).at(3)*(tree.BigRIPSBeam_zet[beam.at(ii)] -
                                                 cutval.at(ii).at(1)), 2)) < 1) {
                    // Applying the elliptic cut for 85Ge

                    fillvals.at(ii).at(4) = tree.BigRIPSBeam_beta[beam.at(ii)];
                    fillvals.at(ii).at(5) = tree.BigRIPSBeam_aoq[beam.at(ii)];

                    for (ulong k = 0; k < _culpritdiag.at(0).size(); k++) {
                        for (int j = 0; j <= corrcount; j++) {
                            _culpritdiag.at(ii).at(k).at(j).Fill(
                                fillvals.at(ii).at(5+j), fillvals.at(ii).at(k));
                        }
                    }
                }
            }
        }
        progress.increaseevent();
    }
    // Step 3 rejoining the histogram
    unitemutex.lock();
    //printf("Thread %i now rejoining higherorder 2D Histograms...\n", threadno);
    for(ulong i=0; i<culpritdiag.size();i++){
        for(ulong j=0; j<culpritdiag.at(0).size();j++){
            for(ulong k=0; k<culpritdiag.at(0).at(0).size();k++){
                culpritdiag.at(i).at(j).at(k).Add(&_culpritdiag.at(i).at(j).at(k));
            }
        }
    }
    unitemutex.unlock();

    progressbar::reset();
}

void higherorder::analyse(const std::vector<std::string> &input, TFile *output) {

    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    printf("Higher Order correction with %i threads.\n", threads);

    vector<string> keys{"F5X","F5A", "F3X", "F3A", "F9A", "F9X",
                        "F11A", "F11X", "BigRIPSBeam.aoq",
                        "BigRIPSBeam.beta", "BigRIPSBeam.zet",};
    for(auto &i: tree) i.setloopkeys(keys);

    // Decide which cut to use based on total event number
    cutval = setting::getHOcutval();
    p1 = setting::getHOparameters();

    // Initialize all the diagrams
    for(ulong i=0; i<arrname.size();i++){ // Loop F7, F11
        culpritdiag.emplace_back();

        for(ulong j=0; j<arrname.at(0).size();j++){ // Loop F3X,F5X,...
            culpritdiag.back().emplace_back();

            for(int k=0;k<=corrcount; k++){ // Loop Correction number
                string arr = arrname.at(i).at(j) + to_string(k); // Generate array name
                string arrn = arrtitle.at(i).at(j) + " Corr: " + to_string(k);
                double ymax = 150;
                double ymin = -100;
                if(j == 4) { // Beta diagrams need smaller bins
                    ymax = 0.7;
                    ymin = 0.6 -0.2*i; // outgoing have lower velocities
                }
                culpritdiag.back().back().emplace_back(
                           arr.c_str(),arrn.c_str(),600,2.5,2.8, 400,ymin,ymax);
                culpritdiag.back().back().back().SetOption("colz");
                culpritdiag.back().back().back().GetXaxis()->SetTitle("A/Q");
                culpritdiag.back().back().back().GetYaxis()->SetTitle(
                                                   arrname.at(i).at(j).c_str());
            }
        }
    }

    printf("Successfully generated Histograms for higher order...\n");

    vector<thread> th;
    for(int i=0; i<threads; i++){
        vector<int> ranges = {(int)(i*goodevents.size()/threads),
                              (int)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&higherorder::innerloop, this,
                               std::ref(tree.at(i)), ref(ranges)));
    }

    for(auto &i: th) i.detach();

    progressbar finishcondition;
    while(progressbar::ongoing()) finishcondition.draw();

    // After analysis get correlations represented by a linear fit
    TF1 corrlinfit("Linear Fit", linfit, cutval[0][1]-cutfrac*cutval[0][3],
                              cutval[0][1]+cutfrac*cutval[0][3], 2);
    corrlinfit.SetParNames("absolute", "linear");
    for (auto &elem : culpritdiag) {
        projections.emplace_back();

        for(auto &elem2: elem){
            projections.back().emplace_back();

            for(int i_corr=0; i_corr<=corrcount;i_corr++){
                projections.back().back().push_back(elem2.at(i_corr).ProfileY());
                projections.back().back().back()->Fit(&corrlinfit, "Q");
            }
        }
    }

    // Generate root file structure and then write out all histograms
    for (auto &i_fold: folders) for(auto &j: i_fold) output->mkdir(j.c_str());

    for(ulong i=0; i<projections.size();i++){  // Pre, Post
        for(ulong k=0; k<projections.at(0).at(0).size(); k++){ // Corr 0,1,2,

            string tempfolder = folders.at(i).at(k) + "/Profiles";
            output->mkdir(tempfolder.c_str());

            for(ulong j=0; j<projections.at(0).size();j++){ // F3X, F5X, ...
                output->cd(tempfolder.c_str());
                projections.at(i).at(j).at(k)->Write();

                output->cd(folders.at(i).at(k).c_str());
                culpritdiag.at(i).at(j).at(k).Write();
            }
        }
    }

    output->cd("");
    printf("Finished with higher order corrections!\n");

    //Clean up tree
    for(auto &i: projections) for(auto &j: i) for(auto &k:j) k->Delete();
}