//
// Created by afrotscher on 8/2/18.
//

#include "corrections/higherorder.h"
#include "libconstant.h"

using namespace std;

thread higherorder::innerloop(treereader *tree, std::vector<std::atomic<bool>>
        &goodevents, std::vector<int> range) {
    // Step 1: Cloning histograms
    vector<vector<vector<TH2D>>> _culpritdiag;
    for(auto &i: culpritdiag){
        vector<vector<TH2D>> temp2D;
        for(auto &j:i){
            vector<TH2D> temp1D;
            for(auto &k:j) temp1D.emplace_back(TH2D(k));
            temp2D.emplace_back(temp1D);
        }
        _culpritdiag.emplace_back(temp2D);
    }

    // Step 2: Preparing Variables
    const int downscale = (int)((range.at(1)-range.at(0))/100.);
    int threadno = range.at(0)/(range.at(1)-range.at(0));
    int i = range.at(0); // counting variable

    vector<vector<double>> fillvals(2,vector<double>(9,0)); // Fill dependent variable and

    while(i<range.at(1)){
        if(goodevents.at(i)) { // We absolutely need CCSC cuts for HOC
            tree->getevent(i);
            fillvals.at(0).at(0) = tree->F3X;
            fillvals.at(0).at(1) = tree->F3A;
            fillvals.at(0).at(2) = tree->F5X;
            fillvals.at(0).at(3) = tree->F5A;
            fillvals.at(1).at(0) = tree->F9X;
            fillvals.at(1).at(1) = tree->F9A;
            fillvals.at(1).at(2) = tree->F11X;
            fillvals.at(1).at(3) = tree->F11A;

            // Fill corrected Values
            fillvals.at(0).at(6) = tree->BigRIPSBeam_aoq[beam.at(0)] + p1.F7absF5X0 -
                                   (p1.F7absF5X + tree->F5X * p1.F7linF5X);
            fillvals.at(0).at(7) = fillvals.at(0).at(6) + tree->F5A * p1.F7linF5A;
            fillvals.at(0).at(8) = fillvals.at(0).at(7) + tree->F3X * p1.F7linF3X;

            fillvals.at(1).at(6) = tree->BigRIPSBeam_aoq[beam.at(1)] + p1.F11absF9X0 -
                                   (p1.F11absF9X + tree->F9X * p1.F11linF9X);
            fillvals.at(1).at(7) = fillvals.at(1).at(6) - tree->F9A * p1.F11linF9A;
            fillvals.at(1).at(8) = fillvals.at(1).at(7) - tree->F11A * p1.F11linF11A;

            // Fill pre and post events
            for (uint ii = 0; ii < beam.size(); ii++) {
                if ((pow(1./cutval.at(ii).at(2)*(tree->BigRIPSBeam_aoq[beam.at(ii)] -
                                                 cutval.at(ii).at(0)), 2) +
                     pow(1./cutval.at(ii).at(3)*(tree->BigRIPSBeam_zet[beam.at(ii)] -
                                                 cutval.at(ii).at(1)), 2)) < 1) {
                    // Applying the elliptic cut for 85Ge

                    fillvals.at(ii).at(4) = tree->BigRIPSBeam_beta[beam.at(ii)];
                    fillvals.at(ii).at(5) = tree->BigRIPSBeam_aoq[beam.at(ii)];

                    for (uint k = 0; k < _culpritdiag.at(0).size(); k++) {
                        for (uint j = 0; j <= corrcount; j++) {
                            _culpritdiag.at(ii).at(k).at(j).Fill(
                                fillvals.at(ii).at(5+j), fillvals.at(ii).at(k));
                        }
                    }
                }
            }
        }
        i++;
        if(!((i-range.at(0))%downscale)){
            consolemutex.lock();
            progressbar(i-range.at(0),range.at(1)-range.at(0),threadno);
            consolemutex.unlock();
        }
    }
    // Step 3 rejoining the histogram
    unitemutex.lock();
    printf("Thread %i now rejoining higherorder 2D Histograms...\n", threadno);
    for(uint i=0; i<culpritdiag.size();i++){
        for(uint j=0; j<culpritdiag.at(0).size();j++){
            for(uint k=0; k<culpritdiag.at(0).at(0).size();k++){
                culpritdiag.at(i).at(j).at(k).Add(
                        new TH2D(_culpritdiag.at(i).at(j).at(k)));
            }
        }
    }
    unitemutex.unlock();
}

void higherorder::analyse(std::vector<treereader *> tree, TFile *output) {
    threads = (int)tree.size();
    printf("Higher Order correction with %i threads.\n", threads);

    vector<string> keys{"F5X","F5A", "F3X", "F3A", "F9A", "F9X",
                        "F11A", "F11X", "BigRIPSBeam.aoq",
                        "BigRIPSBeam.beta", "BigRIPSBeam.zet",};
    for(auto &i: tree) i->setloopkeys(keys);

    // Decide which cut to use based on total event number
    if(runinfo::transsize == goodevents.size()){
        cutval = nancytrans::cutval;
        p1 = nancytrans::hoparame;
    }
    else if(runinfo::emptysize == goodevents.size()){
        cutval = nancyempty::cutval;
        p1 = nancyempty::hoparame;
    }
    else {
        cutval = nancy::cutval;
        p1 = nancy::hoparame;
    }

    // Initialize all the diagrams
    for(uint i=0; i<arrname.size();i++){ // Loop F7, F11
        vector<vector<TH2D>> temp2d;
        for(uint j=0; j<arrname.at(0).size();j++){ // Loop F3X,F5X,...
            vector<TH2D> temp1d;
            for(uint k=0;k<=corrcount; k++){ // Loop Correction number
                string arr = arrname.at(i).at(j) + to_string(k); // Generate array name
                string arrn = arrtitle.at(i).at(j) + " Corr: " + to_string(k);
                double ymax = 150;
                double ymin = -100;
                if(j == 4) { // Beta diagrams need smaller bins
                    ymax = 0.7;
                    ymin = 0.6 -0.2*i; // outgoing have lower velocities
                }
                temp1d.emplace_back(TH2D(arr.c_str(),arrn.c_str(),600,2.5,2.8,
                                         400,ymin,ymax));
                temp1d.back().SetOption("colz");
                temp1d.back().GetXaxis()->SetTitle("A/Q");
                temp1d.back().GetYaxis()->SetTitle(arrname.at(i).at(j).c_str());
            }
            temp2d.push_back(temp1d);
        }
        culpritdiag.push_back(temp2d);
    }

    printf("Successfully generated Histograms for higher order...\n");

    vector<thread> th;
    for(uint i=0; i<threads; i++){
        vector<int> ranges = {i*goodevents.size()/threads,
                              (i+1)*goodevents.size()/threads-1};
        th.emplace_back(thread(&higherorder::innerloop, this, tree.at(i),
                               ref(goodevents),ranges));
    }

    for(auto &i: th) i.join();

    // After analysis get correlations represented by a linear fit
    auto corrlinfit = new TF1("Linear Fit", linfit, cutval[0][1]-cutfrac*cutval[0][3],
                              cutval[0][1]+cutfrac*cutval[0][3], 2);
    corrlinfit->SetParNames("absolute", "linear");
    for (auto &elem : culpritdiag) {
        vector<vector<TProfile*>> proftemp2d;
        for(auto &elem2: elem){
            vector<TProfile*> proftemp1d;
            for(uint i_corr=0; i_corr<=corrcount;i_corr++){
                proftemp1d.push_back(elem2.at(i_corr).ProfileY());
                proftemp1d.back()->Fit("Linear Fit","Q");
            }
            proftemp2d.push_back(proftemp1d);
        }
        projections.push_back(proftemp2d);
    }

    // Generate root file structure and then write out al histograms
    for (auto &i_fold: folders) for(auto &j: i_fold) output->mkdir(j.c_str());

    for(uint i=0; i<projections.size();i++){  // Pre, Post
        for(uint k=0; k<projections.at(0).at(0).size(); k++){ // Corr 0,1,2,

            string tempfolder = folders.at(i).at(k) + "/Profiles";
            output->mkdir(tempfolder.c_str());

            for(uint j=0; j<projections.at(0).size();j++){ // F3X, F5X, ...
                output->cd(tempfolder.c_str());
                projections.at(i).at(j).at(k)->Write();

                output->cd(folders.at(i).at(k).c_str());
                culpritdiag.at(i).at(j).at(k).Write();
            }
        }
    }

    output->cd("");
    printf("Finished with higher order corrections!\n");
}