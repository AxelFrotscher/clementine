#include <cutclasses/chargestatecut.h>
#include "MakeAllTree_78Ni.hh"
#include "cutclasses/triggercut.h"
#include "cutclasses/ppaccut.h"
#include "cutclasses/plasticcut.h"
#include "corrections/higherorder.h"
#include "cutclasses/ICcut.h"
#include "histograms.hh"
#include "histogram_cuts.hh"
#include "libconstant.h"

using namespace std;

calibpar p1;
mutex consolemutex;
mutex writemutex;

void dalicalib(treereader *tree, TFile *output){
    // This Method aims to calibrate the 187 detectors of DALI
    printf("Now beginning the Calibration of the NaI crystals... \n");

    vector<string> keys{"DALINaI", "DALINaI.fADC", "DALINaI.id"};
    tree->setloopkeys(keys);

    TH2D gammadetectors("dalispectra", "Spectrum of each Gamma Detector",
                        186,0,186,4096,0,4096);
    gammadetectors.SetOption("colz");
    gammadetectors.GetXaxis()->SetTitle("Detector Number");
    gammadetectors.GetYaxis()->SetTitle("ADC Channel");

    int numdet =0;

    // Progress Bar setup
    int currevt=0; // counting variable
    Long64_t totevents = tree->NumEntries();
    const int downscale = 500; // every n-th event

    while(tree->singleloop()){
        numdet = tree->DALINaI_;
        for(int i=0; i<numdet;i++){
            if(tree->DALINaI_fADC[i]) gammadetectors.Fill(tree->DALINaI_id[i],
                tree->DALINaI_fADC[i]);
        }
        currevt++;
        if(!(currevt%downscale)) progressbar(currevt,totevents,0);
    }

    output->mkdir("DALI");
    output->cd("DALI");
    gammadetectors.Write();
    output->cd("");
    printf("\nFinished DALI Calibration.\n");
}

void pidth(treereader *tree, vector<vector<TH2D>> &PID, const vector<double> &cutval,
           const vector<double> &targetval, vector<uint> &reactionval,
           const vector<int> &range, const vector<atomic<bool>> &goodevents, const int id,
           vector<TH1D> &reactf9){
    //Thread worker for PID plot
    //uint eventcounter =0;
    const int downscale = 50000; // every n-th event
    consolemutex.lock();
    cout << "Created thread " << range.at(0) << endl;
    consolemutex.unlock();

    vector<vector<double>> valinc;     // Store temporary beam values
    for(uint eventcounter=range.at(0); eventcounter<range.at(1); eventcounter++){
        tree->getevent(eventcounter);
        if(goodevents.at(eventcounter)){
            double beamaoqcorr = tree->BigRIPSBeam_aoq[0] + p1.F7absF5X0 -
                                 (p1.F7absF5X+tree->F5X*p1.F7linF5X) +
                                 tree->F5A*p1.F7linF5A +
                                 tree->F3X*p1.F7linF3X;
            //cout << "Thread: " << id << " Cor AOQ: " << beamaoqcorr << endl;
            double beamaoqcorr2 = tree->BigRIPSBeam_aoq[4] +  p1.F11absF9X0-
                                  (p1.F11absF9X+tree->F9X*p1.F11linF9X) -
                                  tree->F9A*p1.F11linF9A - tree->F11A*p1.F11linF11A;

            //Loop over all elements in the tree
            valinc.push_back({tree->BigRIPSBeam_aoq[0],tree->BigRIPSBeam_aoq[1]});
            valinc.push_back({tree->BigRIPSBeam_zet[0],tree->BigRIPSBeam_zet[1]});
            if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
                // Cut Particles that have variating aoq or zet
                PID.at(0).at(0).Fill(tree->BigRIPSBeam_aoq[0],
                                     tree->BigRIPSBeam_zet[0]);
                PID.at(1).at(0).Fill(beamaoqcorr, tree->BigRIPSBeam_zet[0]);
            }

            valinc.push_back({tree->BigRIPSBeam_aoq[2],tree->BigRIPSBeam_aoq[3],
                              tree->BigRIPSBeam_aoq[4]});
            valinc.push_back({tree->BigRIPSBeam_zet[2],tree->BigRIPSBeam_zet[3],
                              tree->BigRIPSBeam_zet[4]});
            if(closeness(valinc.at(2)) && closeness(valinc.at(3))){
                PID.at(0).at(1).Fill(tree->BigRIPSBeam_aoq[4],
                                     tree->BigRIPSBeam_zet[4]);
                PID.at(1).at(1).Fill(beamaoqcorr2, tree->BigRIPSBeam_zet[4]);
            }
            valinc.clear();

            // We now fill the cut data (cut by ellipsoid)
            if((pow(1./cutval.at(2)*(beamaoqcorr-cutval.at(0)),2) +
                pow(1/cutval.at(3)*(tree->BigRIPSBeam_zet[0]-cutval.at(1)),2))<1){
                PID.at(0).at(2).Fill(tree->BigRIPSBeam_aoq[0],
                                     tree->BigRIPSBeam_zet[0]);
                PID.at(1).at(2).Fill(beamaoqcorr, tree->BigRIPSBeam_zet[0]);
                PID.at(0).at(3).Fill(tree->BigRIPSBeam_aoq[4],
                                     tree->BigRIPSBeam_zet[4]);
                PID.at(1).at(3).Fill(beamaoqcorr2, tree->BigRIPSBeam_zet[4]);
                reactionval.at(0) = reactionval.at(0)+1;

                // Fill F7 value with PID
                reactf9.at(0).Fill(tree->F5X);

                // Second ellipsoid for cross section
                if((pow(1./targetval.at(2)*(beamaoqcorr2-targetval.at(0)),2) +
                    pow(1./targetval.at(3)*(tree->BigRIPSBeam_zet[4]-targetval.at(1)),2))<1) {
                    reactionval.at(1) = reactionval.at(1) + 1;
                    // Investigate F7 position of (p,2p) ions (off center effects)
                    reactf9.at(1).Fill(tree->F5X);
                }
            }
        }
        if(!(eventcounter%downscale)){
            consolemutex.lock();
            progressbar(eventcounter-range.at(0),range.at(1)-range.at(0),id);
            consolemutex.unlock();
        }
    }
}

void makepid(const vector<string> input, TFile *output,
             const vector<atomic<bool>> &goodevents){
    // 5 beams, 2incoming, 3outgoing
    printf("Making PID now...\n");
    // Progress Bar setup
    const uint totevents = goodevents.size();
    const int threadno= 20;

    vector<TChain*> chain;
    for(int i=0; i<threadno; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    vector <string> keys{"BigRIPSBeam.aoq", "BigRIPSBeam.zet", "F3X", "F3A",
                         "F5X", "F5A", "F9X", "F9A", "F11X", "F11A",
                         "BigRIPSIC.fCalMeVSqSum"};
    for(auto &i:tree) i->setloopkeys(keys);

    vector<vector<TH2D>> PID {{
        TH2D("pidinc", "PID Incoming F3-F7",  300,2.45,2.9, 200,30,50),
        TH2D("pidout", "PID Outgoing F8-F11", 300,2.45,2.9, 200,30,50),
        TH2D("pidincut","PID Inc cut  F8-F11",300,2.45,2.9,200,30,50),
        TH2D("pidincutout", "PID Out w/ in cut", 300,2.45,2.9,200,30,50)},
       {TH2D("pidinccorr", "PID Incoming F3-F7",  300,2.45,2.9, 200,30,50),
        TH2D("pidoutcorr", "PID Outgoing F8-F11", 300,2.45,2.9, 200,30,50),
        TH2D("pidincutcorr","PID Inc cut  F8-F11",300,2.45,2.9,200,30,50),
        TH2D("pidincutoutcorr", "PID Out w/ in cut", 300,2.45,2.9,200,30,50)}
    };

    vector<TH1D> reactF5 {
        TH1D("F5beam", "F5 beam profile", 1100,-200,2000),
        TH1D("F5react", "F5-position of reacted particles", 1100,-200,2000)
    };
    for (auto &elem: reactF5){
        elem.GetXaxis()->SetTitle("x [mm]");
        elem.GetYaxis()->SetTitle("N");
    }

    for(auto &elem: PID){
        for(auto &uelem: elem){
            uelem.SetOption("colz");
            uelem.GetXaxis()->SetTitle("A/Q");
            uelem.GetYaxis()->SetTitle("Z");
        }
    }

    vector<double> incval; // cut on incoming particles (F7)
    vector<double> targetval; // second cut to detected particles (F11)

    if(runinfo::transsize == goodevents.size()){
        incval = nancytrans::incval;
        targetval = nancytrans::targetval;
    }
    else if(runinfo::emptysize == goodevents.size()){
        incval = nancyempty::incval;
        targetval = nancyempty::targetval;
    }
    else{
        incval = nancy::incval;
        targetval = nancy::targetval;
    }

    //Setup Crossection trigger:
    vector<uint> reactioncounter{0,0};

    // Multiply the file structure to have many threads
    vector<vector<vector<TH2D>>> PIDthread;

    // Set up data ranges
    vector<vector<int>> range;

    // Setup crosssection vector:
    vector<vector<uint>> reactionthread;
    vector<vector<TH1D>> reactF9;

    for(int i=0; i<threadno; i++){
        PIDthread.push_back(PID);
        reactionthread.push_back(reactioncounter);
        range.push_back({i*totevents/threadno, (i+1)*totevents/threadno-1});
        reactF9.push_back(reactF5);
    }
    range.back().back() = totevents;

    cout << "Created all thread structure for Threads: " << threadno << endl;
    vector<thread> th;
    for(int i=0; i<threadno;i++){
        th.emplace_back(
            thread(pidth, tree.at(i), ref(PIDthread.at(i)), ref(incval),
                   ref(targetval), ref(reactionthread.at(i)), ref(range.at(i)),
                   ref(goodevents), i, ref(reactF9.at(i))));
    }

    for(auto &t : th) t.join();

    // Rejoin the data structure
    for(auto &i: reactionthread){
        reactioncounter.at(1) = reactioncounter.at(1) + i.at(1);
        reactioncounter.at(0) = reactioncounter.at(0) + i.at(0);
    }

    for(uint i=0; i<PID.size();i++){
        for(uint j=0; j<PID.at(0).size(); j++){
            for(auto &hist:PIDthread)
                PID.at(i).at(j).Add(new TH2D(hist.at(i).at(j)));
        }
    }
    for(auto &i: reactF9){
        reactF5.at(0).Add(new TH1D(i.at(0)));
        reactF5.at(1).Add(new TH1D(i.at(1)));
    }

    // Using custom function to scale and get output
    const vector<double> acceptancerange{0.,75.}; // mm
    int dividend = 0; // dividend/divisor
    double divisor = 0;
    int binlow  = reactF5.at(0).FindBin(acceptancerange.at(0));
    int binhigh = reactF5.at(0).FindBin(acceptancerange.at(1));
    double chisq =0;

    // Calculating scaling factor alpha
    for(int i=binlow; i < binhigh; i++){
        dividend += reactF5.at(0).GetBinContent(i);
        divisor  += pow(reactF5.at(0).GetBinContent(i),2)/
                (double)reactF5.at(1).GetBinContent(i);
    }
    double alpha =dividend/divisor;

    // Calculating corresponding chi-square
    for(int i=binlow; i<binhigh; i++){
        chisq += pow(alpha*reactF5.at(0).GetBinContent(i)-
                 reactF5.at(1).GetBinContent(i),2)/
                (double)reactF5.at(1).GetBinContent(i);
    }
    chisq /= binhigh - binlow -1;
    reactF5.push_back(reactF5.at(0));
    reactF5.back().Scale(alpha);
    reactF5.back().SetTitle("F5 scaled beam profile");

    double curpart = reactF5.at(1).Integral()/reactF5.at(2).Integral();
    printf("\n Scaling Factor: %f, with chisq/ndof %f, transmission: %f !\n",
           alpha, chisq,curpart);

    printf("Finished making PIDs!\n");
    vector<string> folders{"PID/Uncorrected","PID/Corrected","PID/investigate"};
    for (auto &i :folders) output->mkdir(i.c_str());

    for(uint i=0; i<folders.size(); i++){
        output->cd(folders.at(i).c_str());
        if(i<2) for(auto &elem: PID.at(i)) elem.Write();
        else for(auto &elem: reactF5) elem.Write();
    }
    output->cd("");

    double crosssection = 1./0.433*reactioncounter.at(1)/reactioncounter.at(0)/0.7754;
    double cserror = crosssection*pow(1./reactioncounter.at(0) +
                                      1./reactioncounter.at(1),0.5);
    printf("Inclusive 111Nb(p,2p)110Zr sigma is: %f +- %f b\n", crosssection,
           cserror);

    printf("Raw: in %i out %i ratio %f %%\n",reactioncounter.at(0),reactioncounter.at(1),
           100.*reactioncounter.at(1)/reactioncounter.at(0));
}

void makehistograms(const vector<string> input) {
    const int threadno = 10;
    cout << "Making " << threadno << " new Threads..." << endl;

    vector<TChain*> chain;
    for(int i=0; i<threadno; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    // Generating an outputfile that matches names with the input file
    string output = "build/output/" + input.at(0).substr(34,9) + "hist" +
                    to_string(input.size()) + ".root";

    auto outputfile = new TFile(output.c_str(), "RECREATE");
    if(!outputfile->IsOpen()) __throw_invalid_argument("Output file not valid");

    cout << "Beginning reconstruction of " << chain.at(0)->GetEntries()
         << " Elements." << endl;
    // Store events that cannot be used
    vector<atomic<bool>> goodevents((uint)chain.at(0)->GetEntries());
    for(auto &i:goodevents) i.exchange(true);

    // Determine run type:
    if(runinfo::transsize == chain.at(0)->GetEntries()){
        printf("!!! Analysing an transmission run !!!\n");
    }
    if(runinfo::emptysize == chain.at(0)->GetEntries()){
        printf("!!! Analysing an empty target run !!!\n");
    }

    triggercut(tree, goodevents);
    iccut(tree,goodevents,outputfile);
    plasticcut(tree, goodevents, outputfile);
    ppaccut(tree, goodevents, outputfile);
    ccsc(tree,goodevents,outputfile);
    higherorder(tree, goodevents, outputfile);

    vector<thread> th;
    //th.emplace_back(thread(plastics, tree.at(0), outputfile, ref(goodevents)));
    //th.emplace_back(thread(chargestatecut, tree.at(1), outputfile, ref(goodevents)));
    //th.emplace_back(thread(ionisationchamber, tree.at(2), outputfile, ref(goodevents)));
    //th.emplace_back(thread(ppacs, tree.at(4), outputfile, ref(goodevents)));
    //th.emplace_back(thread(highordercorrection, tree.at(3), outputfile, ref(goodevents)));
    th.emplace_back(thread(targetcut,tree.at(5),outputfile, ref(goodevents)));
    //th.emplace_back(thread(triggercut, tree.at(6),outputfile, ref(goodevents)));

    for(auto &i: th) i.join();

    /*if(options.at(6)) chargestatecut(alt2dtree, outputfile, goodevents);
    if (options.at(0)) plastics(alt2dtree, outputfile, goodevents);
    if (options.at(2)) ionisationchamber(alt2dtree, outputfile, goodevents);
    if (options.at(1)) ppacs(alt2dtree, outputfile, goodevents);
    if (options.at(3)) highordercorrection(alt2dtree, outputfile, goodevents);*/

    // Get Z vs. A/Q
    makepid(input, outputfile, goodevents);
    printf("Made PID histograms in %s\n", output.c_str());
    //Get ADC Spectra for DALI
    //dalicalib(alt4dtree, outputfile);

    cout << "Run has " <<100.* accumulate(goodevents.begin(),goodevents.end(),0)
                          /goodevents.size() << " % good Elements" << endl;

    outputfile->Close();
}
