//
// Created by afrotscher on 8/7/18.
//

#include <libconstant.h>
#include "PID/pid.h"
#include "histogram_cuts.hh"
#include <thread>
#include "progress.h"
#include "TF1.h"
#include <algorithm>

using namespace std;

void PID::innerloop(treereader *tree, std::vector<std::atomic<bool>>
                    &goodevents, std::vector<uint> range) {
    // Step 1: duplicate the data structure
    vector<TH1D> _reactF5;
    vector<vector<TH2D>> _PIDplot;

    for(auto &i: reactF5) _reactF5.emplace_back(TH1D(i));
    for(auto &i: PIDplot){
        vector<TH2D> temp;
        for(auto &j:i) temp.emplace_back(TH2D(j));
        _PIDplot.emplace_back(temp);
    }

    uint threadno = range.at(0)/(range.at(1)-range.at(0));
    progressbar progress(range.at(1)-range.at(0), threadno);

    vector<vector<double>> valinc;     // Store temporary beam values

    for(uint eventcounter=range.at(0); eventcounter<range.at(1); eventcounter++){
        if(goodevents.at(eventcounter)){
            tree->getevent(eventcounter);
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
                _PIDplot.at(0).at(0).Fill(tree->BigRIPSBeam_aoq[0],
                                     tree->BigRIPSBeam_zet[0]);
                _PIDplot.at(1).at(0).Fill(beamaoqcorr, tree->BigRIPSBeam_zet[0]);
            }

            valinc.push_back({tree->BigRIPSBeam_aoq[2],tree->BigRIPSBeam_aoq[3],
                              tree->BigRIPSBeam_aoq[4]});
            valinc.push_back({tree->BigRIPSBeam_zet[2],tree->BigRIPSBeam_zet[3],
                              tree->BigRIPSBeam_zet[4]});
            if(closeness(valinc.at(2)) && closeness(valinc.at(3))){
                _PIDplot.at(0).at(1).Fill(tree->BigRIPSBeam_aoq[4],
                                     tree->BigRIPSBeam_zet[4]);
                _PIDplot.at(1).at(1).Fill(beamaoqcorr2, tree->BigRIPSBeam_zet[4]);
            }
            valinc.clear();

            // We now fill the cut data (cut by ellipsoid)
            if((pow(1./incval.at(2)*(beamaoqcorr-incval.at(0)),2) +
                pow(1/incval.at(3)*(tree->BigRIPSBeam_zet[0]-incval.at(1)),2))<1){
                _PIDplot.at(0).at(2).Fill(tree->BigRIPSBeam_aoq[0],
                                     tree->BigRIPSBeam_zet[0]);
                _PIDplot.at(1).at(2).Fill(beamaoqcorr, tree->BigRIPSBeam_zet[0]);
                _PIDplot.at(0).at(3).Fill(tree->BigRIPSBeam_aoq[4],
                                     tree->BigRIPSBeam_zet[4]);
                _PIDplot.at(1).at(3).Fill(beamaoqcorr2, tree->BigRIPSBeam_zet[4]);
                reactionpid1++;

                // Fill F7 value with PID
                _reactF5.at(0).Fill(tree->F5X);

                // Second ellipsoid for cross section
                if((pow(1./targetval.at(2)*(beamaoqcorr2-targetval.at(0)),2) +
                    pow(1./targetval.at(3)*(tree->BigRIPSBeam_zet[4]-targetval.at(1)),2))<1) {
                    reactionpid2++;
                    // Investigate F7 position of (p,2p) ions (off center effects)
                    _reactF5.at(1).Fill(tree->F5X);

                    //Check double cut particles
                    _PIDplot.at(0).at(4).Fill(tree->BigRIPSBeam_aoq[4],
                                              tree->BigRIPSBeam_zet[4]);
                    _PIDplot.at(1).at(4).Fill(beamaoqcorr2, tree->BigRIPSBeam_zet[4]);
                }
            }
        }
        progress.increaseevent();
    }

    // Step 3: rejoining data structure
    unitemutex.lock();
    for(uint i=0;i<reactF5.size();i++) reactF5.at(i).Add(new TH1D(_reactF5.at(i)));
    for(uint i=0;i<PIDplot.size();i++){
        for(uint j=0;j<PIDplot.at(0).size();j++){
            PIDplot.at(i).at(j).Add(new TH2D(_PIDplot.at(i).at(j)));
        }
    }
    unitemutex.unlock();

    progress.reset();
}

void PID::analyse(std::vector <std::string> input, TFile *output) {
    // Set Parameters according to reaction specified in string
    PID::reactionparameters();

    // Generate folder names and check for previous analysises
    vector<string> folders{"PID/"+reaction+"/Uncorrected",
                           "PID/"+reaction+"/Corrected",
                           "PID/"+reaction+"/investigate"};
    for (auto &i :folders){
        // Do not analyse the same thing twice
        if (output->GetDirectory(i.c_str())) return;
        output->mkdir(i.c_str());
    }

    printf("Making PID with %i threads now...\n", threads);

    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto &h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    vector<string> keys{"BigRIPSBeam.aoq", "BigRIPSBeam.zet", "F3X", "F3A",
                         "F5X", "F5A", "F9X", "F9A", "F11X", "F11A",
                         "BigRIPSIC.fCalMeVSqSum"};
    for(auto &i:tree) i->setloopkeys(keys);

    //Constructing all the histograms
    PID::histogramsetup();

    //Making threads

    vector<thread> th;
    for(uint i=0; i<threads;i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&PID::innerloop, this,
                               tree.at(i),ref(goodevents), ranges));
    }

    for(auto &t : th) t.detach();

    progressbar finishcondition;
    while(finishcondition.ongoing()) finishcondition.draw();

    // Calculate the transmission and the crosssection
    PID::offctrans();
    PID::crosssection();

    for(uint i=0; i<folders.size(); i++) {
        output->cd(folders.at(i).c_str());
        if (i < 2) for (auto &elem: PIDplot.at(i)) elem.Write();
        else {
            for (auto &elem: reactF5) elem.Write();
            fitplot.Write();
        }
        output->cd("");
    }
}

void PID::offctrans() {
    // Calculating off-center transmission
    //int threshold = 200; // 100cts for reacted beam
    /*double dividend = 0; // dividend/divisor
    double divisor = 0;
    int binlow  = reactF5.at(0).FindBin(acceptancerange.at(0));
    int binhigh = reactF5.at(0).FindBin(acceptancerange.at(1));
    //printf("Chisqure minimization from bin %i to bin %i\n", binlow, binhigh);
    vector<double> chisq ={0,0,0};

    // Calculating scaling factor alpha
    for(int i=binlow; i < binhigh; i++){
        if(reactF5.at(1).GetBinContent(i)){
            dividend += reactF5.at(0).GetBinContent(i);
            divisor  += pow(reactF5.at(0).GetBinContent(i),2)/
                        (double)reactF5.at(1).GetBinContent(i);
        }
    }*/

    //double alpha =dividend/divisor; // Scaling factor F5(beam) F5(reacted)

    // Calculating corresponding chi-square, and test minimum
    /*for(int i=binlow; i<binhigh; i++)
        if(reactF5.at(1).GetBinContent(i)) {
            double e = 0.05; // testvalue
            chisq.at(0) += pow(alpha * reactF5.at(0).GetBinContent(i) -
                               reactF5.at(1).GetBinContent(i), 2) /
                           (double) reactF5.at(1).GetBinContent(i);
            chisq.at(1) += pow((1-e)*alpha * reactF5.at(0).GetBinContent(i) -
                               reactF5.at(1).GetBinContent(i), 2) /
                           (double) reactF5.at(1).GetBinContent(i);
            chisq.at(2) += pow((1+e)*alpha * reactF5.at(0).GetBinContent(i) -
                               reactF5.at(1).GetBinContent(i), 2) /
                           (double) reactF5.at(1).GetBinContent(i);
        }*/

    // Pushing out the new scaled beam values
    //for(auto &i: chisq) i /= binhigh - binlow -1;
    /*reactF5.push_back(reactF5.at(0));
    reactF5.back().Scale(alpha);
    reactF5.back().SetTitle("F5 scaled beam profile");*/

    // Calculating raw beam ratio:
    reactF5.push_back(reactF5.at(1));
    for (auto &i:reactF5) i.Sumw2();
    reactF5.back().Divide(new TH1D(reactF5.at(0)));
    reactF5.back().SetTitle("F5 beam profile ratio");
    reactF5.back().SetName("F5ratio");

    // New fit method
    const int minrange =6;
    const uint startbin = (uint)(reactF5.back().FindFirstBinAbove(0));
    const uint stopbin = (uint)(reactF5.back().FindLastBinAbove(0));
    for(uint i=startbin; i<=stopbin; i++){ // Startbin loop
        vector<TF1*> temp;
        for(uint j=i+minrange; j<=stopbin; j++){ // Length loop
            auto corrfit = new TF1(Form("Fit %i-%i", i,j),constfit,
                    reactF5.back().GetBinCenter(i),
                    reactF5.back().GetBinCenter(j),1);
            reactF5.back().Fit(corrfit, "RQE");
            if(corrfit->GetNDF() >= (minrange - 1))
                fitplot.SetBinContent(i,j-i,corrfit->GetChisquare()/
                                            corrfit->GetNDF());
            temp.push_back(corrfit);
        }
        fitstyle.push_back(temp);
    }
    // Get best Fit
    for(uint i= (uint)(fitplot.GetNbinsY()-1); i>minrange; i--){
        for(uint j=startbin; j<fitplot.GetNbinsX(); j++){
            if((fitplot.GetBinContent(j,i) > 0) &&
               (fitplot.GetBinContent(j,i) < maxchisq)){

                reactF5.push_back(reactF5.at(0));
                reactF5.back().Scale(fitstyle.at(j-startbin).at(i-minrange)->GetParameter(0));
                reactF5.back().SetTitle("F5 scaled beam profile");

                // Write off center trans, it cannot exceed 1 however
                offcentertransmission = min(reactF5.at(1).Integral()/reactF5.back().Integral(),1.);

                offcentertransmissionerror =
                        fitstyle.at(j-startbin).at(i-minrange)->GetParError(0)/
                        fitstyle.at(j-startbin).at(i-minrange)->GetParameter(0)*
                        offcentertransmission;

                cout << fitstyle.at(j-startbin).at(i-minrange)->GetName()
                     << " Chisq: " << fitplot.GetBinContent(j,i)
                     << " Offcentertransmission: " << offcentertransmission
                     << " +- " << offcentertransmissionerror << endl;

                return;
            }
        }
    }

    // Off-center cut works only for sufficient statistics:
    /*if((goodevents.size() == runinfo::emptysize) ||
            (goodevents.size() == runinfo::transsize)){
        printf("Not doing off-center transmission. No physics run.\n");
        return;
    }
    else if(reactF5.at(1).Integral() > threshold){
        printf("\n Scaling Factor: %f, with chisq/ndof %f, transmission: %f. "
            "Scaling factor from Area %f !\n",
            alpha, chisq.at(0),reactF5.at(1).Integral()/reactF5.at(2).Integral(),
            reactF5.at(1).Integral(binlow,binhigh)/
            reactF5.at(0).Integral(binlow,binhigh));
        offcentertransmission = reactF5.at(1).Integral()/reactF5.at(2).Integral();
        return;
    }
    else{
        printf("Too low statistics for offcenter effects (%.1f)\n",
                reactF5.at(1).Integral());
        return;
    }*/
}

void PID::crosssection() {
    // Calculate the final crosssection for this run
    const double numberdensity = 0.433; // atoms/cm2 for 10cm LH2
    const double numberdensityerror = 0.00924; // relative error
    const double tottransmission = 0.7754; // empty target transmission for centered beam
    const double tottransmissionerror = 0.01874; // associated relative error

    const double crosssection = 1./numberdensity*reactionpid2/
                          reactionpid1/tottransmission/offcentertransmission;
    const double cserror = crosssection*
                     pow(1./reactionpid1 + 2./reactionpid2+
                         pow(offcentertransmissionerror/offcentertransmission,2)+
                         pow(tottransmissionerror,2)+
                         pow(numberdensityerror,2), 0.5);  // Error on Transm.

    printf("Inclusive N=%.2f Z=%.2f to N=%.2f Z=%.2f sigma is: %f +- %f b\n",
           incval.at(0)*incval.at(1), incval.at(1), targetval.at(0)*targetval.at(1),
           targetval.at(1), crosssection, cserror);

    printf("Raw: in %u out %i ratio %.3f %%\n",reactionpid1.load(),
           reactionpid2.load(),
           100.*reactionpid2/reactionpid1.load());
}

void PID::reactionparameters() {
    // Setup cut values
    switch(goodevents.size()){
        case runinfo::transsize:{
            incval = nancytrans::incval;
            targetval = nancytrans::targetval;
            reaction = "transmission";
            break;
        }
        case runinfo::emptysize:{
            incval = nancyempty::incval;
            targetval = nancyempty::targetval;
            reaction = "empty";
            break;
        }
        default:{
            if(reaction == "111NbPPN"){
                incval = nancy::incval111Nb;
                targetval = nancy::targetval110Nb;
                binning = 50;
                acceptancerange = vector<double>{0,70};
            }
            else if(reaction == "111NbPP2N"){
                incval = nancy::incval111Nb;
                targetval = nancy::targetval109Nb;
                binning = 50;
                acceptancerange = vector<double>{15,70};
            }
            else if(reaction == "111NbP2P"){
                incval    = nancy::incval111Nb;
                targetval = nancy::targetval110Zr;
                binning   = 40;
                acceptancerange = vector<double>{-50,0};
            }
            else if(reaction == "110NbPPN"){
                incval = nancy::incval110Nb;
                targetval = nancy::targetval109Nb;
                binning = 100;
                acceptancerange = vector<double>{-10,70};
            }
            else if(reaction == "110NbP2P"){
                incval = nancy::incval110Nb;
                targetval = nancy::targetval109Zr;
                binning = 40;
                acceptancerange = vector<double>{-15,40};
            }
            else __throw_invalid_argument("Invalid reaction !\n");

            break;
        }
    }
}

void PID::histogramsetup() {
    // Nasty histogram setup routine
    const vector<uint> y{200,36,46}; // y boundaries
    const uint xbin = 400; // number of x-bins
    const vector<double> x{2.55,2.85};

    vector<TH2D> temp1, temp2;
    temp1.emplace_back(TH2D("pidinc", "PID Incoming F3-F7",  xbin,x[0],x[1], y[0],y[1],y[2]));
    temp1.emplace_back(TH2D("pidout", "PID Outgoing F8-F11", xbin,x[0],x[1], y[0],y[1],y[2]));
    temp1.emplace_back(TH2D("pidincut","PID Inc cut  F3-F7",xbin,x[0],x[1], y[0],y[1],y[2]));
    temp1.emplace_back(TH2D("pidincutout", "PID Out w/ in cut", xbin,x[0],x[1],y[0],y[1],y[2]));
    temp1.emplace_back(TH2D("pidincutoutcut", "PID Out cut w/ in cut", xbin,x[0],x[1],y[0],y[1],y[2]));

    temp2.emplace_back(TH2D("pidinccorr", "PID Incoming F3-F7",  xbin,x[0],x[1], y[0],y[1],y[2]));
    temp2.emplace_back(TH2D("pidoutcorr", "PID Outgoing F8-F11", xbin,x[0],x[1], y[0],y[1],y[2]));
    temp2.emplace_back(TH2D("pidincutcorr","PID Inc cut  F3-F7",xbin,x[0],x[1],y[0],y[1],y[2]));
    temp2.emplace_back(TH2D("pidincutoutcorr", "PID Out w/ in cut", xbin,x[0],x[1],y[0],y[1],y[2]));
    temp2.emplace_back(TH2D("pidincutoutcutcorr", "PID Out cut w/ in cut", xbin,x[0],x[1],y[0],y[1],y[2]));
    PIDplot.emplace_back(temp1);
    PIDplot.emplace_back(temp2);

    reactF5.emplace_back(TH1D("F5beam", "F5 beam profile", binning,-100,100));
    reactF5.emplace_back(TH1D("F5react", "F5-position of reacted particles", binning,-100,100));

    fitplot = TH2D("chisqfit","reduced #chi^{2}-fitrange", binning-1,1,binning, binning-1,1,binning);
    fitplot.GetXaxis()->SetTitle("Starting Bin");
    fitplot.GetYaxis()->SetTitle("Number of Bins");
    fitplot.SetOption("colz");
    fitplot.SetMaximum(maxchisq);


    for (auto &elem: reactF5){
        elem.GetXaxis()->SetTitle("x [mm]");
        elem.GetYaxis()->SetTitle("N");
    }

    for(auto &elem: PIDplot){
        for(auto &uelem: elem){
            uelem.SetOption("colz");
            uelem.GetXaxis()->SetTitle("A/Q");
            uelem.GetYaxis()->SetTitle("Z");
        }
    }
}