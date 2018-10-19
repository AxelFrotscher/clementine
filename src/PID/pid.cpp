//
// Created by afrotscher on 8/7/18.
//

#include <libconstant.h>
#include "PID/pid.h"
#include "histogram_cuts.hh"
#include "txtwriter.h"
#include <thread>
#include "progress.h"
#include "TF1.h"
#include <algorithm>
#include "zdssetting.h"
#include <sstream>

using namespace std;

void PID::innerloop(treereader *tree, std::vector<std::atomic<bool>>
                    &goodevents, std::vector<uint> range) {
    // Step 1: duplicate the data structure
    vector<TH1D> _reactF5;
    vector<TH2D> _reactF7PPAC;
    vector<vector<TH2D>> _PIDplot;

    for(auto &i: reactF5) _reactF5.emplace_back(TH1D(i));
    for(auto &i: reactF7PPAC) _reactF7PPAC.emplace_back(TH2D(i));
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
                                 (p1.F7absF5X+tree->F5X*p1.F7linF5X) -
                                 tree->F5A*p1.F7linF5A -
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

                    // Fill Information to the beam profile/shape/energy
                    // 17 == F7PPAC-2B
                    _reactF7PPAC.at(0).Fill(tree->BigRIPSPPAC_fX[17],
                                            tree->BigRIPSPPAC_fY[17]);
                    _reactF7PPAC.at(1).Fill(
                        1E3*atan2((tree->BigRIPSPPAC_fX[17]-tree->BigRIPSPPAC_fX[15]), 945.),
                        1E3*atan2((tree->BigRIPSPPAC_fY[17]-tree->BigRIPSPPAC_fY[15]), 945.));
                    _reactF7PPAC.at(2).Fill(tree->BigRIPSIC_fCalMeVSqSum[0],
                                            tree->BigRIPSIC_fCalMeVSqSum[1]);
                }
            }
        }
        progress.increaseevent();
    }

    // Step 3: rejoining data structure
    unitemutex.lock();
    for(uint i=0;i<reactF5.size();i++) reactF5.at(i).Add(new TH1D(_reactF5.at(i)));
    for(uint i=0;i<reactF7PPAC.size();i++) reactF7PPAC.at(i).Add(new TH2D(_reactF7PPAC.at(i)));
    for(uint i=0;i<PIDplot.size();i++){
        for(uint j=0;j<PIDplot.at(0).size();j++){
            PIDplot.at(i).at(j).Add(new TH2D(_PIDplot.at(i).at(j)));
        }
    }
    unitemutex.unlock();

    progress.reset();
}

void PID::analyse(const std::vector <std::string> &input, TFile *output) {
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
                         "BigRIPSIC.fCalMeVSqSum", "BigRIPSPPAC.fX",
                         "BigRIPSPPAC.fY"};
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
            for(auto &elem: reactF5) elem.Write();
            fitplot.Write();
            for(auto &elem: reactF7PPAC) elem.Write();
            if(bestfit) bestfit->Write();
        }
        output->cd("");
    }
}

void PID::offctrans() {
    // Calculating raw beam ratio:
    reactF5.push_back(reactF5.at(1));
    for (auto &i:reactF5) i.Sumw2();
    reactF5.back().Divide(new TH1D(reactF5.at(0)));
    reactF5.back().SetName("F5ratio");

    if((goodevents.size() == runinfo::emptysize.at(0)) ||
       (goodevents.size() == runinfo::transsize.at(0))){
        printf("Not doing off-center transmission. No physics run.\n");
        return;
    }

    if(reactF5.at(1).Integral() < 75){
        printf("Not doing off-center transmission. Lack of statistics %.2f cts\n",
               reactF5.at(1).Integral());
        return;
    }

    // Get weighted mean of ratio:
    const int startbin = (uint)(reactF5.back().FindFirstBinAbove(0));
    const int stopbin = (uint)(reactF5.back().FindLastBinAbove(0));
    double sum =0, totalweight =0;
    for(int i=startbin; i<stopbin;i++){
        if((bool)reactF5.back().GetBinContent(i)) {
            sum += reactF5.back().GetBinContent(i) / pow(reactF5.back().GetBinError(i), 2);
            totalweight += pow(reactF5.back().GetBinError(i), -2);
        }
    }

    // New fit method
    const int minrange = 4;
    const double mean = sum/totalweight;
    for(int i=startbin; i<=stopbin; i++){ // Startbin loop
        vector<TF1*> temp;
        for(int j=i+minrange; j<=stopbin; j++){ // Length loop
            auto corrfit = new TF1(Form("Fit %i-%i", i,j),constfit,
                    reactF5.back().GetBinCenter(i),
                    reactF5.back().GetBinCenter(j),1);
            reactF5.back().Fit(corrfit, "NRQ");
            if(corrfit->GetNDF() >= (minrange - 1) && corrfit->GetParameter(0) > 0.95*mean)
                fitplot.SetBinContent(i,j-i,corrfit->GetChisquare()/
                                            corrfit->GetNDF());
            temp.push_back(corrfit);
        }
        fitstyle.push_back(temp);
    }

    // Write out mean to histogram
    reactF5.back().SetTitle(Form("F5 beam profile ratio [mean %.4f (E-3)]",
                            1E3*mean));

    // Get best Fit
    vector<int> backupfit{0,0};
    double minchisq = 10E8; // current minimum chi-sq.

    for(int i= (fitplot.GetNbinsY()-1); i>=minrange; i--){ // switch for effect (current normal)
        for(int j=startbin; j<fitplot.GetNbinsX(); j++){
            //Backup preparation:
            if((fitplot.GetBinContent(j,i) < minchisq) &&
                    (fitplot.GetBinContent(j,i) > 0)){
               minchisq =  fitplot.GetBinContent(j,i);
               backupfit = vector<int>{j,i};
            }

            // Main sequence
            if((fitplot.GetBinContent(j,i) > 0) &&
               (fitplot.GetBinContent(j,i) < maxchisq)){

                reactF5.push_back(reactF5.at(0));
                reactF5.back().Scale(fitstyle.at(j-startbin).at(i-minrange)->GetParameter(0));
                reactF5.back().SetTitle("F5 scaled beam profile");
                bestfit = fitstyle.at(j-startbin).at(i-minrange);

                // Write off center trans, it cannot exceed 1 however
                offcentertransmission = min(reactF5.at(1).Integral()/
                                            reactF5.back().Integral(),1.);

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

    // In case there was no sufficient fit, the lowest chisq is taken...
    printf("Only bad fits available, using lowest chisq. %f\n", minchisq);

    if(!backupfit.at(0)){
        printf("Fuck this shit. Not a single chisq. fit was ok. skipping. Mean %f\n",mean);
        return;
    }

    reactF5.push_back(reactF5.at(0));
    reactF5.back().Scale(fitstyle.at(backupfit[0]-startbin).at(backupfit[1]-minrange)->GetParameter(0));
    reactF5.back().SetTitle("F5 scaled beam profile");

    bestfit = fitstyle.at(backupfit[0]-startbin).at(backupfit[1]-minrange);

    // Write off center trans, it cannot exceed 1 however
    offcentertransmission = min(reactF5.at(1).Integral()/
                                reactF5.back().Integral(),1.);

    offcentertransmissionerror =
            fitstyle.at(backupfit[0]-startbin).at(backupfit[1]-minrange)->GetParError(0)/
            fitstyle.at(backupfit[0]-startbin).at(backupfit[1]-minrange)->GetParameter(0)*
            offcentertransmission;

    cout << fitstyle.at(backupfit[0]-startbin).at(backupfit[1]-minrange)->GetName()
         << " Chisq: " << fitplot.GetBinContent(backupfit[0],backupfit[1])
         << " Offcentertransmission: " << offcentertransmission
         << " \u00b1 " << offcentertransmissionerror << endl;
}

void PID::crosssection() {
    // Calculate the final crosssection for this run
    const double numberdensity = 0.433; // atoms/cm2 for 10cm LH2
    const double numberdensityerror = 0.00924; // relative error
    const double tottransmission = 0.7839; // empty target transmission for centered beam
    const double tottransmissionerror = 0.01865; // associated relative error

    const double crosssection = 1./numberdensity*reactionpid2/
                          reactionpid1/tottransmission/offcentertransmission;
    const double cserror = crosssection*
                     pow(1./reactionpid1 + 2./reactionpid2+
                         pow(offcentertransmissionerror/offcentertransmission,2)+
                         pow(tottransmissionerror,2)+
                         pow(numberdensityerror,2), 0.5);  // Error on Transm.

    // make cross section string:
    setting set;
    stringstream stringout;
    stringout << reaction << " \u03C3: " << setprecision(4) << 1E3*crosssection
        << " \u00b1 " << 1E3*cserror <<  "mb, Offc. Transmission: "
        << offcentertransmission << " \u00b1 " << offcentertransmissionerror
        << " CTS: " << reactF5.at(1).Integral();
    if(set.isemptyortrans())
        stringout << " Ratio: " << 100.*reactionpid2/reactionpid1.load() <<" \u00b1 "
                  << 100.*reactionpid2/reactionpid1.load()*
                     pow(1./reactionpid2+ 1./reactionpid1.load(),0.5) << " %";
    cout  << stringout.str() << endl;

    txtwriter txt;
    txt.addline(stringout.str());

    cout << "Raw: In " << reactionpid1.load() << " out " << reactionpid2.load()
         << " ratio " << 100.*reactionpid2/reactionpid1.load() << "% " << endl;
}

void PID::reactionparameters() {
    // Setup cut values
    setting set;
    if(set.isemptyortrans()){
        incval = set.getPIDincutvalue();
        targetval = set.getPIDoutcutvalue();
        reaction = set.getmodename();

        return;
    }

    binning = 50;
    if(reaction == "111NbPPN"){
        incval = nancy::incval111Nb;
        targetval = nancy::targetval110Nb; }
    else if(reaction == "111NbPP2N"){
        incval = nancy::incval111Nb;
        targetval = nancy::targetval109Nb; }
    else if(reaction == "111NbP2P"){
        incval    = nancy::incval111Nb;
        targetval = nancy::targetval110Zr;
        binning   = 80; }
    else if(reaction == "110NbPPN"){
        incval = nancy::incval110Nb;
        targetval = nancy::targetval109Nb;
        binning = 100; }
    else if(reaction == "110NbP2P"){
        incval = nancy::incval110Nb;
        targetval = nancy::targetval109Zr;
        binning = 40; }
    else if(reaction == "110MoP3P"){
        incval = nancy::incval110Mo;
        targetval = nancy::targetval108Zr; }
    else if(reaction == "111MoP3P"){
        incval = nancy::incval111Mo;
        targetval = nancy::targetval109Zr; }
    else if(reaction == "112MoP3P"){
        incval = nancy::incval112Mo;
        targetval = nancy::targetval110Zr; }
    else if(reaction == "113TcP3P"){
        incval = nancy::incval113Tc;
        targetval = nancy::targetval111Nb; }
    else if(reaction == "112TcP3P"){
        incval = nancy::incval112Tc;
        targetval = nancy::targetval110Nb; }
    else if(reaction == "114TcP3P"){
        incval = nancy::incval114Tc;
        targetval = nancy::targetval112Nb; }
    else if(reaction == "113MoP3P"){
        incval = nancy::incval113Mo;
        targetval = nancy::targetval111Zr; }
    else if(reaction == "90SeP2P"){
        incval = nancy::incval90Se;
        targetval = nancy::targetval89As; }
    else if(reaction == "90SeP3P"){
        incval = nancy::incval90Se;
        targetval = nancy::targetval88Ge; }
    else if(reaction == "89SeP2P"){
        incval = nancy::incval89Se;
        targetval = nancy::targetval88As; }
    else if(reaction == "89SeP3P"){
        incval = nancy::incval89Se;
        targetval = nancy::targetval87Ge; }
    else if(reaction == "88AsP2P"){
        incval = nancy::incval88As;
        targetval = nancy::targetval87Ge; }
    else if(reaction == "89AsP2P"){
        incval = nancy::incval89As;
        targetval = nancy::targetval88Ge; }
    else if(reaction == "89AsP3P"){
        incval = nancy::incval89As;
        targetval = nancy::targetval87Ga; }
    else if(reaction == "93BrP2P"){
        incval = nancy::incval93Br;
        targetval = nancy::targetval92Se; }
    else if(reaction == "93BrP3P"){
        incval = nancy::incval93Br;
        targetval = nancy::targetval91As; }
    else if(reaction == "94BrP2P"){
        incval = nancy::incval94Br;
        targetval = nancy::targetval93Se; }
    else if(reaction == "94BrP3P"){
        incval = nancy::incval94Br;
        targetval = nancy::targetval92As; }
    else if(reaction == "95BrP2P"){
        incval = nancy::incval95Br;
        targetval = nancy::targetval94Se; }
    else if(reaction == "95BrP3P"){
        incval = nancy::incval95Br;
        targetval = nancy::targetval93As; }
    else if(reaction == "94KrP2P"){
        incval = nancy::incval94Kr;
        targetval = nancy::targetval93Br; }
    else if(reaction == "94KrP3P"){
        incval = nancy::incval94Kr;
        targetval = nancy::targetval92Se; }
    else if(reaction == "95KrP2P"){
        incval = nancy::incval95Kr;
        targetval = nancy::targetval94Br; }
    else if(reaction == "95KrP3P"){
        incval = nancy::incval95Kr;
        targetval = nancy::targetval93Se; }
    else if(reaction == "96KrP2P"){
        incval = nancy::incval96Kr;
        targetval = nancy::targetval95Br; }
    else if(reaction == "96KrP3P"){
        incval = nancy::incval96Kr;
        targetval = nancy::targetval94Se; }
    else if(reaction == "97RbP2P"){
        incval = nancy::incval97Rb;
        targetval = nancy::targetval96Kr; }
    else if(reaction == "97RbP3P"){
        incval = nancy::incval97Rb;
        targetval = nancy::targetval95Br; }
    else if(reaction == "99RbP2P"){                            // fourth Setting
        incval = nancy::incval99Rb;
        targetval = nancy::targetval98Kr; }
    else if(reaction == "99RbP3P"){
        incval = nancy::incval99Rb;
        targetval = nancy::targetval97Br; }
    else if(reaction == "100RbP2P"){
        incval = nancy::incval100Rb;
        targetval = nancy::targetval99Kr; }
    else if(reaction == "100RbP3P"){
        incval = nancy::incval100Rb;
        targetval = nancy::targetval98Br; }
    else if(reaction == "100SrP2P"){
        incval = nancy::incval100Sr;
        targetval = nancy::targetval99Rb; }
    else if(reaction == "100SrP3P"){
        incval = nancy::incval100Sr;
        targetval = nancy::targetval98Kr; }
    else if(reaction == "101SrP2P"){
        incval = nancy::incval101Sr;
        targetval = nancy::targetval100Rb; }
    else if(reaction == "101SrP3P"){
        incval = nancy::incval101Sr;
        targetval = nancy::targetval99Kr; }
    else if(reaction == "102SrP2P"){
        incval = nancy::incval102Sr;
        targetval = nancy::targetval101Rb; }
    else if(reaction == "102SrP3P"){
        incval = nancy::incval102Sr;
        targetval = nancy::targetval100Kr; }
    else if(reaction == "102YP2P"){
        incval = nancy::incval102Y;
        targetval = nancy::targetval101Sr; }
    else if(reaction == "102YP3P"){
        incval = nancy::incval102Y;
        targetval = nancy::targetval100Rb; }
    else if(reaction == "103YP2P"){
        incval = nancy::incval103Y;
        targetval = nancy::targetval102Sr; }
    else if(reaction == "103YP3P"){
        incval = nancy::incval103Y;
        targetval = nancy::targetval101Rb; }
    else __throw_invalid_argument(Form("Invalid reaction %s!\n", &reaction[0]));
}

void PID::histogramsetup() {
    // Nasty histogram setup routine
    setting set;
    const vector<uint> y = set.getZrange(); // y boundaries
    const vector<double> x{400, 2.55, 2.85};  // xbins, lower x, upper x

    vector<TH2D> temp1, temp2;
    vector<vector<string>> t1s = {{"pidinc", "PID Incoming F3-F7"},
        {"pidout", "PID Outgoing F8-F11"}, {"pidincut","PID Inc cut  F3-F7"},
        {"pidincutout", "PID Out w/ in cut"},
        {"pidincutoutcut", "PID Out cut w/ in cut"}};
    vector<vector<string>> t2s = {{"pidinccorr", "PID Incoming F3-F7"},
        {"pidoutcorr", "PID Outgoing F8-F11"},{"pidincutcorr","PID Inc cut  F3-F7"},
        {"pidincutoutcorr", "PID Out w/ in cut"},
        {"pidincutoutcutcorr", "PID Out cut w/ in cut"}};

    for(auto &i : t1s) temp1.emplace_back(i.at(0).c_str(),i.at(1).c_str(),
                                          x[0],x[1],x[2],y[0],y[1],y[2]);
    for(auto &i : t2s) temp2.emplace_back(i.at(0).c_str(),i.at(1).c_str(),
                                          x[0],x[1],x[2],y[0],y[1],y[2]);

    PIDplot.emplace_back(temp1);
    PIDplot.emplace_back(temp2);

    reactF5.emplace_back(TH1D("F5beam", "F5 beam profile", binning,-100,100));
    reactF5.emplace_back(TH1D("F5react", "F5-position of reacted particles", binning,-100,100));

    reactF7PPAC.emplace_back(TH2D("F7pos", "PID F7 beamshape", 200,-40,40,200,-40,40));
    reactF7PPAC.emplace_back(TH2D("F7ang", "PID F7 beam angular shape", 100,-100,100,100,-100,100));
    reactF7PPAC.emplace_back(TH2D("F7en", "PID F7 energy distribution", 300,300,900,300,200,800));

    fitplot = TH2D("chisqfit","reduced #chi^{2}-fitrange", binning-1,1,binning, binning-1,1,binning);
    fitplot.GetXaxis()->SetTitle("Starting Bin");
    fitplot.GetYaxis()->SetTitle("Number of Bins");
    fitplot.SetOption("colz");
    fitplot.SetMaximum(2.*maxchisq);

    for (auto &elem: reactF5){
        elem.GetXaxis()->SetTitle("x [mm]");
        elem.GetYaxis()->SetTitle("N");
    }

    for(auto &elem: PIDplot){
        for(auto &uelem: elem){
            uelem.SetOption("colz");
            uelem.GetXaxis()->SetTitle("A/Q");
            uelem.GetYaxis()->SetTitle("Z");
            uelem.SetMinimum(1);
        }
    }

    reactF7PPAC.at(0).GetXaxis()->SetTitle("F7X [mm]");
    reactF7PPAC.at(0).GetYaxis()->SetTitle("F7Y [mm]");
    reactF7PPAC.at(1).GetXaxis()->SetTitle("F7A [mrad]");
    reactF7PPAC.at(1).GetYaxis()->SetTitle("F7B [mrad]");
    reactF7PPAC.at(2).GetXaxis()->SetTitle("E|F7 [MeV]");
    reactF7PPAC.at(2).GetYaxis()->SetTitle("E|F11 [MeV]");

    for(auto &i: reactF7PPAC){
        i.SetOption("colz");
        i.SetMinimum(1);
    }
}