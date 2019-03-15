//
// Created by afrotscher on 8/7/18.
//

#include "minos.h"
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

using std::vector, std::string, std::atomic, std::thread, std::cout, std::endl,
      std::stringstream, std::to_string, std::__throw_invalid_argument, std::min;

void PID::innerloop(treereader *tree, treereader *minostree,
                    vector<vector<atomic<bool>>> &goodevents, vector<uint> range) {
    // Step 1: duplicate the data structure
    vector<TH1D>         _reactF5;
    vector<TH2D>         _minosresults;
    vector<TH1D>         _minos1dresults;
    vector<vector<TH2D>> _reactPPAC;
    vector<vector<TH2D>> _PIDplot;
    // MINOS with single events
    vector<TH2D> _minossingleevent;

    for(auto &i: reactF5) _reactF5.emplace_back(TH1D(i));

    for(auto &i: reactPPAC){
        vector<TH2D> temp;
        for(auto &j:i) temp.emplace_back(TH2D(j));
        _reactPPAC.emplace_back(temp);
    }

    for(auto &i: PIDplot){
        vector<TH2D> temp;
        for(auto &j:i) temp.emplace_back(TH2D(j));
        _PIDplot.emplace_back(temp);
    }

    for(auto &i: minosresults) _minosresults.emplace_back(TH2D(i));
    for(auto &i: minos1dresults) _minos1dresults.emplace_back(TH1D(i));

    double maxbrho = 10; // Tm, higher than all my values
    if(incval.size() == 8){
        if(     reaction.find("P2P") != string::npos) maxbrho = incval.at(5);
        else if(reaction.find("P3P") != string::npos) maxbrho = incval.at(7);
    }

    uint threadno = range.at(0)/(range.at(1)-range.at(0));
    progressbar progress(range.at(1)-range.at(0), threadno);

    vector<vector<double>> valinc;     // Store temporary beam values

    const vector<int> ppacFpositions{12,17,25,33}; //F5-2B, F7-2B, F9-2B, F11-2B
    const vector<int> ppacangledistance{650,945,700,500}; // in mm

    for(uint eventcounter=range.at(0);eventcounter<range.at(1);eventcounter++){
        if(goodevents.at(eventcounter).at(0)){
            tree->getevent(eventcounter);
            double beamaoqcorr = tree->BigRIPSBeam_aoq[0] + p1.F7absF5X0 -
                                 (p1.F7absF5X+tree->F5X*p1.F7linF5X) -
                                 tree->F5A*p1.F7linF5A -
                                 tree->F3X*p1.F7linF3X;

            double beamaoqcorr2 = tree->BigRIPSBeam_aoq[4] +  p1.F11absF9X0-
                                  (p1.F11absF9X+tree->F9X*p1.F11linF9X) -
                                  tree->F9A*p1.F11linF9A -
                                  tree->F11A*p1.F11linF11A;

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

            // F7-F11 part
            if(closeness(valinc.at(2)) && closeness(valinc.at(3)) &&
               goodevents.at(eventcounter).at(1)){
                _PIDplot.at(0).at(1).Fill(tree->BigRIPSBeam_aoq[4],
                                     tree->BigRIPSBeam_zet[4]);
                _PIDplot.at(1).at(1).Fill(beamaoqcorr2,tree->BigRIPSBeam_zet[4]);
            }
            valinc.clear();

            // We now fill the cut data (cut by ellipsoid and brho)

            if((pow(1./incval.at(2)*(beamaoqcorr-incval.at(0)),2) +
                pow(1/incval.at(3)*(tree->BigRIPSBeam_zet[0]-incval.at(1)),2))
                < 1 ){
                _PIDplot.at(0).at(2).Fill(tree->BigRIPSBeam_aoq[0],
                                     tree->BigRIPSBeam_zet[0]);
                _PIDplot.at(1).at(2).Fill(beamaoqcorr,tree->BigRIPSBeam_zet[0]);

                if(goodevents.at(eventcounter).at(1)) { // F7-F11 part
                    _PIDplot.at(0).at(3).Fill(tree->BigRIPSBeam_aoq[4],
                                              tree->BigRIPSBeam_zet[4]);
                    _PIDplot.at(1).at(3).Fill(beamaoqcorr2,
                                              tree->BigRIPSBeam_zet[4]);
                }
                // cut on energy(Brho) to avoid offcenter Transmission
                if(tree->BigRIPSRIPS_brho[1] < maxbrho) reactionpid1++;

                // Fill F7 value with PID
                _reactF5.at(0).Fill(tree->F5X);

                // Second ellipsoid for cross section F7-F11 part
                if((pow(1./targetval.at(2)*(beamaoqcorr2-targetval.at(0)),2) +
                    pow(1./targetval.at(3)*
                    (tree->BigRIPSBeam_zet[4]-targetval.at(1)),2))<1
                    && goodevents.at(eventcounter).at(1)) {
                    // cut on energy(Brho) to avoid offcenter Transmission
                    if(tree->BigRIPSRIPS_brho[1] < maxbrho) reactionpid2++;
                    // Investigate F7 position of (p,2p) ions (off center effects)
                    _reactF5.at(1).Fill(tree->F5X);

                    //Check double cut particles
                    _PIDplot.at(0).at(4).Fill(tree->BigRIPSBeam_aoq[4],
                                              tree->BigRIPSBeam_zet[4]);
                    _PIDplot.at(1).at(4).Fill(beamaoqcorr2,
                                              tree->BigRIPSBeam_zet[4]);

                    // Fill Information to the beam profile/shape/energy
                    // 17 == F7PPAC-2B, 15 == F7PPAC-1B
                    for(int i=0; i<_reactPPAC.size();i++){
                        _reactPPAC.at(i).at(0).Fill(
                                tree->BigRIPSPPAC_fX[ppacFpositions.at(i)],
                                tree->BigRIPSPPAC_fY[ppacFpositions.at(i)]);
                        _reactPPAC.at(i).at(1).Fill(
                            1E3*atan2(tree->BigRIPSPPAC_fX[ppacFpositions.at(i)]-
                                      tree->BigRIPSPPAC_fX[ppacFpositions.at(i)-2],
                                      ppacangledistance.at(i)),
                            1E3*atan2(tree->BigRIPSPPAC_fY[ppacFpositions.at(i)]-
                                      tree->BigRIPSPPAC_fY[ppacFpositions.at(i)-2],
                                      ppacangledistance.at(i)));
                        _reactPPAC.at(i).at(2).Fill(
                                     tree->BigRIPSPPAC_fX[ppacFpositions.at(i)],
                                     tree->BigRIPSRIPS_brho[1]);
                    }
                    // do the correlation between F5 and F9
                    fitplot.at(1).Fill(tree->F5X, tree->F9X);

                    // make the mighty minos analysis
                    if(reaction.find("P0P") == string::npos){
                        minostree->getevent(eventcounter);
                        minosana analysis(minostree->Trackamount, minostree->Tshaping,
                                          minostree->TimeBinElec, minostree->DelayTrig,
                                          minostree->VDrift, minostree->minostrackxy,
                                          minostree->Minoscalibvalues, minostree->minostime,
                                          minostree->MinosClustX, minostree->MinosClustY,
                                          minostree->MinosClustQ, (int)threadno,
                                          _minossingleevent);
                        TMinosPass minres = analysis.analyze();
                        _minosresults.at(0).Fill(minres.thetaz1, minres.thetaz2);
                        _minosresults.at(1).Fill(minres.trackNbr,minres.trackNbr_final);
                        _minos1dresults.at(0).Fill(minres.z_vertex);
                        _minos1dresults.at(1).Fill(minres.phi_vertex);
                    }
                }
            }
        }
        progress.increaseevent();
    }

    // Step 3: rejoining data structure
    unitemutex.lock();
    for(uint i=0;i<reactF5.size();i++) reactF5.at(i).Add(new TH1D(_reactF5.at(i)));

    for(uint i=0;i<minosresults.size();i++)
        minosresults.at(i).Add(new TH2D(_minosresults.at(i)));

    for(uint i=0; i<minos1dresults.size(); i++)
        minos1dresults.at(i).Add(new TH1D(_minos1dresults.at(i)));

    for(uint i=0;i<reactPPAC.size();i++){
        for(uint j=0; j<reactPPAC.at(0).size(); j++){
            reactPPAC.at(i).at(j).Add(new TH2D(_reactPPAC.at(i).at(j)));
        }
    }

    for(uint i=0;i<PIDplot.size();i++){
        for(uint j=0;j<PIDplot.at(0).size();j++){
            PIDplot.at(i).at(j).Add(new TH2D(_PIDplot.at(i).at(j)));
        }
    }
    minossingleevent.insert(minossingleevent.end(), _minossingleevent.begin(),
                            _minossingleevent.end());

    unitemutex.unlock();

    progress.reset();
}

void PID::analyse(const std::vector <std::string> &input, TFile *output) {
    // Set Parameters according to reaction specified in string
    PID::reactionparameters();

    // Generate folder names and check for previous analysises
    vector<string> folders{
        "PID/"+reaction+"/Uncorrected", "PID/"+reaction+"/Corrected",
        "PID/"+reaction+"/investigate", "PID/"+reaction+"/investigate/F5",
        "PID/"+reaction+"/investigate/F7", "PID/"+reaction+"/investigate/F9",
        "PID/"+reaction+"/investigate/F11", "PID/"+reaction+"/chargestate",
        "PID/"+reaction+"/MINOS", "PID/"+reaction+"/MINOS/events"};

    for (auto &i :folders){
        // Do not analyse the same thing twice
        if (output->GetDirectory(i.c_str())) return;
        output->mkdir(i.c_str());
    }

    printf("Making PID with %i threads now...\n", threads);

    /// Make regular TChain for fast input
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
                        "BigRIPSPPAC.fX", "BigRIPSPPAC.fY", "BigRIPSRIPS.brho"};
    for(auto &i:tree) i->setloopkeys(keys);

    /// Make slow chain for MINOS readout
    vector<TChain*> minoschain;
    for(int i=0; i<threads; i++){
        minoschain.emplace_back(new TChain("tree"));
        for(auto &h: input) minoschain.back()->Add(h.c_str());
    }

    vector<treereader*> minostree;
    for(auto *i:minoschain){
        minostree.emplace_back(new treereader(i));
    }

    vector<string> minoskeys{"VDrift", "MinosClustX", "MinosClustY",
                             "MinosClustQ", "DelayTrig", "Trackamount",
                             "Minoscalibvalues", "TimeBinElec", "Tshaping",
                             "minostrackxy", "minostime"};
    for(auto &i: minostree) i->setloopkeys(minoskeys);

    //Constructing all the histograms
    PID::histogramsetup();

    //Making threads
    vector<thread> th;
    for(uint i=0; i<threads;i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&PID::innerloop, this,
                               tree.at(i), minostree.at(i), ref(goodevents), ranges));
    }

    for(auto &t : th) t.detach();

    progressbar finishcondition;
    while(finishcondition.ongoing()) finishcondition.draw();

    // Calculate the transmission and the crosssection
    PID::offctrans();
    PID::chargestatecut();
    PID::crosssection();
    PID::brhoprojections();

    for(uint i=0; i<2; i++) {
        output->cd(folders.at(i).c_str());
        for (auto &elem: PIDplot.at(i)) elem.Write();
        output->cd("");
    }

    output->cd(folders.at(2).c_str());
    for(auto &i: fitplot) i.Write();

    output->cd(folders.at(3).c_str());
    for(auto &elem: reactF5) elem.Write();
    if(bestfit) bestfit->Write();

    for(int i=0; i<reactPPAC.size();i++){
        output->cd(folders.at(3+i).c_str());
        for(auto &elem: reactPPAC.at(i)) elem.Write();
    }

    for(int i=0; i<brhoprojection.at(0).size(); i++){
        output->cd(folders.at(3+i).c_str());
        brhoprojection.at(0).at(i).Write();
        brhoprojection.at(1).at(i).Write();
    }

    output->cd(folders.at(7).c_str());
    for(auto &elem: chargestate) elem.Write();

    output->cd(folders.at(8).c_str());
    for(auto &elem: minosresults) elem.Write();
    for(auto &elem: minos1dresults) elem.Write();

    output->cd(folders.at(9).c_str());
    for(auto &elem: minossingleevent) elem.Write();
}

void PID::offctrans() {
    // Calculating raw beam ratio:
    reactF5.push_back(reactF5.at(1));
    for (auto &i:reactF5) i.Sumw2();
    reactF5.back().Divide(new TH1D(reactF5.at(0)));
    reactF5.back().SetName("F5ratio");

    setting set;
    if(set.isemptyortrans()){
        printf("Not doing off-center transmission. No physics run.\n");
        return;
    }

    if(reactF5.at(1).Integral() < 75){
        printf("Not doing off-center transmission. Lack of statistics %.2f"
               "cts\n", reactF5.at(1).Integral());
        return;
    }

    // Get weighted mean of ratio:
    const int startbin = (uint)(reactF5.back().FindFirstBinAbove(0));
    const int stopbin = (uint)(reactF5.back().FindLastBinAbove(0));
    double sum =0, totalweight =0;
    for(int i=startbin; i<stopbin;i++){
        if((bool)reactF5.back().GetBinContent(i)) {
            sum += reactF5.back().GetBinContent(i) /
                   pow(reactF5.back().GetBinError(i), 2);
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
                fitplot.at(0).SetBinContent(i,j-i,corrfit->GetChisquare()/
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

    for(int i= (fitplot.at(0).GetNbinsY()-1); i>=minrange; i--){ // switch for effect (current normal)
        for(int j=startbin; j<fitplot.at(0).GetNbinsX(); j++){
            //Backup preparation:
            if((fitplot.at(0).GetBinContent(j,i) < minchisq) &&
                    (fitplot.at(0).GetBinContent(j,i) > 0)){
               minchisq =  fitplot.at(0).GetBinContent(j,i);
               backupfit = vector<int>{j,i};
            }

            // Main sequence
            if((fitplot.at(0).GetBinContent(j,i) > 0) &&
               (fitplot.at(0).GetBinContent(j,i) < maxchisq)){

                reactF5.push_back(reactF5.at(0));
                reactF5.back().Scale(
                        fitstyle.at(j-startbin).at(i-minrange)->GetParameter(0));
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
                     << " Chisq: " << fitplot.at(0).GetBinContent(j,i)
                     << " Offcentertransmission: " << offcentertransmission
                     << " +- " << offcentertransmissionerror << endl;

                return;
            }
        }
    }

    // In case there was no sufficient fit, the lowest chisq is taken...
    printf("Only bad fits available, using lowest chisq. %f\n", minchisq);

    if(!backupfit.at(0)){
        printf("Not a single chisq. fit was ok. skipping. Mean %f\n",mean);
        return;
    }

    reactF5.push_back(reactF5.at(0));
    reactF5.back().Scale(fitstyle.at(backupfit[0]-startbin)
                                 .at(backupfit[1]-minrange)->GetParameter(0));
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
         << " Chisq: " << fitplot.at(0).GetBinContent(backupfit[0],backupfit[1])
         << " Offcentertransmission: " << offcentertransmission
         << " \u00b1 " << offcentertransmissionerror << endl;
}

void PID::crosssection() {
    // Find out reaction type
    double tottransmission = 1;
    if(incval.size() == 8){
        if(     reaction.find("P2P") != string::npos) tottransmission = incval.at(4);
        else if(reaction.find("P3P") != string::npos) tottransmission = incval.at(6);
    }

    // Calculate the final crosssection for this run
    const double numberdensity = 0.433; // atoms/cm2 for 10cm LH2
    const double numberdensityerror = 0.00924; // relative error
    const double tottransmissionerror = 0.05; // conservative relative error

    // Reduce CS contribution by Brho cut
    if(reactF5.at(1).GetEntries())
        chargestatevictims *= reactionpid2.load()/reactF5.at(1).GetEntries();

    double crosssection = std::max(0., 1./numberdensity*(reactionpid2-chargestatevictims)/
                          reactionpid1/tottransmission);
    // 1% absolute error + 5% relative error
    double cserror = crosssection*pow(1./reactionpid1 + 2./reactionpid2+
                         pow(tottransmissionerror+ 0.02/tottransmission,2)+
                         pow(numberdensityerror,2), 0.5);  // Error on Transm.
    if(isnan(cserror)) cserror = 0;
    if(isnan(crosssection)) crosssection = 0;

    // make cross section string:
    setting set;
    stringstream stringout;

    stringout.precision(2);
    stringout << reaction << " \u03C3: " << std::scientific << 1E3*crosssection
              << " \u00b1 " << 1E3*cserror <<  "mb, CTS: " << std::defaultfloat
              << (reactionpid2.load()-chargestatevictims) << " of "
              << reactionpid1.load()/1. << ", T = " << tottransmission << " B\u03F1 = "
              << reactionpid2.load()/reactF5.at(1).GetEntries() << " CSC: "
              << chargestatevictims;
    if(set.isemptyortrans())
        stringout << " Ratio: " << 100.*(reactionpid2-chargestatevictims)/reactionpid1.load()
                  <<" \u00b1 "  << 100.*(reactionpid2-chargestatevictims)/reactionpid1.load()*
                     pow(1./reactionpid2+ 1./reactionpid1.load(),0.5) << " %";
    cout  << stringout.str() << endl;

    txtwriter txt;
    txt.addline(stringout.str());

    cout << "Ratio " << 100.*(reactionpid2-chargestatevictims)/reactionpid1.load()
         << "% " << endl;

    // Add Cross Section to the TGraph;
    tcross.SetPoint(tcross.GetN(), ncross, 1E3*crosssection);
    tcross.SetPointError(tcross.GetN()-1, 0, 1E3*cserror);
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
    if(reaction == "111NbPPN"){      incval = nancy::incval111Nb; targetval = nancy::targetval110Nb; }
    else if(reaction == "111NbPP2N"){incval = nancy::incval111Nb; targetval = nancy::targetval109Nb; }
    else if(reaction == "111NbP2P"){ incval = nancy::incval111Nb; targetval = nancy::targetval110Zr;
                                     binning   = 80; }
    else if(reaction == "110NbPPN"){ incval = nancy::incval110Nb; targetval = nancy::targetval109Nb;
                                     binning = 100; }
    else if(reaction == "110NbP2P"){ incval = nancy::incval110Nb; targetval = nancy::targetval109Zr;
                                     binning = 40; }
    else if(reaction == "110NbP0P"){ incval = nancy::incval110Nb; targetval = nancy::targetval110Nb; }
    else if(reaction == "110MoP3P"){ incval = nancy::incval110Mo; targetval = nancy::targetval108Zr; }
    else if(reaction == "111MoP3P"){ incval = nancy::incval111Mo; targetval = nancy::targetval109Zr; }
    else if(reaction == "112MoP3P"){ incval = nancy::incval112Mo; targetval = nancy::targetval110Zr; }
    else if(reaction == "113TcP3P"){ incval = nancy::incval113Tc; targetval = nancy::targetval111Nb; }
    else if(reaction == "112TcP3P"){ incval = nancy::incval112Tc; targetval = nancy::targetval110Nb; }
    else if(reaction == "114TcP3P"){ incval = nancy::incval114Tc; targetval = nancy::targetval112Nb; }
    else if(reaction == "113MoP3P"){ incval = nancy::incval113Mo; targetval = nancy::targetval111Zr; }
    else if(reaction == "110MoP2P"){ incval = nancy::incval110Mo; targetval = nancy::targetval109Nb; }
    else if(reaction == "111MoP2P"){ incval = nancy::incval111Mo; targetval = nancy::targetval110Nb; }
    else if(reaction == "112MoP2P"){ incval = nancy::incval112Mo; targetval = nancy::targetval111Nb; }
    else if(reaction == "90SeP2P"){  incval = nancy::incval90Se;  targetval = nancy::targetval89As; }
    else if(reaction == "90SeP3P"){  incval = nancy::incval90Se;  targetval = nancy::targetval88Ge; }
    else if(reaction == "89SeP2P"){  incval = nancy::incval89Se;  targetval = nancy::targetval88As; }
    else if(reaction == "89SeP3P"){  incval = nancy::incval89Se;  targetval = nancy::targetval87Ge; }
    else if(reaction == "88AsP2P"){  incval = nancy::incval88As;  targetval = nancy::targetval87Ge; }
    else if(reaction == "89AsP2P"){  incval = nancy::incval89As;  targetval = nancy::targetval88Ge; }
    else if(reaction == "89AsP3P"){  incval = nancy::incval89As;  targetval = nancy::targetval87Ga; }
    else if(reaction == "88GeP0P"){  incval = nancy::incval88Ge;  targetval = nancy::targetval88Ge; }
    else if(reaction == "89AsP0P"){  incval = nancy::incval89As;  targetval = nancy::targetval89As; }
    else if(reaction == "93BrP2P"){  incval = nancy::incval93Br;  targetval = nancy::targetval92Se; }
    else if(reaction == "93BrP3P"){  incval = nancy::incval93Br;  targetval = nancy::targetval91As; }
    else if(reaction == "94BrP2P"){  incval = nancy::incval94Br;  targetval = nancy::targetval93Se; }
    else if(reaction == "94BrP3P"){  incval = nancy::incval94Br;  targetval = nancy::targetval92As; }
    else if(reaction == "95BrP2P"){  incval = nancy::incval95Br;  targetval = nancy::targetval94Se; }
    else if(reaction == "95BrP3P"){  incval = nancy::incval95Br;  targetval = nancy::targetval93As; }
    else if(reaction == "94KrP2P"){  incval = nancy::incval94Kr;  targetval = nancy::targetval93Br; }
    else if(reaction == "94KrP3P"){  incval = nancy::incval94Kr;  targetval = nancy::targetval92Se; }
    else if(reaction == "95KrP2P"){  incval = nancy::incval95Kr;  targetval = nancy::targetval94Br; }
    else if(reaction == "95KrP3P"){  incval = nancy::incval95Kr;  targetval = nancy::targetval93Se; }
    else if(reaction == "96KrP2P"){  incval = nancy::incval96Kr;  targetval = nancy::targetval95Br; }
    else if(reaction == "96KrP3P"){  incval = nancy::incval96Kr;  targetval = nancy::targetval94Se; }
    else if(reaction == "97RbP2P"){  incval = nancy::incval97Rb;  targetval = nancy::targetval96Kr; }
    else if(reaction == "97RbP3P"){  incval = nancy::incval97Rb;  targetval = nancy::targetval95Br; }
    else if(reaction == "99RbP2P"){  incval = nancy::incval99Rb;  targetval = nancy::targetval98Kr; }
    else if(reaction == "99RbP3P"){  incval = nancy::incval99Rb;  targetval = nancy::targetval97Br; }
    else if(reaction == "100RbP2P"){ incval = nancy::incval100Rb; targetval = nancy::targetval99Kr; }
    else if(reaction == "100RbP3P"){ incval = nancy::incval100Rb; targetval = nancy::targetval98Br; }
    else if(reaction == "100SrP2P"){ incval = nancy::incval100Sr; targetval = nancy::targetval99Rb; }
    else if(reaction == "100SrP3P"){ incval = nancy::incval100Sr; targetval = nancy::targetval98Kr; }
    else if(reaction == "101SrP2P"){ incval = nancy::incval101Sr; targetval = nancy::targetval100Rb; }
    else if(reaction == "101SrP3P"){ incval = nancy::incval101Sr; targetval = nancy::targetval99Kr; }
    else if(reaction == "102SrP2P"){ incval = nancy::incval102Sr; targetval = nancy::targetval101Rb; }
    else if(reaction == "102SrP3P"){ incval = nancy::incval102Sr; targetval = nancy::targetval100Kr; }
    else if(reaction == "102YP2P"){  incval = nancy::incval102Y;  targetval = nancy::targetval101Sr; }
    else if(reaction == "102YP3P"){  incval = nancy::incval102Y;  targetval = nancy::targetval100Rb; }
    else if(reaction == "103YP2P"){  incval = nancy::incval103Y;  targetval = nancy::targetval102Sr; }
    else if(reaction == "103YP3P"){  incval = nancy::incval103Y;  targetval = nancy::targetval101Rb; }
    else if(reaction == "66MnP2P"){  incval = nancy::incval66Mn;  targetval = nancy::targetval65Cr;}
    else if(reaction == "66MnP3P"){  incval = nancy::incval66Mn;  targetval = nancy::targetval64V;}
    else if(reaction == "67MnP2P"){  incval = nancy::incval67Mn;  targetval = nancy::targetval66Cr;}
    else if(reaction == "67MnP3P"){  incval = nancy::incval67Mn;  targetval = nancy::targetval65V;}
    else if(reaction == "67FeP2P"){  incval = nancy::incval67Fe;  targetval = nancy::targetval66Mn;}
    else if(reaction == "67FeP3P"){  incval = nancy::incval67Fe;  targetval = nancy::targetval65Cr;}
    else if(reaction == "68FeP2P"){  incval = nancy::incval68Fe;  targetval = nancy::targetval67Mn;}
    else if(reaction == "68FeP3P"){  incval = nancy::incval68Fe;  targetval = nancy::targetval66Cr;}
    else if(reaction == "68CoP2P"){  incval = nancy::incval68Co;  targetval = nancy::targetval67Fe;}
    else if(reaction == "68CoP3P"){  incval = nancy::incval68Co;  targetval = nancy::targetval66Mn;}
    else if(reaction == "69CoP2P"){  incval = nancy::incval69Co;  targetval = nancy::targetval68Fe;}
    else if(reaction == "69CoP3P"){  incval = nancy::incval69Co;  targetval = nancy::targetval67Mn;}
    else if(reaction == "70CoP2P"){  incval = nancy::incval70Co;  targetval = nancy::targetval69Fe;}
    else if(reaction == "70CoP3P"){  incval = nancy::incval70Co;  targetval = nancy::targetval68Mn;}
    else if(reaction == "70NiP2P"){  incval = nancy::incval70Ni;  targetval = nancy::targetval69Co;}
    else if(reaction == "70NiP3P"){  incval = nancy::incval70Ni;  targetval = nancy::targetval68Fe;}
    else if(reaction == "71NiP2P"){  incval = nancy::incval71Ni;  targetval = nancy::targetval70Co;}
    else if(reaction == "71NiP3P"){  incval = nancy::incval71Ni;  targetval = nancy::targetval69Fe;}
    else if(reaction == "75ZnP2P"){  incval = nancy::incval75Zn;  targetval = nancy::targetval74Cu;}
    else if(reaction == "75ZnP3P"){  incval = nancy::incval75Zn;  targetval = nancy::targetval73Ni;}
    else if(reaction == "76ZnP2P"){  incval = nancy::incval76Zn;  targetval = nancy::targetval75Cu;}
    else if(reaction == "76ZnP3P"){  incval = nancy::incval76Zn;  targetval = nancy::targetval74Ni;}
    else if(reaction == "74CuP2P"){  incval = nancy::incval74Cu;  targetval = nancy::targetval73Ni;}
    else if(reaction == "74CuP3P"){  incval = nancy::incval74Cu;  targetval = nancy::targetval72Co;}
    else if(reaction == "75CuP2P"){  incval = nancy::incval75Cu;  targetval = nancy::targetval74Ni;}
    else if(reaction == "75CuP3P"){  incval = nancy::incval75Cu;  targetval = nancy::targetval73Co;}
    else if(reaction == "72NiP2P"){  incval = nancy::incval72Ni;  targetval = nancy::targetval71Co;}
    else if(reaction == "72NiP3P"){  incval = nancy::incval72Ni;  targetval = nancy::targetval70Fe;}
    else if(reaction == "73NiP2P"){  incval = nancy::incval73Ni;  targetval = nancy::targetval72Co;}
    else if(reaction == "73NiP3P"){  incval = nancy::incval73Ni;  targetval = nancy::targetval71Fe;}
    else if(reaction == "74NiP2P"){  incval = nancy::incval74Ni;  targetval = nancy::targetval73Co;}
    else if(reaction == "74NiP3P"){  incval = nancy::incval74Ni;  targetval = nancy::targetval72Fe;}
    else if(reaction == "82GeP2P"){  incval = nancy::incval82Ge;  targetval = nancy::targetval81Ga;}
    else if(reaction == "82GeP3P"){  incval = nancy::incval82Ge;  targetval = nancy::targetval80Zn;}
    else if(reaction == "83GeP2P"){  incval = nancy::incval83Ge;  targetval = nancy::targetval82Ga;}
    else if(reaction == "83GeP3P"){  incval = nancy::incval83Ge;  targetval = nancy::targetval81Zn;}
    else if(reaction == "80GaP2P"){  incval = nancy::incval80Ga;  targetval = nancy::targetval79Zn;}
    else if(reaction == "80GaP3P"){  incval = nancy::incval80Ga;  targetval = nancy::targetval78Cu;}
    else if(reaction == "81GaP2P"){  incval = nancy::incval81Ga;  targetval = nancy::targetval80Zn;}
    else if(reaction == "81GaP3P"){  incval = nancy::incval81Ga;  targetval = nancy::targetval79Cu;}
    else if(reaction == "82GaP2P"){  incval = nancy::incval82Ga;  targetval = nancy::targetval81Zn;}
    else if(reaction == "82GaP3P"){  incval = nancy::incval82Ga;  targetval = nancy::targetval80Cu;}
    else if(reaction == "78ZnP2P"){  incval = nancy::incval78Zn;  targetval = nancy::targetval77Cu;}
    else if(reaction == "78ZnP3P"){  incval = nancy::incval78Zn;  targetval = nancy::targetval76Ni;}
    else if(reaction == "79ZnP2P"){  incval = nancy::incval79Zn;  targetval = nancy::targetval78Cu;}
    else if(reaction == "79ZnP3P"){  incval = nancy::incval79Zn;  targetval = nancy::targetval77Ni;}
    else if(reaction == "80ZnP2P"){  incval = nancy::incval80Zn;  targetval = nancy::targetval79Cu;}
    else if(reaction == "80ZnP3P"){  incval = nancy::incval80Zn;  targetval = nancy::targetval78Ni;}
    else if(reaction == "81ZnP2P"){  incval = nancy::incval81Zn;  targetval = nancy::targetval80Cu;}
    //else if(reaction == "81ZnP3P"){  incval = nancy::incval81Zn;  targetval = nancy::targetval76Ni;}
    else if(reaction == "77CuP2P"){  incval = nancy::incval77Cu;  targetval = nancy::targetval76Ni;}
    else if(reaction == "77CuP3P"){  incval = nancy::incval77Cu;  targetval = nancy::targetval75Co;}
    else if(reaction == "78CuP2P"){  incval = nancy::incval78Cu;  targetval = nancy::targetval77Ni;}
    else if(reaction == "78CuP3P"){  incval = nancy::incval78Cu;  targetval = nancy::targetval76Co;}
    else if(reaction == "79CuP2P"){  incval = nancy::incval79Cu;  targetval = nancy::targetval78Ni;}
    //else if(reaction == "79CuP3P"){  incval = nancy::incval79Cu;  targetval = nancy::targetval75Co;}

    else __throw_invalid_argument(Form("Invalid reaction %s!\n", &reaction[0]));
}

void PID::chargestatecut(){
    setting set;
    if(set.isemptyortrans()){
        printf("Not doing charge state cut. No physics run.\n");
        return;
    }

    // Magic to get the charge state contribution
    vector<string> chemelem ={"R",
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
            "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
            "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc"};
    string numbers = "0123456789";
    int startindex = (int)reaction.find_first_not_of(numbers);
    projectileN = std::stoi(reaction.substr(0,startindex));
    ncross = projectileN; // Set Mass number for plotting

    // Figure out reaction:
    bool p3pbool = true;
    if( reaction.find("P2P") != string::npos) p3pbool = false;

    //Figure out starting Q
    int projectileZ = 0;
    string temp = reaction;
    temp.erase(temp.length()-3); // Erasing the P2/3P-part
    temp.erase(0,reaction.substr(0,startindex).size()); // Erase Number
    for(int i=1; i<chemelem.size(); i++){
        if(temp == chemelem.at(i)){ projectileZ = i; break;}
    }

    //cout << "Reaction: " << reaction << " Start number: " << projectileN << " Z " << projectileZ << endl;

    // Apply Reaction to values
    if (p3pbool){
         projectileN -= 2;
         projectileZ -= 2; }
    else{projectileN -= 1;
         projectileZ -= 1;}


    double aoq = projectileN/(double)projectileZ; // AoQ of daughter Nucleus
    double aoqoffset = aoq - targetval.at(0);

    cout <<" A/Q for daughter: " << aoq << ", off by: " << aoqoffset << endl;

    for(int i=-4; i<0; i++){
        double chargestateaoq = (projectileN+i)/(double)(projectileZ-1);
        if(abs(chargestateaoq-aoq) < 2*targetval.at(2)){ // Ellipsoids are touching
            cout << "For Q-1, Z = " << projectileN+i << " distance is " << chargestateaoq-aoq<< endl;
            // Calculate center of blob whose charge state is interfering
            double aoqintruder = (projectileN+i)/(double)projectileZ - aoqoffset;

            // Calculate intersecting ellipse of the part thats intruding
            double aoqintersect = aoqintruder + chargestateaoq - aoq;

            // Define upper and lower boundaries for ellipse overlap check
            vector<int> binlow{0,0,0};
            PIDplot.at(1).at(3).GetBinXYZ(
                PIDplot.at(1).at(3).FindBin(aoqintruder-targetval.at(2),
                                            targetval.at(1)-targetval.at(3)),
                binlow.at(0), binlow.at(1), binlow.at(2));
            vector<int> binhigh{0,0,0};
            PIDplot.at(1).at(3).GetBinXYZ(
                PIDplot.at(1).at(3).FindBin(aoqintruder+targetval.at(2),
                                            targetval.at(1)+targetval.at(3)),
                binhigh.at(0), binhigh.at(1), binhigh.at(2));

            //cout << "X-bin Width " << binhigh.at(0) << " - " << binlow.at(0)
            //     << ". Y-Bin Width " << binhigh.at(1) << " - " << binlow.at(1) << endl;

            // Prepare cut diagram.
            chargestate.emplace_back(TH2D(PIDplot.at(1).at(3)));
            chargestate.back().Reset();
            chargestate.back().SetName(Form("ChargestateN%i.%i",projectileN+i, projectileZ));
            chargestate.back().SetTitle(Form("CS contr. {}^{%i}%s{}^{%i+} Calc A/Q %.4f",
                                             projectileN +i, chemelem.at(projectileZ).c_str(),
                                             projectileZ-1,aoqintruder));

            // Loop over all bins and check for double inside ellipse
            for(int j=binlow.at(0); j<= binhigh.at(0); j++){ // x-loop
                double xcenter = PIDplot.at(1).at(3).GetXaxis()->GetBinCenter(j);
                for(int k=binlow.at(1); k<=binhigh.at(1); k++){
                    double ycenter = PIDplot.at(1).at(3).GetYaxis()->GetBinCenter(k);

                    if((pow(1./targetval.at(2)*(aoqintruder-xcenter),2) +
                        pow(1./targetval.at(3)*(targetval.at(1)-ycenter),2))<1 &&
                       (pow(1./targetval.at(2)*(aoqintersect-xcenter),2) +
                        pow(1./targetval.at(3)*(targetval.at(1)-ycenter),2))<1){
                        //cout << "Charge State X: " << j << " Y: " << k << " with: "
                        //     << PIDplot.at(1).at(3).GetBinContent(j,k) << " cts." << endl;
                        chargestate.back().SetBinContent(j,k,PIDplot.at(1).at(3).GetBinContent(j,k)+1E-9);
                    }
                }
            }
        }
    }

    // See all charge state contributions with 1% of initial value
    for(auto &elem: chargestate) chargestatevictims += elem.Integral()*0.01;
}

void PID::histogramsetup() {
    // Nasty histogram setup routine
    setting set;
    const vector<uint> y = set.getZrange(); // y boundaries
    const vector<double> x{500, 2.45, 2.85};  // xbins, lower x, upper x

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
    reactF5.emplace_back(TH1D("F5react", "F5-position of reacted particles",
                              binning,-100,100));

    const int brhoslice = 320;
    const vector<double> bl = {6.5,7.3};
    vector<TH1D> brt1, brt2;
    brt1.emplace_back(TH1D("F5brhop", "F5X projection B#rho",brhoslice,bl[0],bl[1]));
    brt1.emplace_back(TH1D("F7brhop", "F7X projection B#rho",brhoslice,bl[0],bl[1]));
    brt1.emplace_back(TH1D("F9brhop", "F9X projection B#rho",brhoslice,bl[0],bl[1]));
    brt1.emplace_back(TH1D("F11brhop", "F11X projection B#rho #sigma",brhoslice,bl[0],bl[1]));
    brt2.emplace_back(TH1D("F5brhostd", "F5X projection B#rho #sigma",brhoslice,bl[0],bl[1]));
    brt2.emplace_back(TH1D("F7brhostd", "F7X projection B#rho #sigma",brhoslice,bl[0],bl[1]));
    brt2.emplace_back(TH1D("F9brhostd", "F9X projection B#rho #sigma",brhoslice,bl[0],bl[1]));
    brt2.emplace_back(TH1D("F11brhostd", "F11X projection B#rho #sigma",brhoslice,bl[0],bl[1]));
    brhoprojection.emplace_back(brt1);
    brhoprojection.emplace_back(brt2);

    for(int j=0; j<4; j++){
        vector<TH2D> temp;
        string no = "F" + to_string(5+2*j); // form F5+F7+F9+F11
        int k =1; // Scaling factor F9
        if(j==2 || j==0) k=3; // make space wider for F9 and F5
        temp.emplace_back(TH2D((no+"pos").c_str(), ("PID "+no+" beamshape").c_str(),
                               250,-40*k,40*k,200,-30*k,30*k));
        temp.emplace_back(TH2D((no+"ang").c_str(), ("PID "+no+" beam angular shape").c_str(),
                               100,-50,50,100,-50,50));
        temp.emplace_back(TH2D((no+"brho").c_str(), ("PID "+no+" B#rho Distribution").c_str(),
                               250,-40*k,40*k,brhoslice,bl[0],bl[1]));
        reactPPAC.emplace_back(temp);
    }


    minosresults.emplace_back(TH2D("theta", "#theta correlation",
                                   90,0,90,90,0,90));
    minosresults.emplace_back(TH2D("tracknbr", "Track No. vs. Final Track No.",
                                   10,-0.5,9.5,10,-0.5,9.5));

    minosresults.at(0).GetXaxis()->SetTitle("#theta_{1} #circ");
    minosresults.at(0).GetYaxis()->SetTitle("#theta_{2} / #circ");
    minosresults.at(1).GetXaxis()->SetTitle("Track Number");
    minosresults.at(1).GetYaxis()->SetTitle("Final Track Number");

    minos1dresults.emplace_back(TH1D("zdistr", "Reaction distribution", 100,-70,130));
    minos1dresults.emplace_back(TH1D("phidistr", "Reaction angle distribution", 90, 0,180));

    minos1dresults.at(0).GetXaxis()->SetTitle("z / mm");
    minos1dresults.at(0).GetYaxis()->SetTitle("N");
    minos1dresults.at(1).GetXaxis()->SetTitle("#phi / #circ");
    minos1dresults.at(1).GetYaxis()->SetTitle("N");

    for(auto &i: minosresults) i.SetOption("colz");

    fitplot.emplace_back(TH2D("chisqfit","reduced #chi^{2}-fitrange", binning-1,1,
                              binning, binning-1,1,binning));
    fitplot.at(0).GetXaxis()->SetTitle("Starting Bin");
    fitplot.at(0).GetYaxis()->SetTitle("Number of Bins");
    fitplot.at(0).SetMaximum(2.*maxchisq);
    fitplot.emplace_back(TH2D("F9XF5X", "PID F9F5-X correlation", binning,-100,
                              100, 100,-120,120));
    fitplot.at(1).GetXaxis()->SetTitle("F5X");
    fitplot.at(1).GetYaxis()->SetTitle("F9X");
    for (auto &i: fitplot) i.SetOption("colz");

    for (auto &elem: reactF5){
        elem.GetXaxis()->SetTitle("x [mm]");
        elem.GetYaxis()->SetTitle("N");
    }

    for(int i=0; i<brhoprojection.at(0).size(); i++){
        string no = "F"+ to_string(5+2*i) + "X";
        brhoprojection.at(0).at(i).GetXaxis()->SetTitle("B#rho BigRIPS [Tm]");
        brhoprojection.at(1).at(i).GetXaxis()->SetTitle("B#rho BigRIPS [Tm]");
        brhoprojection.at(0).at(i).GetYaxis()->SetTitle(no.c_str());
        brhoprojection.at(1).at(i).GetYaxis()->SetTitle((no + " #sigma").c_str());
    }

    for(auto &elem: PIDplot){
        for(auto &uelem: elem){
            uelem.SetOption("colz");
            uelem.GetXaxis()->SetTitle("A/Q");
            uelem.GetYaxis()->SetTitle("Z");
            uelem.SetMinimum(1);
        }
    }

    for(int i=0; i<reactPPAC.size(); i++){
        string no = "F" + to_string(5+2*i);
        reactPPAC.at(i).at(0).GetXaxis()->SetTitle((no+"X [mm]").c_str());
        reactPPAC.at(i).at(0).GetYaxis()->SetTitle((no+"Y [mm]").c_str());
        reactPPAC.at(i).at(1).GetXaxis()->SetTitle((no+"A [mrad]").c_str());
        reactPPAC.at(i).at(1).GetYaxis()->SetTitle((no+"B [mrad]").c_str());
        reactPPAC.at(i).at(2).GetXaxis()->SetTitle((no+"X [mm]").c_str());
        reactPPAC.at(i).at(2).GetYaxis()->SetTitle("B#rho BigRIPS [Tm]");
    }

    for(auto &i: reactPPAC){
        for(auto &j :i){
            j.SetOption("colz");
            j.SetMinimum(1);
        }
    }
}

void PID::brhoprojections(){
    // Fill histograms of i.e. F5X(Brho), needs to be done afterwards
    for(int i=0; i<brhoprojection.at(0).size(); i++){ // Loop over all Focal Planes
        for(int j=1; j<reactPPAC.at(i).at(2).GetNbinsY(); j++){ // Loop over all bins
            TH1D* temp = reactPPAC.at(i).at(2).ProjectionX(Form("_pfx%i",j),j,j,"e");
            double mean = temp->GetMean();
            double stddev = temp->GetStdDev();
            brhoprojection.at(0).at(i).SetBinContent(j, mean);
            brhoprojection.at(0).at(i).SetBinError(j, temp->GetMeanError());
            brhoprojection.at(1).at(i).SetBinContent(j, stddev);
            brhoprojection.at(1).at(i).SetBinError(j, temp->GetStdDevError());
        }
    }
}