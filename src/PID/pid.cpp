//
// Created by afrotscher on 8/7/18.
//

#include "PID/pid.h"
#include "TF1.h"
#include "histogram_cuts.hh"
#include "minos.h"
#include "progress.h"
#include "txtwriter.h"
#include "zdssetting.h"
#include <algorithm>
#include <libconstant.h>
#include <sstream>
#include <thread>
#include <variant>

using std::vector, std::string, std::atomic, std::thread, std::cout, std::endl,
      std::stringstream, std::to_string, std::__throw_invalid_argument,
      std::min;

void PID::innerloop(treereader &tree, treereader &minostree, vector<uint> range,
                    const bool minosanalyse) {
    /// Step 1: duplicate the data structure
    decltype(reactF5)          _reactF5(reactF5);
    decltype(minosresults)     _minosresults(minosresults);
    decltype(minos1dresults)   _minos1dresults(minos1dresults);
    decltype(reactPPAC)        _reactPPAC(reactPPAC);
    decltype(PIDplot)          _PIDplot(PIDplot);
    decltype(minos3dresults)   _minos3dresults(minos3dresults);
    // MINOS with single events
    decltype(minossingleevent) _minossingleevent(minossingleevent);

    double brho_max = 10; // Tm, higher than all my values
    double brho_min = 6; // Tm, lower than all my values

    if(     reaction.find("P2P") != string::npos){
      brho_min = incval.at(5);
      brho_max = incval.at(6);
    }
    else if(reaction.find("P3P") != string::npos){
      brho_max = incval.at(8);
      brho_max = incval.at(9);
    }

    const uint threadno = range.at(0)/(range.at(1)-range.at(0));
    progressbar progress(range.at(1)-range.at(0), threadno);

    vector<vector<double>> valinc;     // Store temporary beam values

    //F5-2B, F7-2B, F9-2B, F11-2B
    const std::array<int, 4> ppacFpositions{12,17,25,33};
    const std::array<int, 4> ppacangledistance{650,945,700,500}; // in mm

    /// 2. Do the event loop
    for(int eventcounter=range.at(0); eventcounter<range.at(1); eventcounter++){
        /// Don't take sorted out events
        progress.increaseevent();
        if(!goodevents.at(eventcounter).at(0)) continue;

        tree.getevent(eventcounter);
        const double beamaoqcorr = tree.BigRIPSBeam_aoq[0] + p1.F7absF5X0 -
                (p1.F7absF5X+tree.F5X*p1.F7linF5X) -
                tree.F5A*p1.F7linF5A -
                tree.F3X*p1.F7linF3X;

        const double beamaoqcorr2 = tree.BigRIPSBeam_aoq[4] +  p1.F11absF9X0-
                (p1.F11absF9X+tree.F9X*p1.F11linF9X) -
                tree.F9A*p1.F11linF9A -
                tree.F11A*p1.F11linF11A;

        //Loop over all elements in the tree
        valinc.push_back({tree.BigRIPSBeam_aoq[0],tree.BigRIPSBeam_aoq[1]});
        valinc.push_back({tree.BigRIPSBeam_zet[0],tree.BigRIPSBeam_zet[1]});
        if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
            // Cut Particles that have variating aoq or zet
            _PIDplot.at(0).at(0).Fill(tree.BigRIPSBeam_aoq[0],
                    tree.BigRIPSBeam_zet[0]);
            _PIDplot.at(1).at(0).Fill(beamaoqcorr, tree.BigRIPSBeam_zet[0]);
        }

        valinc.push_back({tree.BigRIPSBeam_aoq[2],tree.BigRIPSBeam_aoq[3],
                          tree.BigRIPSBeam_aoq[4]});
        valinc.push_back({tree.BigRIPSBeam_zet[2],tree.BigRIPSBeam_zet[3],
                          tree.BigRIPSBeam_zet[4]});

        // F7-F11 part
        if(closeness(valinc.at(2)) && closeness(valinc.at(3)) &&
            goodevents.at(eventcounter).at(1)){
            _PIDplot.at(0).at(1).Fill(tree.BigRIPSBeam_aoq[4],
                    tree.BigRIPSBeam_zet[4]);
            _PIDplot.at(1).at(1).Fill(beamaoqcorr2,tree.BigRIPSBeam_zet[4]);
        }
        valinc.clear();

        // We now fill the cut data (cut by ellipsoid and brho)
        if((pow(1./incval.at(2)*(beamaoqcorr-incval.at(0)),2) +
            pow(1/incval.at(3)*(tree.BigRIPSBeam_zet[0]-incval.at(1)),2))
            >= 1 ) continue;

        _PIDplot.at(0).at(2).Fill(tree.BigRIPSBeam_aoq[0],
                tree.BigRIPSBeam_zet[0]);
        _PIDplot.at(1).at(2).Fill(beamaoqcorr,tree.BigRIPSBeam_zet[0]);

        if(goodevents.at(eventcounter).at(1)) { // F7-F11 part
            _PIDplot.at(0).at(3).Fill(tree.BigRIPSBeam_aoq[4],
                    tree.BigRIPSBeam_zet[4]);
            _PIDplot.at(1).at(3).Fill(beamaoqcorr2, tree.BigRIPSBeam_zet[4]);
        }
        // cut on energy(Brho) to avoid offcenter Transmission
        if(tree.BigRIPSRIPS_brho[1] < brho_max &&
           tree.BigRIPSRIPS_brho[1] > brho_min) reactionpid1++;

        // Fill F5 value with PID
        _reactF5.at(0).Fill(tree.F5X);

        // Second ellipsoid for cross section F7-F11 part
        if(!((pow(1./targetval.at(2)*(beamaoqcorr2-targetval.at(0)),2) +
              pow(1./targetval.at(3)*(tree.BigRIPSBeam_zet[4]-
                                      targetval.at(1)),2))<1
            && goodevents.at(eventcounter).at(1))) continue;

        // cut on energy(Brho) to avoid offcenter Transmission
        if(tree.BigRIPSRIPS_brho[1] < brho_max &&
           tree.BigRIPSRIPS_brho[1] > brho_min) reactionpid2++;
        // Investigate F7 position of (p,2p) ions (off center effects)
        _reactF5.at(1).Fill(tree.F5X);

        //Check double cut particles
        _PIDplot.at(0).at(4).Fill(tree.BigRIPSBeam_aoq[4],
                                  tree.BigRIPSBeam_zet[4]);
        _PIDplot.at(1).at(4).Fill(beamaoqcorr2, tree.BigRIPSBeam_zet[4]);

        // Fill Information to the beam profile/shape/energy
        // 17 == F7PPAC-2B, 15 == F7PPAC-1B
        for(int i=0; i<_reactPPAC.size();i++){
            _reactPPAC.at(i).at(0).Fill(
                    tree.BigRIPSPPAC_fX[ppacFpositions.at(i)],
                    tree.BigRIPSPPAC_fY[ppacFpositions.at(i)]);
            _reactPPAC.at(i).at(1).Fill(
                    1E3*atan2(tree.BigRIPSPPAC_fX[ppacFpositions.at(i)]-
                              tree.BigRIPSPPAC_fX[ppacFpositions.at(i)-2],
                              ppacangledistance.at(i)),
                    1E3*atan2(tree.BigRIPSPPAC_fY[ppacFpositions.at(i)]-
                              tree.BigRIPSPPAC_fY[ppacFpositions.at(i)-2],
                              ppacangledistance.at(i)));
            _reactPPAC.at(i).at(2).Fill(
                    tree.BigRIPSPPAC_fX[ppacFpositions.at(i)],
                    tree.BigRIPSRIPS_brho[1]);
        }
        // do the correlation between F5 and F9
        fitplot.at(1).Fill(tree.F5X, tree.F9X);

        // make the mighty minos analysis
        if(reaction.find("P0P") == string::npos && minosanalyse){
            minostree.getevent(eventcounter);
            minosana analysis(minostree.Trackamount, minostree.Tshaping,
                    minostree.TimeBinElec, minostree.DelayTrig,
                    minostree.VDrift, minostree.minostrackxy,
                    minostree.Minoscalibvalues, minostree.minostime,
                    minostree.MinosClustX, minostree.MinosClustY,
                    minostree.MinosClustQ, (int)threadno, _minossingleevent);
            TMinosPass minres = analysis.analyze();

            if(minres.trackNbr_final == 2){
                _minosresults.at(0).Fill(minres.thetaz[0], minres.thetaz[1]);
                _minosresults.at(3).Fill(minres.chargeweight[0],
                                         minres.chargeweight[1]);
                _minos1dresults.at(1).Fill(minres.phi_vertex.at(0));
            }

            _minosresults.at(1).Fill(minres.trackNbr,minres.trackNbr_final);

            if(minres.trackNbr_final == 3){
                _minosresults.at(2).Fill(minres.lambda2d.at(0),
                                         minres.lambda2d.at(1));
                _minos3dresults.at(0).Fill(minres.thetaz[0], minres.thetaz[1],
                                           minres.thetaz[2]);
                for(auto &i: minres.phi_vertex) _minos1dresults.at(5).Fill(i);
            }

            if(minres.trackNbr_final > 1 ){
                // at least two tracks need to survive
                _minos1dresults.at(0).Fill(minres.z_vertex);
            }

            // Get Vertex data
            if(reaction.find("P2P") != string::npos &&
               minres.trackNbr_final == 2) { // P2P-mode && 2 Tracks
              _minos1dresults.at(2).Fill(minres.vertexdist[0]);

              // Get all theta angles as 1d
              for(auto &i: minres.thetaz) _minos1dresults.at(4).Fill(i);
            }

            if(reaction.find("P3P") != string::npos &&
               minres.trackNbr_final == 3) { // P3P-mode && 3 Tracks
              for (auto &i : minres.vertexdist)
                _minos1dresults.at(2).Fill(i);

              // Get all theta angles as 1d
              for(auto &i: minres.thetaz) _minos1dresults.at(4).Fill(i);
            }

            // Get Lambda Errors
            for(auto &i: minres.lambda2dE) _minos1dresults.at(3).Fill(i);
        }
    }

    /// Step 3: rejoining data structure
    unitemutex.lock();
    for(uint i=0;i<reactF5.size();i++) reactF5.at(i).Add(&_reactF5.at(i));

    for(uint i=0;i<minosresults.size();i++)
        minosresults.at(i).Add(&_minosresults.at(i));

    for(uint i=0; i<minos1dresults.size(); i++)
        minos1dresults.at(i).Add(&_minos1dresults.at(i));

    for(uint i=0; i<minos3dresults.size(); i++)
        minos3dresults.at(i).Add(&_minos3dresults.at(i));

    for(uint i=0;i<reactPPAC.size();i++){
        for(uint j=0; j<reactPPAC.at(0).size(); j++){
            reactPPAC.at(i).at(j).Add(&_reactPPAC.at(i).at(j));
        }
    }

    for(uint i=0;i<PIDplot.size();i++){
        for(uint j=0;j<PIDplot.at(0).size();j++){
            PIDplot.at(i).at(j).Add(&_PIDplot.at(i).at(j));
        }
    }
    minossingleevent.insert(minossingleevent.end(), _minossingleevent.begin(),
                            _minossingleevent.end());

    unitemutex.unlock();
    progressbar::reset();
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

    if(setting::getminos())
        threads = std::min(25, std::max((int)sqrt(goodevents.size())/450,2));

    printf("Making PID with %i threads now...\n", threads);

    /// Make regular TChain for fast input
    vector<treereader> tree;
    tree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) tree.emplace_back(input);

    vector<string> keys{"BigRIPSBeam.aoq", "BigRIPSBeam.zet", "F3X", "F3A",
                        "F5X", "F5A", "F9X", "F9A", "F11X", "F11A",
                        "BigRIPSPPAC.fX", "BigRIPSPPAC.fY", "BigRIPSRIPS.brho"};
    for(auto &i:tree) i.setloopkeys(keys);

    /// Make slow chain for

    vector<treereader> minostree;
    minostree.reserve(threads); // MUST stay as reallocation will call d'tor
    for(int i=0; i<threads; i++) minostree.emplace_back(input);

    vector<string> minoskeys{"VDrift", "MinosClustX", "MinosClustY",
                             "MinosClustQ", "DelayTrig", "Trackamount",
                             "Minoscalibvalues", "TimeBinElec", "Tshaping",
                             "minostrackxy", "minostime"};
    for(auto &i: minostree) i.setloopkeys(minoskeys);

    //Constructing all the histograms
    PID::histogramsetup();

    //Making threads
    vector<thread> th;
    th.reserve(threads); // Make space for $threads Threads
    for(uint i=0; i<threads;i++){
        vector<uint> ranges = {(uint)(i*goodevents.size()/threads),
                               (uint)((i+1)*goodevents.size()/threads-1)};
        th.emplace_back(thread(&PID::innerloop, this, std::ref(tree.at(i)),
                               std::ref(minostree.at(i)), ranges,
                               setting::getminos()));
    }

    for(auto &t : th) t.detach();

    progressbar finishcondition;
    while(progressbar::ongoing()) finishcondition.draw();

    // Calculate the transmission and the cross section
    PID::offctrans();
    PID::chargestatecut();
    PID::crosssection();
    PID::brhoprojections();

    for(int i=0; i<2; i++) {
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
        for(auto &j : brhoprojection) j.at(i).Write();
    }

    output->cd(folders.at(7).c_str());
    for(auto &elem: chargestate) elem.Write();

    output->cd(folders.at(8).c_str());
    for(auto &elem: minosresults) elem.Write();
    for(auto &elem: minos1dresults) elem.Write();
    for(auto &elem: minos3dresults) elem.Write();

    output->cd(folders.at(9).c_str());
    for(auto &elem: minossingleevent) elem.Write();

    //Clean up tree
    for(auto &i: fitstyle) for(auto &j:i) j->Delete();
}

void PID::offctrans() {
    // Calculating raw beam ratio:
    reactF5.push_back(reactF5.at(1));
    for (auto &i:reactF5) i.Sumw2();
    reactF5.back().Divide(&reactF5.at(0));
    reactF5.back().SetName("F5ratio");

    if(setting::isemptyortrans()){
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
            auto corrfit= new TF1(Form("Fit %i-%i", i,j),constfit,
                    reactF5.back().GetBinCenter(i),
                    reactF5.back().GetBinCenter(j),1);
            reactF5.back().Fit(corrfit, "NRQ");
            if(corrfit->GetNDF() >= (minrange - 1) &&
               corrfit->GetParameter(0) > 0.95*mean)
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

    // switch for effect (current normal)
    for(int i= (fitplot.at(0).GetNbinsY()-1); i>=minrange; i--){
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
    if(     reaction.find("P2P") != string::npos)
            tottransmission = incval.at(4);
    else if(reaction.find("P3P") != string::npos)
            tottransmission = incval.at(7);

    // Calculate the final crosssection for this run
    const double numberdensity = 0.433; // atoms/cm2 for 10cm LH2
    const double numberdensityerror = 0.00924; // relative error
    const double tottransmissionerror = 0.05; // conservative relative error
    const double rel_num_dens_crap = 0.046; // relative amount of atoms F7<->D8

    // Reduce CS contribution by Brho cut
    if(reactF5.at(1).GetEntries() > 0)
        chargestatevictims *= reactionpid2.load()/reactF5.at(1).GetEntries();

    double crosssection = std::max(0.,
                             1./numberdensity*(reactionpid2-chargestatevictims)/
                             reactionpid1/tottransmission*
                             (1-rel_num_dens_crap/(1+rel_num_dens_crap)));
    // 1% absolute error + 5% relative error
    double cserror = sqrt(
        pow(sqrt((double)reactionpid2)/(numberdensity*reactionpid1*tottransmission),2)+
        pow(sqrt(chargestatevictims)/(numberdensity*reactionpid1*tottransmission),2)+
        pow(crosssection/sqrt((double)reactionpid1),2) +
        pow((tottransmissionerror + 0.02*tottransmission)*crosssection, 2) +
        pow(numberdensityerror*crosssection, 2)
        );

    if(isnan(cserror)) cserror = 0;
    if(isnan(crosssection)) crosssection = 0;

    // make cross section string:
    stringstream stringout;

    stringout.precision(2);
    stringout << reaction << " \u03C3: " << std::scientific << 1E3*crosssection
              << " \u00b1 " << 1E3*cserror <<  "mb, CTS: " << std::defaultfloat
              << (reactionpid2.load()-chargestatevictims) << " of "
              << reactionpid1.load()/1. << ", T = " << tottransmission << " B\u03F1 = "
              << reactionpid2.load()/reactF5.at(1).GetEntries() << " CSC: "
              << chargestatevictims;
    if(setting::isemptyortrans())
        stringout
            << " Ratio: " << 100.*(reactionpid2-chargestatevictims)/reactionpid1.load()
            <<" \u00b1 "  << 100.*(reactionpid2-chargestatevictims)/reactionpid1.load()*
                pow(1./reactionpid2+ 1./reactionpid1.load(),0.5) << " %";
    cout  << stringout.str() << endl;

    txtwriter::addline(stringout.str());

    cout << "Ratio " << 100.*(reactionpid2-chargestatevictims)/reactionpid1.load()
         << "% " << endl;

    // Add Cross Section to the TGraph;
    tcross.SetPoint(tcross.GetN(), ncross, 1E3*crosssection);
    tcross.SetPointError(tcross.GetN()-1, 0, 1E3*cserror);
}

void PID::reactionparameters() {
    // Setup cut values
    if(setting::isemptyortrans()){
        incval = setting::getPIDincutvalue();
        targetval = setting::getPIDoutcutvalue();
        reaction = setting::getmodename();
        return;
    }

    using namespace nancy;

    vector<std::array< std::variant<string, a4, a10>,3>> master{
        {"111NbPPN", incval111Nb, targetval110Nb},
        {"111NbPP2N",incval111Nb, targetval109Nb},
        {"111NbP2P", incval111Nb, targetval110Zr},
        {"110NbPPN", incval110Nb, targetval109Nb},
        {"110NbP2P", incval110Nb, targetval109Zr},
        {"110NbP0P", incval110Nb, targetval110Nb},
        {"110MoP3P", incval110Mo, targetval108Zr},
        {"111MoP3P", incval111Mo, targetval109Zr},
        {"112MoP3P", incval112Mo, targetval110Zr},
        {"113TcP3P", incval113Tc, targetval111Nb},
        {"112TcP3P", incval112Tc, targetval110Nb},
        {"114TcP3P", incval114Tc, targetval112Nb},
        {"113MoP3P", incval113Mo, targetval111Zr},
        {"110MoP2P", incval110Mo, targetval109Nb},
        {"111MoP2P", incval111Mo, targetval110Nb},
        {"112MoP2P", incval112Mo, targetval111Nb},
        {"90SeP2P", incval90Se, targetval89As},
        {"90SeP3P", incval90Se, targetval88Ge},
        {"89SeP2P", incval89Se, targetval88As},
        {"89SeP3P", incval89Se, targetval87Ge},
        {"88AsP2P", incval88As, targetval87Ge},
        {"89AsP2P", incval89As, targetval88Ge},
        {"89AsP3P", incval89As, targetval87Ga},
        {"88GeP0P", incval88Ge, targetval88Ge},
        {"89AsP0P", incval89As, targetval89As},
        {"93BrP2P", incval93Br, targetval92Se},
        {"93BrP3P", incval93Br, targetval91As},
        {"94BrP2P", incval94Br, targetval93Se},
        {"94BrP3P", incval94Br, targetval92As},
        {"95BrP2P", incval95Br, targetval94Se},
        {"95BrP3P", incval95Br, targetval93As},
        {"94KrP2P", incval94Kr, targetval93Br},
        {"94KrP3P", incval94Kr, targetval92Se},
        {"95KrP2P", incval95Kr, targetval94Br},
        {"95KrP3P", incval95Kr, targetval93Se},
        {"96KrP2P", incval96Kr, targetval95Br},
        {"96KrP3P", incval96Kr, targetval94Se},
        {"97RbP2P", incval97Rb, targetval96Kr},
        {"97RbP3P", incval97Rb, targetval95Br},
        {"99RbP2P", incval99Rb, targetval98Kr},
        {"99RbP3P", incval99Rb, targetval97Br},
        {"100RbP2P", incval100Rb, targetval99Kr},
        {"100RbP3P", incval100Rb, targetval98Br},
        {"100SrP2P", incval100Sr, targetval99Rb},
        {"100SrP3P", incval100Sr, targetval98Kr},
        {"101SrP2P", incval101Sr, targetval100Rb},
        {"101SrP3P", incval101Sr, targetval99Kr},
        {"102SrP2P", incval102Sr, targetval101Rb},
        {"102SrP3P", incval102Sr, targetval100Kr},
        {"102YP2P", incval102Y, targetval101Sr},
        {"102YP3P", incval102Y, targetval100Rb},
        {"103YP2P", incval103Y, targetval102Sr},
        {"103YP3P", incval103Y, targetval101Rb},
        {"66MnP2P", incval66Mn, targetval65Cr},
        {"66MnP3P", incval66Mn, targetval64V},
        {"67MnP2P", incval67Mn, targetval66Cr},
        {"67MnP3P", incval67Mn, targetval65V},
        {"67FeP2P", incval67Fe, targetval66Mn},
        {"67FeP3P", incval67Fe, targetval65Cr},
        {"68FeP2P", incval68Fe, targetval67Mn},
        {"68FeP3P", incval68Fe, targetval66Cr},
        {"68CoP2P", incval68Co, targetval67Fe},
        {"68CoP3P", incval68Co, targetval66Mn},
        {"69CoP2P", incval69Co, targetval68Fe},
        {"69CoP3P", incval69Co, targetval67Mn},
        {"70CoP2P", incval70Co, targetval69Fe},
        {"70CoP3P", incval70Co, targetval68Mn},
        {"70NiP2P", incval70Ni, targetval69Co},
        {"70NiP3P", incval70Ni, targetval68Fe},
        {"71NiP2P", incval71Ni, targetval70Co},
        {"71NiP3P", incval71Ni, targetval69Fe},
        {"75ZnP2P", incval75Zn, targetval74Cu},
        {"75ZnP3P", incval75Zn, targetval73Ni},
        {"76ZnP2P", incval76Zn, targetval75Cu},
        {"76ZnP3P", incval76Zn, targetval74Ni},
        {"74CuP2P", incval74Cu, targetval73Ni},
        {"74CuP3P", incval74Cu, targetval72Co},
        {"75CuP2P", incval75Cu, targetval74Ni},
        {"75CuP3P", incval75Cu, targetval73Co},
        {"72NiP2P", incval72Ni, targetval71Co},
        {"72NiP3P", incval72Ni, targetval70Fe},
        {"73NiP2P", incval73Ni, targetval72Co},
        {"73NiP3P", incval73Ni, targetval71Fe},
        {"74NiP2P", incval74Ni, targetval73Co},
        {"74NiP3P", incval74Ni, targetval72Fe},
        {"82GeP2P", incval82Ge, targetval81Ga},
        {"82GeP3P", incval82Ge, targetval80Zn},
        {"83GeP2P", incval83Ge, targetval82Ga},
        {"83GeP3P", incval83Ge, targetval81Zn},
        {"80GaP2P", incval80Ga, targetval79Zn},
        {"80GaP3P", incval80Ga, targetval78Cu},
        {"81GaP2P", incval81Ga, targetval80Zn},
        {"81GaP3P", incval81Ga, targetval79Cu},
        {"82GaP2P", incval82Ga, targetval81Zn},
        {"82GaP3P", incval82Ga, targetval80Cu},
        {"78ZnP2P", incval78Zn, targetval77Cu},
        {"78ZnP3P", incval78Zn, targetval76Ni},
        {"79ZnP2P", incval79Zn, targetval78Cu},
        {"79ZnP3P", incval79Zn, targetval77Ni},
        {"80ZnP2P", incval80Zn, targetval79Cu},
        {"80ZnP3P", incval80Zn, targetval78Ni},
        {"81ZnP2P", incval81Zn, targetval80Cu},
        {"77CuP2P", incval77Cu, targetval76Ni},
        {"77CuP3P", incval77Cu, targetval75Co},
        {"78CuP2P", incval78Cu, targetval77Ni},
        {"78CuP3P", incval78Cu, targetval76Co},
        {"79CuP2P", incval79Cu, targetval78Ni}};

    for(auto &i: master){
        if(reaction == std::get<string>(i.at(0))){
            incval    = std::get<a10>(i.at(1));
            targetval = std::get<a4>(i.at(2));
            return;
        }
    }

     __throw_invalid_argument(Form("Invalid reaction %s!\n", &reaction[0]));
}

void PID::chargestatecut(){
    if(setting::isemptyortrans()){
        printf("Not doing charge state cut. No physics run.\n");
        return;
    }

    // Magic to get the charge state contribution
    vector<string> chemelem ={"R",
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
            "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
            "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
            "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc"};
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

    //cout << "Reaction: " << reaction << " Start number: "
    //     << projectileN << " Z " << projectileZ << endl;

    // Apply Reaction to values
    if (p3pbool){
         projectileN -= 2;
         projectileZ -= 2;
    }
    else {
      projectileN -= 1;
      projectileZ -= 1;
    }

    double aoq = projectileN/(double)projectileZ; // AoQ of daughter Nucleus
    double aoqoffset = aoq - targetval.at(0);

    cout <<" A/Q for daughter: " << aoq << ", off by: " << aoqoffset << endl;

    for(int i=-4; i<0; i++){
        double chargestateaoq = (projectileN+i)/(double)(projectileZ-1);
        // Ellipsoids are touching
        if(abs(chargestateaoq-aoq) < 2*targetval.at(2)){
            cout << "For Q-1, Z = " << projectileN+i << " distance is "
                 << chargestateaoq-aoq<< endl;
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

            //cout << "X-bin " << binhigh.at(0) << " - " << binlow.at(0)
            //     << " Y-Bin " << binhigh.at(1) << " - " << binlow.at(1) << endl;

            // Prepare cut diagram.
            chargestate.emplace_back(PIDplot.at(1).at(3));
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
                        chargestate.back().SetBinContent(j,k,PIDplot.at(1).at(3).GetBinContent(j,k));
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
    const vector<uint> y = setting::getZrange(); // y boundaries
    const std::array<double, 3> x{500, 2.45, 2.85};  // xbins, lower x, upper x

    vector<vector<string>> t1s = {{"pidinc", "PID Incoming F3-F7"},
        {"pidout", "PID Outgoing F8-F11"}, {"pidincut","PID Inc cut  F3-F7"},
        {"pidincutout", "PID Out w/ in cut"},
        {"pidincutoutcut", "PID Out cut w/ in cut"}};
    vector<vector<string>> t2s = {{"pidinccorr", "PID Incoming F3-F7"},
        {"pidoutcorr", "PID Outgoing F8-F11"},{"pidincutcorr","PID Inc cut  F3-F7"},
        {"pidincutoutcorr", "PID Out w/ in cut"},
        {"pidincutoutcutcorr", "PID Out cut w/ in cut"}};

    PIDplot.emplace_back();
    for(auto &i : t1s) PIDplot.at(0).emplace_back(i.at(0).c_str(),i.at(1).c_str(),
                                          x[0],x[1],x[2],y[0],y[1],y[2]);
    PIDplot.emplace_back();
    for(auto &i : t2s) PIDplot.at(1).emplace_back(i.at(0).c_str(),i.at(1).c_str(),
                                          x[0],x[1],x[2],y[0],y[1],y[2]);

    reactF5 = {{"F5beam", "F5 beam profile", binning,-100,100},
               {"F5react", "F5-position of reacted particles", binning,-100,100}};

    const int brhoslice = 320;
    const std::array<double, 2> bl = {6.5,7.3};

    string t = "projection B#rho #sigma";
    brhoprojection =
        {{{"F5brhop", "F5X projection B#rho", brhoslice, bl[0], bl[1]},
          {"F7brhop", "F7X projection B#rho", brhoslice, bl[0], bl[1]},
          {"F9brhop", "F9X projection B#rho", brhoslice, bl[0], bl[1]},
          {"F11brhop", "F11X projection B#rho", brhoslice, bl[0], bl[1]}},
         {{"F5brhostd", Form("F5X %s", t.c_str()), brhoslice, bl[0], bl[1]},
          {"F7brhostd", Form("F7X %s", t.c_str()), brhoslice, bl[0], bl[1]},
          {"F9brhostd", Form("F9X %s", t.c_str()), brhoslice, bl[0], bl[1]},
          {"F11brhostd", Form("F11X %s", t.c_str()), brhoslice, bl[0], bl[1]}}};

    for(int j=0; j<4; j++){
        reactPPAC.emplace_back();
        string no = "F" + to_string(5+2*j); // form F5+F7+F9+F11
        int k =1; // Scaling factor F9
        if(j==2 || j==0) k=3; // make space wider for F9 and F5

        reactPPAC.back() =
            {{(no+"pos").c_str(), ("PID "+no+" beamshape").c_str(),
                                        250, -40.*k, 40.*k, 200, -30.*k, 30.*k},
             {(no+"ang").c_str(), ("PID "+no+" beam angular shape").c_str(),
                                                    100, -50, 50, 100, -50, 50},
             {(no+"brho").c_str(), ("PID "+no+" B#rho Distribution").c_str(),
                                  250, -40.*k, 40.*k, brhoslice, bl[0], bl[1]}};
    }

    minos3dresults = {{"theta", "#theta 3D corr.", 90,0,90,90,0,90,90,0,90}};
    minos3dresults.at(0).GetXaxis()->SetTitle("#theta_{min} / #circ");
    minos3dresults.at(0).GetYaxis()->SetTitle("#theta_{med} / #circ");
    minos3dresults.at(0).GetZaxis()->SetTitle("#theta_{max} / #circ");

    minosresults =
        {{"theta", "#theta correlation", 90, 0, 90, 90, 0, 90},
         {"tracknbr", "Track No. vs. Final Track No.", 10,-0.5,9.5,10,-0.5,9.5},
         {"lambda", "Two smaller 2d angles for p,3p", 60, 0, 120, 90, 0, 180},
         {"minoscharge", "Charge deposition of protons",150,0,1500,150,0,1500}};

    minosresults.at(0).GetXaxis()->SetTitle("#theta_{1} #circ");
    minosresults.at(0).GetYaxis()->SetTitle("#theta_{2} / #circ");
    minosresults.at(1).GetXaxis()->SetTitle("Track Number");
    minosresults.at(1).GetYaxis()->SetTitle("Final Track Number");
    minosresults.at(2).GetXaxis()->SetTitle("smaller angle / #circ");
    minosresults.at(2).GetYaxis()->SetTitle("larger angle / #circ");
    minosresults.at(3).GetXaxis()->SetTitle("Q_{1}/l_{TPC} AU");
    minosresults.at(3).GetYaxis()->SetTitle("Q_{2}/l_{TPC} AU");

    for(auto &i: minosresults) i.SetOption("colz");

    minos1dresults =
        {{"zdistr", "Reaction distribution", 100,-70,130},
         {"phidistr2", "Reaction angle two tracks", 90, 0, 180},
         {"vertexdistr", "Distance between proton tracks", 300, 0, 50},
         {"lambdaE", "Angular Error of tracks", 300,0, 15},
         {"theta", "#theta of all protons", 180,0, 180},
         {"phidistr3", "Reaction angle three tracks", 90, 0, 180}};

    minos1dresults.at(0).GetXaxis()->SetTitle("z / mm");
    minos1dresults.at(1).GetXaxis()->SetTitle("#phi / #circ");
    minos1dresults.at(2).GetXaxis()->SetTitle("#Delta x_{vertex} / mm");
    minos1dresults.at(3).GetXaxis()->SetTitle("#Delta #lambda_{track} / #circ");
    minos1dresults.at(4).GetXaxis()->SetTitle("#theta / #circ");
    minos1dresults.at(5).GetXaxis()->SetTitle("#phi / #circ");

    for(auto &elem: minos1dresults) elem.GetYaxis()->SetTitle("N");

    fitplot =
        {{"chisqfit","reduced #chi^{2}-fitrange", binning-1,1,
                   (double)binning, binning-1,1,(double)binning},
         {"F9XF5X", "PID F9F5-X correlation", binning,-100, 100, 100,-120,120}};

    fitplot.at(0).GetXaxis()->SetTitle("Starting Bin");
    fitplot.at(0).GetYaxis()->SetTitle("Number of Bins");
    fitplot.at(0).SetMaximum(2.*maxchisq);
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
        const string no = "F" + to_string(5+2*i);
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
            TH1D *temp = reactPPAC.at(i).at(2).ProjectionX(Form("_pfx%i",j),j,j,"e");
            double mean = temp->GetMean();
            double stddev = temp->GetStdDev();
            brhoprojection.at(0).at(i).SetBinContent(j, mean);
            brhoprojection.at(0).at(i).SetBinError(j, temp->GetMeanError());
            brhoprojection.at(1).at(i).SetBinContent(j, stddev);
            brhoprojection.at(1).at(i).SetBinError(j, temp->GetStdDevError());
            temp->Delete();
        }
    }
}