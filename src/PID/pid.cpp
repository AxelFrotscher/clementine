//
// Created by afrotscher on 8/7/18.
//

#include <libconstant.h>
#include "PID/pid.h"
#include "thread"
#include "histogram_cuts.hh"
#include "MakeAllTree_78Ni.hh"

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

    const int downscale = (int)((range.at(1)-range.at(0))/100.); // every n-th event
    int threadno = range.at(0)/(range.at(1)-range.at(0));

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
                }
            }
        }
        if(!((eventcounter-range.at(0))%downscale)){
            consolemutex.lock();
            progressbar(eventcounter-range.at(0),range.at(1)-range.at(0),threadno);
            consolemutex.unlock();
        }
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

}

void PID::analyse(std::vector <std::string> input, TFile *output) {
    // Progress Bar setup
    printf("Making PID with %i threads now...\n", threads);

    vector<TChain*> chain;
    for(int i=0; i<threads; i++){
        chain.emplace_back(new TChain("tree"));
        for(auto h: input) chain.back()->Add(h.c_str());
    }

    vector<treereader*> tree;
    for(auto *i:chain){
        tree.emplace_back(new treereader(i));
    }

    vector<string> keys{"BigRIPSBeam.aoq", "BigRIPSBeam.zet", "F3X", "F3A",
                         "F5X", "F5A", "F9X", "F9A", "F11X", "F11A",
                         "BigRIPSIC.fCalMeVSqSum"};
    for(auto &i:tree) i->setloopkeys(keys);

    // Set Parameters according to reaction specified in string
    PID::reactionparameters();

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

    for(auto &t : th) t.join();

    // Calculate the transmission and the crosssection
    double transmission = PID::offctrans();
    PID::crosssection(transmission);

    vector<string> folders{"PID/"+reaction+"/Uncorrected",
                           "PID/"+reaction+"/Corrected",
                           "PID/"+reaction+"/investigate"};
    for (auto &i :folders) output->mkdir(i.c_str());

    for(uint i=0; i<folders.size(); i++){
        output->cd(folders.at(i).c_str());
        if(i<2) for(auto &elem: PIDplot.at(i)) elem.Write();
        else for(auto &elem: reactF5) elem.Write();
    }
    output->cd("");

}

double PID::offctrans() {
    // Calculating off-center transmission
    int threshold = 200; // 100cts for reacted beam
    double dividend = 0; // dividend/divisor
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
    }

    double alpha =dividend/divisor; // Scaling factor F5(beam) F5(reacted)

    // Calculating corresponding chi-square, and test minimum
    for(int i=binlow; i<binhigh; i++)
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
        }

    // Pushing out the new scaled beam values
    for(auto &i: chisq) i /= binhigh - binlow -1;
    reactF5.push_back(reactF5.at(0));
    reactF5.back().Scale(alpha);
    reactF5.back().SetTitle("F5 scaled beam profile");
    // Off-center cut works only for sufficient statistics:
    if(reactF5.at(1).Integral() > threshold){
        printf("\n Scaling Factor: %f, with chisq/ndof %f, transmission: %f. "
            "Scaling factor from Area %f !\n",
            alpha, chisq.at(0),reactF5.at(1).Integral()/reactF5.at(2).Integral(),
            reactF5.at(1).Integral(binlow,binhigh)/
            reactF5.at(0).Integral(binlow,binhigh));
        return reactF5.at(1).Integral()/reactF5.at(2).Integral();
    }
    else{
        printf("Too low statistics for offcenter effects (%.1f)\n",
                reactF5.at(1).Integral());
        return 1;
    }
}

void PID::crosssection(double transmission) {
    // Calculate the final crosssection for this run
    double numberdensity = 0.433; // atoms/cm2 for 10cm LH2
    double tottransmission = 0.7754; // empty target transmission for centered beam
    double tottransmissionerror = 0.01874; // associated error

    double crosssection = 1./numberdensity*reactionpid2/
                          reactionpid1/tottransmission/transmission;
    double cserror = crosssection*
                     pow(1./reactionpid1 + 2./reactionpid2+
                         1./reactF5.at(0).Integral()+
                         pow(tottransmissionerror,2),0.5);  // Error on Transm.

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
            break;
        }
        case runinfo::emptysize:{
            incval = nancyempty::incval;
            targetval = nancyempty::targetval;
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
                acceptancerange = vector<double>{-15,70};
            }
            else __throw_invalid_argument("Invalid reaction !\n");

            break;
        }
    }
}

void PID::histogramsetup() {
    // Nasty histogram setup routine
    vector<TH2D> temp1, temp2;
    temp1.emplace_back(TH2D("pidinc", "PID Incoming F3-F7",  300,2.45,2.9, 200,30,50));
    temp1.emplace_back(TH2D("pidout", "PID Outgoing F8-F11", 300,2.45,2.9, 200,30,50));
    temp1.emplace_back(TH2D("pidincut","PID Inc cut  F8-F11",300,2.45,2.9,200,30,50));
    temp1.emplace_back(TH2D("pidincutout", "PID Out w/ in cut", 300,2.45,2.9,200,30,50));
    temp2.emplace_back(TH2D("pidinccorr", "PID Incoming F3-F7",  300,2.45,2.9, 200,30,50));
    temp2.emplace_back(TH2D("pidoutcorr", "PID Outgoing F8-F11", 300,2.45,2.9, 200,30,50));
    temp2.emplace_back(TH2D("pidincutcorr","PID Inc cut  F8-F11",300,2.45,2.9,200,30,50));
    temp2.emplace_back(TH2D("pidincutoutcorr", "PID Out w/ in cut", 300,2.45,2.9,200,30,50));
    PIDplot.emplace_back(temp1);
    PIDplot.emplace_back(temp2);

    reactF5.emplace_back(TH1D("F5beam", "F5 beam profile", binning,-100,100));
    reactF5.emplace_back(TH1D("F5react", "F5-position of reacted particles", binning,-100,100));

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