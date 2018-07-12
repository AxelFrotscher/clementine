//
// Created by axel on 05.06.18.
//

#include "histogram_cuts.hh"
#include "MakeAllTree_78Ni.hh"
#include "libconstant.h"

using namespace std;

const bool closeness(const vector<double> &d, double sigma){
    // Standard deviation for vector d
    double sum = accumulate(d.begin(),d.end(),0.0);
    double mean = sum /d.size();
    vector<double> diff(d.size());

    transform(d.begin(),d.end(), diff.begin(),
              [mean](double x) { return x - mean; });

    double sq_sum= inner_product(diff.begin(),diff.end(),diff.begin(),0);
    return sigma * mean > sqrt(sq_sum / d.size());
}

double linfit(double *x, double *par){
    // Linear fit function
    return par[0] + par[1]*x[0];
}

void plastics(treereader *tree, TFile *output, vector<bool> &goodevents){
    // This function aims to rebuild the trigger from Q1*Q2 F7plastic
    // therefore we set a limit on sqrt(Q1*Q2) to suppress random noise
    printf("Now beginning with reconstruction of plastic scintillators ...\n");
    const int numplastic = 4;
    const int F7pos = 1;
    const int threshhold = 450; // random noise suppresion value

    vector<string> keys{"BigRIPSPlastic.fQLRaw", "BigRIPSPlastic.fQRRaw",
                        "BigRIPSPlastic.fTLRaw", "BigRIPSPlastic.fTRRaw",
                        "BigRIPSPlastic.fpl"};
    tree->setloopkeys(keys);
    // generate output diagrams
    vector<TH2D> qcorr2D;  // Charge Correlation between 1 & 2
    vector<TH1D> qcorr;    // deposited charge at detector
    vector<TH2D> tqcorr2D; // Charge-Time correlation

    vector<vector<string>> arrayname, arraytitle;
    tree->singleloop(); // Get first event for involved focal planes
    for(uint i=0; i<numplastic; i++){
        arrayname.push_back(vector<string>{
            "f7pltrigQ."  + to_string(tree->BigRIPSPlastic_fpl[i]),
            "PlasticQ2D." + to_string(tree->BigRIPSPlastic_fpl[i]),
            "tQcorr."     + to_string(tree->BigRIPSPlastic_fpl[i])
        });

        arraytitle.push_back(vector<string>{
            "Charge deposition by Plastic F"     + to_string(tree->BigRIPSPlastic_fpl[i]),
            "Charge distribution by Plastic F"   + to_string(tree->BigRIPSPlastic_fpl[i]),
            "Charge Ratio/time difference Pl. F" + to_string(tree->BigRIPSPlastic_fpl[i])
        });

        qcorr2D.emplace_back(TH2D(arrayname.at(i).at(1).c_str(),
                                  arraytitle.at(i).at(1).c_str(),
                                  1500,0,1500, 1500,0,1500));
        qcorr.emplace_back(TH1D(arrayname.at(i).at(0).c_str(),
                                arraytitle.at(i).at(0).c_str(), 1500,0,1500));
        tqcorr2D.emplace_back(TH2D(arrayname.at(i).at(2).c_str(),
                                   arraytitle.at(i).at(2).c_str(),
                                   300,-150,150, 400,-2,2));

        qcorr.back().GetXaxis()->SetTitle("#sqrt{Q_{l}Q_{r}} [ch]");
        qcorr.back().GetYaxis()->SetTitle("N");
        qcorr2D.back().GetXaxis()->SetTitle("Q_{l} [ch]");
        qcorr2D.back().GetYaxis()->SetTitle("Q_{r} [ch]");
        tqcorr2D.back().GetXaxis()->SetTitle("t_{l}-t_{r} [ns]");
        tqcorr2D.back().GetYaxis()->SetTitle("ln(Q_{l}/Q_{r})");

        //Cosmetic changes
        tqcorr2D.back().SetOption("colz");
        qcorr2D.back().SetOption("colz");
    }

    // To avoid multihit-triggers we define a range of acceptance
    vector<vector<int>> range{{480,620},{700,920},{220,330},{270,1510}};

    // Progress Bar setup
    int eventno=0; // counting variable
    uint totevents = goodevents.size();
    const int downscale = 50000; // every n-th event
    // Get Number of Good events before
    int cutcount =0;

    while(tree->singleloop()){
        if(goodevents.at(eventno)){
        if(sqrt(tree->BigRIPSPlastic_fQLRaw[F7pos]*
                tree->BigRIPSPlastic_fQRRaw[F7pos])> threshhold) // F7-cut
            for(uint i=0; i<numplastic; i++){                   // plastic loop
                if((tree->BigRIPSPlastic_fQLRaw[i] >0 )&&
                   (tree->BigRIPSPlastic_fQRRaw[i] >0)){   // 0 charge veto
                    qcorr2D.at(i).Fill(tree->BigRIPSPlastic_fQLRaw[i],
                                       tree->BigRIPSPlastic_fQRRaw[i]);
                    qcorr.at(i).Fill(sqrt(tree->BigRIPSPlastic_fQLRaw[i]*
                                          tree->BigRIPSPlastic_fQRRaw[i]));
                    if((tree->BigRIPSPlastic_fTLRaw[i] >0) &&
                       (tree->BigRIPSPlastic_fTRRaw[i] >0) ){ // 0 time veto
                        tqcorr2D.at(i).Fill(tree->BigRIPSPlastic_fTLRaw[i]-
                                            tree->BigRIPSPlastic_fTRRaw[i],
                                            TMath::Log(tree->BigRIPSPlastic_fQLRaw[i]/
                                            (double)tree->BigRIPSPlastic_fQRRaw[i]));
                    }
                }
            }
            goodeventmutex.lock();
            for(int j=0;j<numplastic;j++){
                goodevents[eventno] =goodevents[eventno]*
                               (sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                                     tree->BigRIPSPlastic_fQRRaw[j]) > range[j][0])*
                               (sqrt(tree->BigRIPSPlastic_fQLRaw[j]*
                                     tree->BigRIPSPlastic_fQRRaw[j]) < range[j][1]);
            }
            cutcount = cutcount + 1-goodevents[eventno];
            goodeventmutex.unlock();
        }
        eventno++;
        if(!(eventno%downscale)){
            consolemutex.lock();
            progressbar(eventno,totevents, 1);
            consolemutex.unlock();
        }
    }
    writemutex.lock();
    printf("\nPlastic Cut out %i Events %f %%\n", cutcount,
           100*cutcount/(double)totevents);

    vector<string> folders{"Plastics/2D", "Plastics/Q1Q2", "Plastics/TQCorr"};
    for (auto str:folders) output->mkdir(str.c_str());

    output->cd("Plastics/2D");
    for(auto histo: qcorr2D) histo.Write();
    output->cd("Plastics/Q1Q2");
    for(auto histo: qcorr) histo.Write();
    output->cd("Plastics/TQCorr");
    for(auto histo: tqcorr2D) histo.Write();
    output->cd("");
    writemutex.unlock();
    printf("Finished Writing plastic histograms! \n");
}

void ppacs(treereader *tree, TFile *output, vector<bool> &goodevents){
    printf("Now beginning with reconstruction of the ppac's ...\n");
    const int numplane = 36;
    const int pl11position = 3; //Plastic at F11 is fourth in array
    vector<string> keys{"BigRIPSPlastic.fTime", "BigRIPSPPAC.fFiredX",
                        "BigRIPSPPAC.fFiredY", "BigRIPSPPAC.fTSumX",
                        "BigRIPSPPAC.fTDiffX", "BigRIPSPPAC.fTSumY",
                        "BigRIPSPPAC.fTDiffY", "BigRIPSPPAC.name"};

    tree->setloopkeys(keys);

    //36 Values per Array (Event)
    TH1D effPPACX("effPPACX", "Efficiency of PPAC X", numplane, 0, numplane);
    TH1D effPPACY("effPPACY", "Efficiency of PPAC Y", numplane, 0, numplane);

    vector<string> arrname = { "PPACXsum","PPACXdiff","PPACYsum", "PPACYdiff"};
    vector<string> arrtitle = {
            "Sum of Signals PPACX", "Difference of Signals PPACX",
            "Sum of Signals PPACY", "Difference of Signals PPACY"};
    vector<vector<TH2D>> sumdiffppac // 1 X,Y 2 Sum,diff (2D NoPPAC, Quantity)
            {{TH2D(arrname.at(0).c_str(),arrtitle.at(0).c_str(),
                   numplane,0,numplane,300,000,300),
              TH2D(arrname.at(1).c_str(),arrtitle.at(1).c_str(),
                   numplane,0,numplane,800,-200,200)},
             {TH2D(arrname.at(2).c_str(),arrtitle.at(2).c_str(),
                   numplane,0,numplane,150,0,150),
              TH2D(arrname.at(3).c_str(),arrtitle.at(3).c_str(),
                   numplane,0,numplane,800,-200,200)}};

    for(auto &hist : sumdiffppac){
        uint toggle =0;
        for(auto &h: hist){
            h.SetOption("colz");
            h.GetXaxis()->SetTitle("PPAC [ch]");
            if(!toggle) h.GetYaxis()->SetTitle("Sum [ns]");
            else h.GetYaxis()->SetTitle("diff [ns]");
            toggle++;
        }
    }

    effPPACX.GetXaxis()->SetTitle("PPAC [ch]");
    effPPACX.GetYaxis()->SetTitle("PPACX(x)/Plastic(F11)");
    effPPACY.GetXaxis()->SetTitle("PPAC [ch]");
    effPPACY.GetYaxis()->SetTitle("PPACY(x)/Pplastic(F11)");

    //Fill names with first event
    tree->singleloop();
    //cout << "T STRING: " << bripsname.operator[](0).Data() << endl;
    for(int i=1; i<=numplane; i++){
        effPPACX.GetXaxis()->SetBinLabel(i, tree->BigRIPSPPAC_name[i-1].Data());
        effPPACY.GetXaxis()->SetBinLabel(i, tree->BigRIPSPPAC_name[i-1].Data());
    }

    // 0=X, 1=Y  || 0...35 Plane No.
    uint total = 0;
    uint fulltotal =1;  //Skipping first event for name generation

    // Define Planes of Importance corresponding to fpl 3,5,7,8,9,11
    const vector<vector<int>> ppacplane{{4,5,6,7},{9,10,11,12},{14,15,16,17},
                                        {18,19,20,21},{22,23,24,25},{30,31,32,33}};

    const vector<vector<vector<int>>> ppacrange{ //{xlow,xup,ylow,yup}
        {{85,105,85,100},{85,100,85,100},{165,180,85,105},{170,185,90,105}},  //Fpl 3
        {{165,185,85,105},{165,185,90,105},{165,175,90,105},{165,180,95,110}},//Fpl 5
        {{160,175,95,110},{165,180,85,100},{75,90,90,100},{85,105,95,110}},   //Fpl 7
        {{170,185,99,110},{165,185,95,110},{165,185,95,110},{165,185,90,105}},//Fpl 8
        {{130,145,60,75},{135,150,60,75},{140,155,70,80},{130,145,50,65}},//Fpl 9
        {{145,160,65,80},{135,150,55,70},{145,160,75,90},{135,150,60,75}} //Fpl11
    };

    vector<bool> temptruth(4, true);

    // Progress Bar setup
    Long64_t totevents = goodevents.size();
    const int downscale = 50000; // every n-th event
    int cutno = 0;

    while(tree->singleloop()){
        if(goodevents.at(fulltotal)){ // Determine cut only on good events
            // Get Efficiency relative to the last Plastic
            if (tree->BigRIPSPlastic_fTime[pl11position] >0){
                for(int i=1; i<=numplane; i++){
                    if(tree->BigRIPSPPAC_fFiredX[i-1])
                        effPPACX.SetBinContent(i, effPPACX.GetBinContent(i)+1);
                    if(tree->BigRIPSPPAC_fFiredY[i-1])
                        effPPACY.SetBinContent(i, effPPACY.GetBinContent(i)+1);
                }
                total++;
            }
            for(int i=0; i<numplane;i++){
                sumdiffppac.at(0).at(0).Fill(i,tree->BigRIPSPPAC_fTSumX[i]);
                sumdiffppac.at(0).at(1).Fill(i,tree->BigRIPSPPAC_fTDiffX[i]);
                sumdiffppac.at(1).at(0).Fill(i,tree->BigRIPSPPAC_fTSumY[i]);
                sumdiffppac.at(1).at(1).Fill(i,tree->BigRIPSPPAC_fTDiffY[i]);
            }

            // Make Cut conditions: each focal Plane needs to have one valid entry
            // := Sum for x and y is in range
            bool temp = true;

            for(uint i =0; i<ppacplane.size();i++){
                for(uint j=0; j<4; j++){ // Loop over all 4 PPAC's per Focal Plane
                    int cpl = ppacplane.at(i).at(j);
                    temptruth.at(j) =
                       (tree->BigRIPSPPAC_fTSumX[cpl] > ppacrange.at(i).at(j).at(0) &&
                        tree->BigRIPSPPAC_fTSumX[cpl] < ppacrange.at(i).at(j).at(1) &&
                        tree->BigRIPSPPAC_fTSumY[cpl] > ppacrange.at(i).at(j).at(2) &&
                        tree->BigRIPSPPAC_fTSumY[cpl] < ppacrange.at(i).at(j).at(3));
                    if(temptruth.at(j)) break; // We got our nice event at this F-point
                }
                // Recursively check each focal plane for at least 1 good Signal
                temp = temp*accumulate(temptruth.begin(), temptruth.end(),0);
            }
            if(!temp){
                goodeventmutex.lock();
                goodevents.at(fulltotal) = false;
                cutno++;
                goodeventmutex.unlock();
            }
        }
        fulltotal++;
        if(!(fulltotal%downscale)){
            consolemutex.lock();
            progressbar(fulltotal, totevents,3);
            consolemutex.unlock();
        }
    }

    effPPACX.Scale(1/(double)total);
    effPPACY.Scale(1/(double)total);

    vector<string> xy = {"X","Y"};
    vector<string> sd = {"Sum", "Diff"};

    writemutex.lock();

    printf("\nPPAC Cut out %i Events %f %%\n", cutno,
           100*cutno/(double)totevents);

    output->mkdir("PPAC/FiredEff");
    output->cd("PPAC/FiredEff");
    effPPACX.Write();
    effPPACY.Write();

    // Generating output structure
    for(uint j=0; j<2; j++){ // Loop over X and Y//
        string folder = "PPAC/" + sd.at(j);
        output->mkdir(folder.c_str());
        output->cd(folder.c_str());
        for(uint k=0; k<2; k++){     // Loop Sum and diff
            sumdiffppac.at(k).at(j).Write();
        }
    }

    output->cd("");
    printf("Finished Writing PPAC histogram!\n");
    writemutex.unlock();
}

void ionisationchamber(treereader *alt2dtree, TFile *output,
                       vector<bool> &goodevents) {
    // This method aims to control the Ionisation Chamber values
    // therefore only certain events are accepted

    printf("Now beginning analysis of the IC's (Fpl 7 & 11)\n");
    //datree.Restart();

    // Explicit Method for 2D arrays

    vector <string> readoutkeys{"BigRIPSIC.nhitchannel", "BigRIPSIC.fADC[32]"};
    alt2dtree->setloopkeys(readoutkeys);

    const int numchannel = 6; // There are 6 channel per IC
    const int numic = 2; // Number of Ionisation Chambers
    vector<vector<TH2D>> comparediag;

    for(int i=1; i<numchannel; i++){
        vector<string> arrname = {"ICratio" + to_string(i) + "to0fpl7",
                                  "ICratio" + to_string(i) + "to0fpl11"};
        string arrtitle = "Peak ADC Ratio " + to_string(i) + " to 0";
        comparediag.push_back({TH2D(arrname.at(0).c_str(),arrtitle.c_str(),
                                    2048,0,16384,2048,0,16384),
                               TH2D(arrname.at(1).c_str(),arrtitle.c_str(),
                                    2048,0,16384,2048,0,16384)});
        for(auto &elem: comparediag.back()){
            elem.SetOption("colz");
            elem.GetXaxis()->SetTitle("ADC0 [ch]");
            elem.GetYaxis()->SetTitle(("ADC"+to_string(i)+ "[ch]").c_str());
        }
    }

    // Progress Bar setup
    uint totalcounter =0; // counting variable
    const Long64_t totevents = goodevents.size();
    uint cutcount =0;
    const int downscale = 50000; // every n-th event

    while(alt2dtree->singleloop()){
        if(goodevents.at(totalcounter)){ // Determine cut only on good events
            if((alt2dtree->BigRIPSIC_nhitchannel[0]*
                alt2dtree->BigRIPSIC_nhitchannel[1]) <16){//Require 4 hits per IC
                goodeventmutex.lock();
                goodevents.at(totalcounter) = false;
                goodeventmutex.unlock();
                cutcount++;
            }
            if((alt2dtree->BigRIPSIC_fADC[0][0] <0) ||
                (alt2dtree->BigRIPSIC_fADC[1][0] <0)) continue;
            for(int i=0; i<(numchannel-1);i++){
                if(alt2dtree->BigRIPSIC_fADC[0][i] >0)
                    comparediag.at(i).at(0).Fill(
                            alt2dtree->BigRIPSIC_fADC[0][0],
                            alt2dtree->BigRIPSIC_fADC[0][i+1]);
                if(alt2dtree->BigRIPSIC_fADC[1][i] >0)
                    comparediag.at(i).at(1).Fill(
                            alt2dtree->BigRIPSIC_fADC[1][0],
                            alt2dtree->BigRIPSIC_fADC[1][i+1]);
            }
        }
        totalcounter++;
        if(!(totalcounter%downscale)){
            consolemutex.lock();
            progressbar(totalcounter,totevents,2);
            consolemutex.unlock();
        }
    }
    writemutex.lock();
    printf("\nIC Cut out %i Events %f %%\n", cutcount,
           100*cutcount/(double)totevents);

    output->mkdir("IC/IC7");
    output->mkdir("IC/IC11");
    output->cd("IC/IC7");
    for(auto &elem: comparediag) elem.at(0).Write();
    output->cd("IC/IC11");
    for(auto &elem: comparediag) elem.at(1).Write();
    output->cd("");
    writemutex.unlock();
}

void chargestatecut(treereader *tree, TFile *output, vector<bool> &goodevents){

    // this method aims to cut charge state changes between F8-9 and F9-11
    printf("Now performing charge state change cut between F8-9 and F9-11...\n");

    // Generate output histogram
    vector<TH2D> cschist{
        TH2D("csc", "Charged state change", 1000,0.8,1.2,1000,3,7),
        TH2D("csccut", "Charged state change cut", 500,0.95,1.05,1000,3,7)};

    for(auto &histo: cschist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("B#rho F8-9");
        histo.GetXaxis()->SetTitle("B#rho F8-9 / B#rho F9-11");
    }

    // Get relevant keys
    vector<string> keys{"BigRIPSBeam.brho"};
    tree->setloopkeys(keys);

    // Progress Bar setup
    uint eventno=0; // counting variable
    uint totevents = goodevents.size();
    const int downscale = 50000; // every n-th event
    // Get Number of Good events before
    int  cutcount =0;
    double brhoratio = 0;

    // Define cut
    TFile cutfile("config/cut.root");
    if(!(cutfile.IsOpen())) __throw_invalid_argument("Could not open cut file "
                                                     "at config/rhocut.root");
    TCutG* mycut;
    if(runinfo::emptysize == goodevents.size()) mycut = (TCutG*)cutfile.Get("brhoempty");
    else mycut = (TCutG*)cutfile.Get("brhocut");
    if(!mycut) __throw_invalid_argument("Could not load cut from file!\n");

    while(tree->singleloop()){
        if(goodevents.at(eventno)){ // Determine cut only on good events
            brhoratio = tree->BigRIPSBeam_brho[3]/tree->BigRIPSBeam_brho[2];
            cschist.at(0).Fill(brhoratio, tree->BigRIPSBeam_brho[2]);

            if(mycut->IsInside(brhoratio,tree->BigRIPSBeam_brho[2])){
                cschist.at(1).Fill(brhoratio,tree->BigRIPSBeam_brho[2]);
            }
            else {
                goodeventmutex.lock();
                goodevents.at(eventno) = false;
                goodeventmutex.unlock();
                cutcount++;
            }
        }
        eventno++;
        if(!(eventno%downscale)){
            consolemutex.lock();
            progressbar(eventno, totevents,0);
            consolemutex.unlock();
        }
    }

    writemutex.lock();
    printf("\nCCSC Cut out %i Events %f %%\n", cutcount,
           100*cutcount/(double)totevents);

    output->mkdir("CSC");
    output->cd("CSC");
    for(auto hist: cschist) hist.Write();
    output->cd("");
    writemutex.unlock();
}

void targetcut(treereader *tree, TFile *output, vector<bool> &goodevents){
    // This method aims to reconstruct the target impinging position at F8
    printf("Now performing the Target cut... \n");

    vector<TH2D> tarhist{
        TH2D("target", "Beam Profile F8", 1000,-50,50,1000,-50,50),
        TH2D("targetcut", "Cut Beam Profile F8", 1000,-50,50,1000,-50,50)};


    for(auto &histo: tarhist){
        histo.SetOption("colz");
        histo.GetYaxis()->SetTitle("x [mm]");
        histo.GetXaxis()->SetTitle("y [mm]");
    }

    // Get relevant keys
    vector<string> keys{"BigRIPSPPAC.fX", "BigRIPSPPAC.fY"};
    tree->setloopkeys(keys);

    const int ppacoffset = 18; // [18] is F8PPAC-1A
    const int f8ppac = 4;
    const int  magnum = -9999;
    const int f8offset = -2; // mm

    // F8 PPAC positions relative to F8 0:PPAcno 1: z(x),z(y)
    const vector<vector<double>> f8z{{-1310.7,-1302.1},{-1273.3,-1281.9},
                                     {-810.7,-802.1},{-773.3,-781.9}};

    vector<vector<double>> slopex(2,vector<double>(0)); //0:x,y 1:val, z
    vector<vector<double>> slopey(2,vector<double>(0)); //0:x,y 1:val, z

    // Progress Bar setup
    uint eventno=0; // counting variable
    uint totevents = goodevents.size();
    const int downscale = 50000; // every n-th event
    // Get Number of Good events before
    int  cutcount =0;
    // Get temporary mean value:
    double meanx = 0, meany = 0, meanxz=0, meanyz =0, fXTar=0, fYTar=0;

    while(tree->singleloop()){
        if(goodevents.at(eventno)){ // Determine cut only on good events
            // 1. Fill valid events
            for(uint i =ppacoffset; i<(ppacoffset+f8ppac); i++) {
                if (abs(tree->BigRIPSPPAC_fX[i] - magnum) > 1){
                    slopex.at(0).push_back(tree->BigRIPSPPAC_fX[i]);
                    slopex.at(1).push_back(f8z.at(i-ppacoffset).at(0));
                }
                if (abs(tree->BigRIPSPPAC_fY[i] - magnum) > 1){
                    slopey.at(0).push_back(tree->BigRIPSPPAC_fY[i]);
                    slopey.at(1).push_back(f8z.at(i-ppacoffset).at(1));
                }
            }
            if((slopex.at(0).size() < 2) || (slopey.at(0).size() <2)){
                // cutting
                cutcount++;
                goodeventmutex.lock();
                goodevents.at(eventno) = false;
                goodeventmutex.unlock();
            }
            else{
                meanx = accumulate(slopex.at(0).begin(), slopex.at(0).end(), 0.0)/
                        slopex.at(0).size();
                meany = accumulate(slopey.at(0).begin(), slopey.at(0).end(), 0.0)/
                        slopey.at(0).size();
                meanxz= accumulate(slopex.at(1).begin(), slopex.at(1).end(), 0.0)/
                        slopex.at(1).size();
                meanyz= accumulate(slopey.at(1).begin(), slopey.at(1).end(), 0.0)/
                        slopey.at(1).size();
                fXTar= meanx - (meanxz-f8offset)*slope(slopex.at(1),slopex.at(0));
                fYTar= meany - (meanyz-f8offset)*slope(slopey.at(1),slopey.at(0));
                tarhist.at(0).Fill(fXTar,fYTar);
                if(pow(fXTar*fXTar+fYTar*fYTar,0.5) < 20){ // target cut
                    tarhist.at(1).Fill(fXTar,fYTar);
                }
                else{
                    cutcount++;
                    goodeventmutex.lock();
                    goodevents.at(eventno) = false;
                    goodeventmutex.unlock();
                }
            }

            slopex.at(0).clear();
            slopey.at(0).clear();
            slopex.at(1).clear();
            slopey.at(1).clear();
        }

        eventno++;
        if(!(eventno%downscale)){
            consolemutex.lock();
            progressbar(eventno, totevents,5);
            consolemutex.unlock();
        }
    }

    writemutex.lock();
    printf("\nTarget Cut out %i Events %f %%\n", cutcount,
           100*cutcount/(double)totevents);

    output->mkdir("Target");
    output->cd("Target");
    for(auto hist: tarhist) hist.Write();
    output->cd("");
    writemutex.unlock();

}
