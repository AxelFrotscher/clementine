#include "histograms.hh"

using namespace std;

const bool closeness(const vector<double> &d, double sigma = 0.1){
    // Standard deviation for vector d
    double sum = accumulate(d.begin(),d.end(),0.0);
    double mean = sum /d.size();
    vector<double> diff(d.size());

    transform(d.begin(),d.end(), diff.begin(), 
              [mean](double x) { return x - mean; });

    double sq_sum= inner_product(diff.begin(),diff.end(),diff.begin(),0);
    return sigma * mean > sqrt(sq_sum / d.size());
}

void makepid(treereader &tree, TFile &output, const vector<bool> &goodevents){
    // 5 beams, 2incoming, 3outgoing
    vector <string> keys{"BigRIPSBeam.aoq", "BigRIPSBeam.zet"};
    tree.setloopkeys(keys);

    vector<TH2D> PID {
        TH2D("pidinc", "PID Incoming F3-F7",  300,2.45,2.9, 200,25,45),
        TH2D("pidout", "PID Outgoing F8-F11", 300,2.45,2.9, 200,25,45),
        TH2D("pidincut","PID Inc cut  F8-F11",300,2.45,2.9,200,25,45),
        TH2D("pidincutout", "PID Out w/ in cut", 300,2.45,2.9,200,25,45)};

    for(auto &elem: PID){
        elem.SetOption("colz");
        elem.GetXaxis()->SetTitle("A/Q");
        elem.GetYaxis()->SetTitle("Z");
    }

    vector<double> cutval{
        2.575, // center x
        32.5,  // center y
        0.008, // radius x
        0.6    // radius y
    };

    int eventcounter =0;
    vector<vector<double>> valinc;     // Store temporary beam values
    while(tree.singleloop()){
        if(goodevents.at(eventcounter)){
        //Loop over all elements in the tree
        valinc.push_back({tree.BigRIPSBeam_aoq[0],tree.BigRIPSBeam_aoq[1],
                          tree.BigRIPSBeam_aoq[2]});
        valinc.push_back({tree.BigRIPSBeam_zet[0],tree.BigRIPSBeam_zet[1],
                          tree.BigRIPSBeam_zet[2]});
        if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
            // Cut Particles that have variating aoq or zet
            PID.at(0).Fill(tree.BigRIPSBeam_aoq[0],tree.BigRIPSBeam_zet[0]);
        }

        valinc.push_back({tree.BigRIPSBeam_aoq[3],tree.BigRIPSBeam_aoq[4]});
        valinc.push_back({tree.BigRIPSBeam_zet[3], tree.BigRIPSBeam_zet[4]});
        if(closeness(valinc.at(2)) && closeness(valinc.at(3))){
            PID.at(1).Fill(tree.BigRIPSBeam_aoq[4],tree.BigRIPSBeam_zet[4]);
        }
        valinc.clear();

        if((pow(1./cutval.at(2)*(tree.BigRIPSBeam_aoq[0]-cutval.at(0)),2) +
           pow(1/cutval.at(3)*(tree.BigRIPSBeam_zet[0]- cutval.at(1)),2)) <1){
            PID.at(2).Fill(tree.BigRIPSBeam_aoq[0],tree.BigRIPSBeam_zet[0]);
            PID.at(3).Fill(tree.BigRIPSBeam_aoq[4],tree.BigRIPSBeam_zet[4]);
        }}
        eventcounter++;
    }
    output.mkdir("PID");
    output.cd("PID");
    for(auto &elem: PID) elem.Write();
    output.cd("");
}

void plastics(treereader &tree, TFile &output, vector<bool> &goodevents){
    // This function aims to rebuild the trigger from Q1*Q2 F7plastic
    // therefore we set a limit on sqrt(Q1*Q2) to suppress random noise
    printf("Now beginning with reconstruction of plastic scintillators ...\n");
    const int numplastic = 4;
    const int F7pos = 1;
    const int threshhold = 450; // random noise suppresion value

    vector<string> keys{"BigRIPSPlastic.fQLRaw", "BigRIPSPlastic.fQRRaw",
                        "BigRIPSPlastic.fTLRaw", "BigRIPSPlastic.fTRRaw",
                        "BigRIPSPlastic.fpl"};
    tree.setloopkeys(keys);
    // generate output diagrams
    vector<TH2D> qcorr2D;  // Charge Correlation between 1 & 2
    vector<TH1D> qcorr;    // deposited charge at detector
    vector<TH2D> tqcorr2D; // Charge-Time correlation
    
    vector<vector<string>> arrayname, arraytitle;
    tree.singleloop(); // Get first event for invovled focal planes
    for(uint i=0; i<numplastic; i++){
        arrayname.push_back(vector<string>{
            "f7pltrigQ."  + to_string(tree.BigRIPSPlastic_fpl[i]),
            "PlasticQ2D." + to_string(tree.BigRIPSPlastic_fpl[i]),
            "tQcorr."     + to_string(tree.BigRIPSPlastic_fpl[i])
        });

        arraytitle.push_back(vector<string>{
            "Charge deposition by Plastic F"     + to_string(tree.BigRIPSPlastic_fpl[i]),
            "Charge distribution by Plastic F"   + to_string(tree.BigRIPSPlastic_fpl[i]),
            "Charge Ratio/time difference Pl. F" + to_string(tree.BigRIPSPlastic_fpl[i])
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
    vector<vector<int>> range{{400,500},{450,750},{213,350},{260,1310}};
    int i=0; // counting variable

    while(tree.singleloop()){
        if(sqrt(tree.BigRIPSPlastic_fQLRaw[F7pos]*
                tree.BigRIPSPlastic_fQRRaw[F7pos])> threshhold)
        for(uint i=0; i<numplastic; i++){
            if((tree.BigRIPSPlastic_fQLRaw[i] >0 )&&
               (tree.BigRIPSPlastic_fQRRaw[i] >0)){
                qcorr2D.at(i).Fill(tree.BigRIPSPlastic_fQLRaw[i],
                                   tree.BigRIPSPlastic_fQRRaw[i]);
                qcorr.at(i).Fill(sqrt(tree.BigRIPSPlastic_fQLRaw[i]*
                                      tree.BigRIPSPlastic_fQRRaw[i]));
                if((tree.BigRIPSPlastic_fTLRaw[i] >0) &&
                   (tree.BigRIPSPlastic_fTRRaw[i] >0) ){
                    tqcorr2D.at(i).Fill(tree.BigRIPSPlastic_fTLRaw[i]-
                                        tree.BigRIPSPlastic_fTRRaw[i],
                        TMath::Log(tree.BigRIPSPlastic_fQLRaw[i]/
                                   (double)tree.BigRIPSPlastic_fQRRaw[i]));
                }
            }
        }
        for(int j=0;j<numplastic;j++){
            goodevents[i] =goodevents[i]*
                       (sqrt(tree.BigRIPSPlastic_fQLRaw[j]*
                             tree.BigRIPSPlastic_fQRRaw[j]) > range[j][0])*
                       (sqrt(tree.BigRIPSPlastic_fQLRaw[j]*
                             tree.BigRIPSPlastic_fQRRaw[j]) < range[j][1]);
        }
        i++;
    }
    
    printf("Writing Plastic 'Histograms to disk ... \n");
    
    output.mkdir("Plastics/2D");
    output.mkdir("Plastics/Q1Q2");
    output.mkdir("Plastics/TQCorr");
    output.cd("Plastics/2D");
    for(auto histo: qcorr2D) histo.Write();
    output.cd("Plastics/Q1Q2");
    for(auto histo: qcorr) histo.Write();
    output.cd("Plastics/TQCorr");
    for(auto histo: tqcorr2D) histo.Write();
    output.cd("");

    // A careful analysis yields a threshhold of sqrt(Q1*Q2)>160
    printf("Finished Writing plastic histograms! \n");
}

void ppacs(treereader &tree, TFile &output, vector<bool> &goodevents){
    printf("Now beginning with reconstruction of the ppac's ...\n");
    const int numplane = 36;
    const int pl11position = 3; //Plastic at F11 is fourth in array
    vector<string> keys{"BigRIPSPlastic.fTime", "BigRIPSPPAC.fFiredX",
                        "BigRIPSPPAC.fFiredY", "BigRIPSPPAC.fTSumX",
                        "BigRIPSPPAC.fTDiffX", "BigRIPSPPAC.fTSumY",
                        "BigRIPSPPAC.fTDiffY", "BigRIPSPPAC.name"};

    tree.setloopkeys(keys);

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
               
    sumdiffppac.at(0).at(0).SetOption("colz");
    sumdiffppac.at(0).at(0).GetXaxis()->SetTitle("PPAC [ch]");
    sumdiffppac.at(0).at(0).GetYaxis()->SetTitle("Sum [ns]");
    sumdiffppac.at(0).at(1).SetOption("colz");
    sumdiffppac.at(0).at(1).GetXaxis()->SetTitle("PPAC [ch]");
    sumdiffppac.at(0).at(1).GetYaxis()->SetTitle("diff [ns]");
    sumdiffppac.at(1).at(0).SetOption("colz");
    sumdiffppac.at(1).at(0).GetXaxis()->SetTitle("PPAC [ch]");
    sumdiffppac.at(1).at(0).GetYaxis()->SetTitle("Sum [ns]");
    sumdiffppac.at(1).at(1).SetOption("colz");
    sumdiffppac.at(1).at(1).GetXaxis()->SetTitle("PPAC [ch]");
    sumdiffppac.at(1).at(1).GetYaxis()->SetTitle("diff [ns]");
    
    effPPACX.GetXaxis()->SetTitle("PPAC [ch]");
    effPPACX.GetYaxis()->SetTitle("PPACX(x)/Plastic(F11)");
    effPPACY.GetXaxis()->SetTitle("PPAC [ch]");
    effPPACY.GetYaxis()->SetTitle("PPACY(x)/Pplastic(F11)");

    //Fill names with first event
    tree.singleloop();
    //cout << "T STRING: " << bripsname.operator[](0).Data() << endl;
    for(int i=1; i<=numplane; i++){
        effPPACX.GetXaxis()->SetBinLabel(i, tree.BigRIPSPPAC_name[i-1].Data());
        effPPACY.GetXaxis()->SetBinLabel(i, tree.BigRIPSPPAC_name[i-1].Data());
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

    printf("Filling in Plastics ... \n");
    while(tree.singleloop()){
        // Get Efficiency relative to the last Plastic
        if (tree.BigRIPSPlastic_fTime[pl11position] >0){
            for(int i=1; i<=numplane; i++){
                //ripsvec.at(0).at(i) 
                if(tree.BigRIPSPPAC_fFiredX[i-1])
                    effPPACX.SetBinContent(i, effPPACX.GetBinContent(i)+1);
                if(tree.BigRIPSPPAC_fFiredY[i-1])
                    effPPACY.SetBinContent(i, effPPACY.GetBinContent(i)+1);
            }
            total++;
        }
        for(int i=0; i<numplane;i++){
            sumdiffppac.at(0).at(0).Fill(i,tree.BigRIPSPPAC_fTSumX[i]);
            sumdiffppac.at(0).at(1).Fill(i,tree.BigRIPSPPAC_fTDiffX[i]);
            sumdiffppac.at(1).at(0).Fill(i,tree.BigRIPSPPAC_fTSumY[i]);
            sumdiffppac.at(1).at(1).Fill(i,tree.BigRIPSPPAC_fTDiffY[i]);
        }

        // Make Cut conditions: each focal Plane needs to have one valid entry
        // := Sum for x and y is in range
        for(uint i =0; i<ppacplane.size();i++){
            for(uint j=0; j<4; j++){ // Loop over all 4 PPAC's per Focal Plane
            temptruth.at(j) =
                (tree.BigRIPSPPAC_fTSumX[ppacplane.at(i).at(j)] > ppacrange.at(i).at(j).at(0) &&
                 tree.BigRIPSPPAC_fTSumX[ppacplane.at(i).at(j)] < ppacrange.at(i).at(j).at(1) &&
                 tree.BigRIPSPPAC_fTSumY[ppacplane.at(i).at(j)] > ppacrange.at(i).at(j).at(2) &&
                 tree.BigRIPSPPAC_fTSumY[ppacplane.at(i).at(j)] < ppacrange.at(i).at(j).at(3));
            }
            // Recursively check each focal plane for at least 1 good Signal
            goodevents.at(fulltotal) = goodevents.at(fulltotal)*
                                       accumulate(temptruth.begin(),
                                                  temptruth.end(),0);
        }
        fulltotal++;
    }

    effPPACX.Scale(1/(double)total);
    effPPACY.Scale(1/(double)total);

    printf("Writing PPAC histograms... \n");
    output.mkdir("PPAC/FiredEff");
    output.cd("PPAC/FiredEff");
    effPPACX.Write();
    effPPACY.Write();

    vector<string> xy = {"X","Y"};
    vector<string> sd = {"Sum", "Diff"};
    
    // Generating output structure
    for(uint j=0; j<2; j++){ // Loop over X and Y//
        string folder = "PPAC/" + sd.at(j);
        output.mkdir(folder.c_str());
        output.cd(folder.c_str());
        for(uint k=0; k<2; k++){     // Loop Sum and diff
            sumdiffppac.at(k).at(j).Write();
        }
    }

    output.cd("");
    printf("Finished Writing PPAC histogram!\n");
}

void ionisationchamber(treereader &alt2dtree, TFile &output,
                       vector<bool> &goodevents) {
    // This method aims to control the Ionisation Chamber values
    // therefore only certain events are accepted

    printf("Now beginning analysis of the IC's (Fpl 7 & 11)\n");
    //datree.Restart();

    // Explicit Method for 2D arrays
    // Create new reader class

    vector <string> readoutkeys{"BigRIPSIC.nhitchannel", "BigRIPSIC.fADC[32]"};
    alt2dtree.setloopkeys(readoutkeys);

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
    uint totalcounter =0;

    while(alt2dtree.singleloop()){
        if((alt2dtree.BigRIPSIC_nhitchannel[0]*
           alt2dtree.BigRIPSIC_nhitchannel[1]) <16){ //Require 4 hits in each IC
            goodevents.at(totalcounter) = false;
        }
        if((alt2dtree.BigRIPSIC_fADC[0][0] <0) ||
           (alt2dtree.BigRIPSIC_fADC[1][0] <0)) continue;
        for(int i=0; i<(numchannel-1);i++){
            if(alt2dtree.BigRIPSIC_fADC[0][i] >0)
                comparediag.at(i).at(0).Fill(
                        alt2dtree.BigRIPSIC_fADC[0][0],
                        alt2dtree.BigRIPSIC_fADC[0][i+1]);
            if(alt2dtree.BigRIPSIC_fADC[1][i] >0)
                comparediag.at(i).at(1).Fill(
                        alt2dtree.BigRIPSIC_fADC[1][0],
                        alt2dtree.BigRIPSIC_fADC[1][i+1]);
        }
        totalcounter++;
    }

    output.mkdir("IC/IC7");
    output.mkdir("IC/IC11");
    output.cd("IC/IC7");
    for(auto &elem: comparediag) elem.at(0).Write();
    output.cd("IC/IC11");
    for(auto &elem: comparediag) elem.at(1).Write();
    output.cd("");
}

void highordercorrection(treereader &alt2dtree, TFile &output){
    // This method aims to determine the higher-order corrections
    // for the matrix elements
    printf("Now beginning with higher order corrections ...\n");
    vector<string> keys{"F3X","F3A", "F5X", "F5A", "BigRIPSBeam.aoq"};
    alt2dtree.setloopkeys(keys);
    
    const int beam = 0; // Evaluate Beam F3-7 (1st element)
    const vector<string> arrname = { "DepF3X","DepF3A","DepF5X", "DepF5A"};
    const vector<string> arrtitle = {
        "Dependence of F3X vs AoQ", "Dependence of F3A vs AoQ", 
        "Dependence of F5X vs AoQ", "Dependence of F5A vs AoQ"};
    
    vector<TH2D> culpritdiag;

    for(uint i=0; i<arrname.size(); i++){
      culpritdiag.emplace_back(TH2D(arrname.at(i).c_str(),arrtitle.at(i).c_str(),
                                    2000,2.2,3.2,500,-100,150));
      culpritdiag.back().SetOption("colz");
      culpritdiag.back().GetXaxis()->SetTitle("A/Q");
      culpritdiag.back().GetYaxis()->SetTitle(arrname.at(i).c_str());
    }

    while(alt2dtree.singleloop()){
        culpritdiag.at(0).Fill(alt2dtree.BigRIPSBeam_aoq[beam], alt2dtree.F3X);
        culpritdiag.at(1).Fill(alt2dtree.BigRIPSBeam_aoq[beam], alt2dtree.F3A);
        culpritdiag.at(2).Fill(alt2dtree.BigRIPSBeam_aoq[beam], alt2dtree.F5X);
        culpritdiag.at(3).Fill(alt2dtree.BigRIPSBeam_aoq[beam], alt2dtree.F5A);
    }
    //gDirec
    output.Delete("Corrections");
    output.mkdir("Corrections");
    output.cd("Corrections");
    for(auto &elem: culpritdiag) elem.Write();
    output.cd("");
    printf("Finished with higher order corrections!\n");
}

void makehistograms(const string input){
    TFile inputfile(input.c_str());
    if(!inputfile.IsOpen()) __throw_invalid_argument("Input file not valid.");

    treereader alt2dtree; // Opens the input file...

    const string output = "build/output/histout.root";
    TFile outputfile(output.c_str(), "RECREATE");
    if(!outputfile.IsOpen()) __throw_invalid_argument("Output file not valid");

    vector<bool> options{
        true,  // Plastics
        true,  // ppac
        true,  // Ionisationchamber
        true,  // highordercorrection
        true   // makepid
    };

    // Store events that cannot be used
    vector<bool> goodevents(1000000, true);

    // Then we read in the tree (lazy)
    TTreeReader mytreereader("tree", &inputfile);

    // Rebuild F7 Trigger combined charge threshhold
    if(options.at(0)) plastics(alt2dtree, outputfile, goodevents);
    
    // Cut the PPAC's
    if(options.at(1)) ppacs(alt2dtree, outputfile, goodevents);

    if(options.at(2)) ionisationchamber(alt2dtree, outputfile, goodevents);
    printf("Finished with PPAC Consistency checks\n");
    // Get Corrections
    if(options.at(3)) highordercorrection(alt2dtree, outputfile);

    // Get Z vs. A/Q
    if(options.at(4)) makepid(alt2dtree, outputfile, goodevents);
    printf("Made PID histograms in %s\n", output.c_str());

    cout << "Runs has " <<accumulate(goodevents.begin(),goodevents.end(),0)
         << " good Elements" << endl;

    //outputfile.Close();
}
