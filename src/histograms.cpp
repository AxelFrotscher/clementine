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

void makepid(TTreeReader &datree, TFile &output, const vector<bool> &goodevents){
    datree.Restart();
    // 5 beams, 2incoming, 3outgoing
    TTreeReaderArray<double> aoqevt (datree, "BigRIPSBeam.aoq");
    TTreeReaderArray<double> zetevt (datree, "BigRIPSBeam.zet");

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
    while(datree.Next()){
        if(goodevents.at(eventcounter)){
        //Loop over all elements in the tree
        valinc.push_back({aoqevt[0],aoqevt[1],aoqevt[2]});
        valinc.push_back({zetevt[0],zetevt[1],zetevt[2]});
        if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
            // Cut Particles that have variating aoq or zet
            PID.at(0).Fill(aoqevt[0],zetevt[0]);
        }

        valinc.push_back({aoqevt[3],aoqevt[4]});
        valinc.push_back({zetevt[3], zetevt[4]});
        if(closeness(valinc.at(2)) && closeness(valinc.at(3))){
            PID.at(1).Fill(aoqevt[4],zetevt[4]);
        }
        valinc.clear();

        if((pow(1./cutval.at(2)*(aoqevt[0]-cutval.at(0)),2) +
           pow(1/cutval.at(3)*(zetevt[0]- cutval.at(1)),2)) <1){
            PID.at(2).Fill(aoqevt[0],zetevt[0]);
            PID.at(3).Fill(aoqevt[4],zetevt[4]);
        }}
        eventcounter++;
    }
    output.mkdir("PID");
    output.cd("PID");
    for(auto &elem: PID) elem.Write();
    output.cd("");
}

void plastics(TTreeReader &datree, TFile &output, vector<bool> &goodevents){
    // This function aims to rebuild the trigger from Q1*Q2 F7plastic
    // therefore we set a limit on sqrt(Q1*Q2) to suppress random noise
    printf("Now beginning with reconstruction of plastic scintillators ...\n");
    const int numplastic = 4;
    const int F7pos = 1;
    const int threshhold = 450; // random noise suppresion value 
    datree.Restart();
    TTreeReaderArray<int> plasticQleft (datree, "BigRIPSPlastic.fQLRaw");
    TTreeReaderArray<int> plasticQright(datree, "BigRIPSPlastic.fQRRaw");    
    TTreeReaderArray<int> plastictleft (datree, "BigRIPSPlastic.fTLRaw");
    TTreeReaderArray<int> plastictright(datree, "BigRIPSPlastic.fTRRaw");    
    TTreeReaderArray<int> plasticnum   (datree, "BigRIPSPlastic.fpl");

    // generate output diagrams
    vector<TH2D> qcorr2D;  // Charge Correlation between 1 & 2
    vector<TH1D> qcorr;    // deposited charge at detector
    vector<TH2D> tqcorr2D; // Charge-Time correlation
    
    vector<vector<string>> arrayname, arraytitle;
    datree.Next(); // Get first event for invovled focal planes
    for(uint i=0; i<numplastic; i++){
        arrayname.push_back(vector<string>{
            "f7pltrigQ."  + to_string(plasticnum[i]), 
            "PlasticQ2D." + to_string(plasticnum[i]),
            "tQcorr."     + to_string(plasticnum[i])
        });

        arraytitle.push_back(vector<string>{
            "Charge deposition by Plastic F"     + to_string(plasticnum[i]),
            "Charge distribution by Plastic F"   + to_string(plasticnum[i]),
            "Charge Ratio/time difference Pl. F" + to_string(plasticnum[i])
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

    while(datree.Next()){
        if(sqrt(plasticQleft[F7pos]*plasticQright[F7pos])> threshhold)
        for(uint i=0; i<numplastic; i++){
            if((plasticQleft[i] >0 )&& (plasticQright[i] >0)){
                qcorr2D.at(i).Fill(plasticQleft[i],plasticQright[i]);
                qcorr.at(i).Fill(sqrt(plasticQleft[i]*plasticQright[i]));
                if((plastictleft[i] >0) && (plastictright[i] >0) ){
                    tqcorr2D.at(i).Fill(plastictleft[i]-plastictright[i],
                        TMath::Log(plasticQleft[i]/(double)plasticQright[i]));
                }
            }
        }
        for(int j=0;j<numplastic;j++){
            goodevents[i] =goodevents[i]*
                       (sqrt(plasticQleft[j]*plasticQright[j]) > range[j][0])*
                       (sqrt(plasticQleft[j]*plasticQright[j]) < range[j][1]);
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

void ppacs(TTreeReader &datree, TFile &output, vector<bool> &goodevents){
    printf("Now beginning with reconstruction of the ppac's ...\n");
    const int numplane = 36;
    const int pl11position = 3; //Plastic at F11 is fourth in array
    datree.Restart();

    TTreeReaderArray<double> triggerpl11(datree,  "BigRIPSPlastic.fTime");
    TTreeReaderArray<bool>   bigripsx   (datree,  "BigRIPSPPAC.fFiredX");
    TTreeReaderArray<bool>   bigripsy   (datree,  "BigRIPSPPAC.fFiredY");
    TTreeReaderArray<double> bigxsum    (datree,  "BigRIPSPPAC.fTSumX");
    TTreeReaderArray<double> bigxdiff   (datree,  "BigRIPSPPAC.fTDiffX");
    TTreeReaderArray<double> bigysum    (datree,  "BigRIPSPPAC.fTSumY");
    TTreeReaderArray<double> bigydiff   (datree,  "BigRIPSPPAC.fTDiffY");
    
    TTreeReaderArray<TString> bripsname(datree, "BigRIPSPPAC.name");
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
    //datree.Next();
    //cout << "T STRING: " << bripsname.operator[](0).Data() << endl;
/*    for(int i=1; i<=numplane; i++){
        effPPACX.GetXaxis()->SetBinLabel(i, bripsname[i-1].Data());
        effPPACY.GetXaxis()->SetBinLabel(i, bripsname[i-1].Data());
    }*/
    
    // 0=X, 1=Y  || 0...35 Plane No.
    uint total = 0;
    uint fulltotal =0;  //Skipping first event for name generation

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
    while(datree.Next()){
        // Get Efficiency relative to the last Plastic
        if (triggerpl11[pl11position] >0){
            for(int i=1; i<=numplane; i++){
                //ripsvec.at(0).at(i) 
                if(bigripsx[i-1]) 
                    effPPACX.SetBinContent(i, effPPACX.GetBinContent(i)+1);
                if(bigripsy[i-1]) 
                    effPPACY.SetBinContent(i, effPPACY.GetBinContent(i)+1);
            }
            total++;
        }
        for(int i=0; i<numplane;i++){
            sumdiffppac.at(0).at(0).Fill(i,bigxsum[i]);
            sumdiffppac.at(0).at(1).Fill(i,bigxdiff[i]);
            sumdiffppac.at(1).at(0).Fill(i,bigysum[i]);
            sumdiffppac.at(1).at(1).Fill(i,bigydiff[i]);
        }

        // Make Cut conditions: each focal Plane needs to have one valid entry
        // := Sum for x and y is in range
        for(uint i =0; i<ppacplane.size();i++){
            for(uint j=0; j<4; j++){ // Loop over all 4 PPAC's per Focal Plane
            temptruth.at(j) =
                (bigxsum[ppacplane.at(i).at(j)] > ppacrange.at(i).at(j).at(0) &&
                 bigxsum[ppacplane.at(i).at(j)] < ppacrange.at(i).at(j).at(1) &&
                 bigysum[ppacplane.at(i).at(j)] > ppacrange.at(i).at(j).at(2) &&
                 bigysum[ppacplane.at(i).at(j)] < ppacrange.at(i).at(j).at(3));
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

void ionisationchamber(TTreeReader &datree, TFile &output,
                       vector<bool> &goodevents) {
    // This method aims to control the Ionisation Chamber values
    // therefore only certain events are accepted

    printf("Now beginning analysis of the IC's (Fpl 7 & 11)\n");
    datree.Restart();

    vector<TTreeReaderArray<int>> icvals{
            TTreeReaderArray<int>(datree, "BigRIPSIC.nhitchannel"),
            TTreeReaderArray<int>(datree, "BigRIPSIC.fADC[32]")};

    TTreeReaderArray<int[32]> testval(datree, "BigRIPSIC.fADC[32]");
    const int numchannel = 6; // There are 6 channel per IC
    const int numic = 2; // Number of Ionisation Chambers
    vector<vector<TH2D>> comparediag;

    for(int i=1; i<numchannel; i++){
        vector<string> arrname = {"ICratio" + to_string(i) + "to0fpl7",
                                  "ICratio" + to_string(i) + "to0fpl11"};
        string arrtitle = "Peak ADC Ratio " + to_string(i) + " to 0";
        comparediag.push_back({TH2D(arrname.at(0).c_str(),arrtitle.c_str(),
                                    2048,0,16384,2048,0,16384),
                               TH2D(arrname.at(1).c_str(),(arrtitle+ "BS").c_str(),
                                    2048,0,16384,2048,0,16384)});
        for(auto &elem: comparediag.back()){
            elem.SetOption("colz");
            elem.GetXaxis()->SetTitle("ADC0 [ch]");
            elem.GetYaxis()->SetTitle(("ADC"+to_string(i)+ "[ch]").c_str());
        }
    }
    uint totalcounter =0;

    while(datree.Next()){
        printf("Getting first nhits %i, %i", testval[1],
               icvals.at(1).At(1));
        if(icvals.at(0).At(0)*icvals.at(0).At(1) <16) //Require 4 hits in each IC
            goodevents.at(totalcounter) = false;
        if((icvals.at(1).At(0) <0) || (icvals.at(1).At(32) <0)) continue;
        for(int i=1; i<numchannel;i++){
            if(icvals.at(1).At(numchannel) >0)
                comparediag.at(0).at(numchannel-1).Fill(icvals.at(1).At(0),
                                                   icvals.at(1).At(numchannel));
            if(icvals.at(1).At(numchannel+32) >0)
                comparediag.at(1).at(numchannel-1).Fill(icvals.at(1).At(32),
                                                icvals.at(1).At(32+numchannel));
        }
        totalcounter++;
    }

    output.mkdir("IC/IC7");
    output.mkdir("IC/IC11");
    output.cd("IC/IC7");
    for(auto &elem: comparediag.at(0)) elem.Write();
    output.cd("IC/IC11");
    for(auto &elem: comparediag.at(1)) elem.Write();
    output.cd("");
}
void highordercorrection(TTreeReader &datree, TFile &output){
    // This method aims to determine the higher-order corrections
    // for the matrix elements
    printf("Now beginning with higher order corrections ...\n");
    datree.Restart();
    vector<TTreeReaderValue<double>> culprit{
        TTreeReaderValue<double> (datree,  "F3X"),
        TTreeReaderValue<double> (datree,  "F3A"),
        TTreeReaderValue<double> (datree,  "F5X"),
        TTreeReaderValue<double> (datree,  "F5A")};

    TTreeReaderArray<double> aoqevt (datree, "BigRIPSBeam.aoq");
    
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
    
    while(datree.Next()){
        for(uint i=0; i<arrname.size();i++){
            culpritdiag.at(i).Fill(aoqevt[beam], *culprit.at(i));
        }
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
    const string output = "build/output/histout.root";
    TFile outputfile(output.c_str(), "RECREATE");
    if(!outputfile.IsOpen()) __throw_invalid_argument("Output file not valid");

    vector<bool> options{
        false,  // Plastics
        false,  // ppac
        true,  // Ionisationchamber
        true,  // highordercorrection
        true   // makepid
    };

    // Store events that cannot be used
    vector<bool> goodevents(1000000, true);

    // Then we read in the tree (lazy)
    TTreeReader mytreereader("tree", &inputfile);

    // Rebuild F7 Trigger combined charge threshhold
    if(options.at(0)) plastics(mytreereader, outputfile, goodevents);
    
    // Cut the PPAC's
    if(options.at(1)) ppacs(mytreereader, outputfile, goodevents);

    if(options.at(2)) ionisationchamber(mytreereader, outputfile, goodevents);

    // Get Corrections
    if(options.at(3)) highordercorrection(mytreereader, outputfile);

    // Get Z vs. A/Q
    if(options.at(4)) makepid(mytreereader, outputfile, goodevents);
    printf("Made PID histograms in %s\n", output.c_str());

    cout << "Runs has " <<accumulate(goodevents.begin(),goodevents.end(),0)
         << " good Elements" << endl;

    outputfile.Close();
}