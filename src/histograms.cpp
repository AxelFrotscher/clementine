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

void addmean(vector<vector<double>> &v){
    // This function takes a 2dim Vector, calculates the mean for each row
    // and pushes back the mean values in a new row
    vector<double> tempmean;
    for (auto &el: v)
        tempmean.push_back(accumulate(el.begin(), el.end(), 0.0)/el.size());

    if(!v.empty()) v.push_back(tempmean);
}

void makepid(TTreeReader &datree, TFile &output){
    datree.Restart();
    // 5 beams, 2incoming, 3outgoing
    TTreeReaderArray<double> aoqevt (datree, "BigRIPSBeam.aoq");
    TTreeReaderArray<double> zetevt (datree, "BigRIPSBeam.zet");

    TH2D PIDincoming("pidinc", "PID Incoming F3-F7", 300,2.45,2.9, 200,25,45);
    TH2D PIDoutgoing("pidout", "PID Outgoing F8-F11", 300,2.45,2.9, 200,25,45);

    PIDincoming.SetOption("colz");
    PIDoutgoing.SetOption("colz");
    PIDincoming.GetXaxis()->SetTitle("A/Q");
    PIDincoming.GetYaxis()->SetTitle("Z");
    PIDoutgoing.GetXaxis()->SetTitle("A/Q");
    PIDoutgoing.GetYaxis()->SetTitle("Z");

    vector<vector<double>> valinc;     // Store temporary beam values
    while(datree.Next()){
        //Loop over all elements in the tree
        valinc.push_back({aoqevt[0],aoqevt[1],aoqevt[2]});
        valinc.push_back({zetevt[0],zetevt[1],zetevt[2]});
        if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
            // Cut Particles that have variating aoq or zet
            //addmean(valinc);
            PIDincoming.Fill(valinc.at(0).at(0),valinc.at(1).at(0));

        }
        valinc.clear();
        valinc.push_back({aoqevt[3],aoqevt[4]});
        valinc.push_back({zetevt[3], zetevt[4]});
        if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
            //addmean(valinc);
            PIDoutgoing.Fill(valinc.at(0).at(1),valinc.at(1).at(1));
        }
        valinc.clear();
    }
    output.mkdir("PID");
    output.cd("PID");
    PIDincoming.Write();
    PIDoutgoing.Write();
    output.cd("");
}

void plastics(TTreeReader &datree, TFile &output){
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
    }
    
    // To avoid multihit-triggers we define a range of acceptance
    vector<vector<int>> range{{400,500},{450,750},{213,350}};
    
    printf("Writing Plastic 'Histograms to disk ... \n");
    
    output.mkdir("Plastics/2D");
    output.mkdir("Plastics/Q1Q2");
    output.mkdir("Plastics/TQCorr");
    output.cd("Plastics/2D");
    for(auto i: qcorr2D) i.Write();
    output.cd("Plastics/Q1Q2");
    for(auto i: qcorr) i.Write();
    output.cd("Plastics/TQCorr");
    for(auto i: tqcorr2D) i.Write();
    output.cd("");

    // A careful analysis yields a threshhold of sqrt(Q1*Q2)>160
    printf("Finished Writing plastic histograms! \n");
}

void ppacs(TTreeReader &datree, TFile &output){
    printf("Now beginning with reconstruction of the ppac's ...\n");
    const int numplane = 36;
    datree.Restart();

    TTreeReaderArray<double> trigger7(datree,  "F7X");
    TTreeReaderArray<bool>   bigripsx(datree,  "BigRIPSPPAC.fFiredX");
    TTreeReaderArray<bool>   bigripsy(datree,  "BigRIPSPPAC.fFiredY");
    TTreeReaderArray<double> bigxsum  (datree,  "BigRIPSPPAC.fTSumX");
    TTreeReaderArray<double> bigxdiff (datree,  "BigRIPSPPAC.fTDiffX");
    TTreeReaderArray<double> bigysum  (datree,  "BigRIPSPPAC.fTSumY");
    TTreeReaderArray<double> bigydiff (datree,  "BigRIPSPPAC.fTDiffY");
    
    //TTreeReaderArray<TString> bripsname(datree, "BigRIPSPPAC.name");
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
    effPPACX.GetYaxis()->SetTitle("PPACX(x)/PPACX(34)");
    effPPACY.GetXaxis()->SetTitle("PPAC [ch]");
    effPPACY.GetYaxis()->SetTitle("PPACY(x)/PPACX(34)");


    //Fill names with first event
    //datree.Next();
    //for(int i=1; i<=numplane; i++){
    //    effPPACX.GetXaxis()->SetBinLabel(i, bripsname[i-1].);
    //    effPPACY.GetXaxis()->SetBinLabel(i, bripsname[i-1].c_str());
    //}
    
    // 0=X, 1=Y  || 0...35 Plane No.
    int total = 0;
    
    while(datree.Next()){
        // Get Efficiency relative to the last PPAC
        if (bigripsx[33]){
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
    
    const int beam = 0; // Evaluate Beam F8-11 (5th element)
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
    const string output = "output/histout.root";
    TFile outputfile(output.c_str(), "RECREATE");
    if(!outputfile.IsOpen()) __throw_invalid_argument("Output file not valid");

    vector<bool> options{
        true,  // Plastics
        true,  // ppacs
        true,  // makepid
        true    // highordercorrection
    };
    
    // Then we read in the tree (lazy)
    TTreeReader mytreereader("tree", &inputfile);

    // Rebuild F7 Trigger combined charge threshhold
    if(options.at(0)) plastics(mytreereader, outputfile);
    
    // Understand Triggers
    if(options.at(1)) ppacs(mytreereader, outputfile);

    // Get Z vs. A/Q
    if(options.at(2)) makepid(mytreereader, outputfile);
    printf("Made PID histograms in %s\n", output.c_str());
    
    // Get Corrections
    if(options.at(3)) highordercorrection(mytreereader, outputfile);

    outputfile.Close();
}