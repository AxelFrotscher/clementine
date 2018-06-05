#include <MakeAllTree_78Ni.hh>
#include "histograms.hh"
#include "histogram_cuts.hh"

using namespace std;

calibpar p1;

void highordercorrection(treereader &tree, TFile &output){
    // This method aims to determine the higher-order corrections
    // for the matrix elements
    printf("Now beginning with higher order corrections ...\n");
    vector<string> keys{"F5X","F5A", "F3X", "F3A", "F9A", "F9X",
                        "F11A", "F11X", "BigRIPSBeam.aoq",
                        "BigRIPSBeam.beta", "BigRIPSBeam.zet",};
    tree.setloopkeys(keys);
    
    const vector<int> beam{0,4}; // Evaluate Beam F3-7 (1st element)
    const int corrcount = 3;     // Number of corrections
    const vector<vector<string>> arrname = {
            { "DepF3X","DepF3A","DepF5X", "DepF5A", "DepBetaF3-7"},
            { "DepF9X", "DepF9A","DepF11X", "DepF11A", "DepBetaF8-11"}};

    const vector<vector<string>> arrtitle = {
            {"Dependence of F3X vs AoQ", "Dependence of F3A vs AoQ",
             "Dependence of F5X vs AoQ", "Dependence of F5A vs AoQ",
             "Dependence of #beta (F7) vs AoQ"},
            {"Dependence of F9X vs AoQ", "Dependence of F9A vs AoQ",
             "Dependence of F11X vs AoQ", "Dependece of F11A vs AoQ",
             "Dependence of #beta (F11) vs AoQ"}};
    
    vector<vector<vector<TH2D>>> culpritdiag;
    vector<vector<double>> cutval{ // for corrections we use 85Ge
        {2.6449,      // center x
         42.0,         // center y
         0.008,        // radius x
         0.6     },    // radius y
        {2.6437,          // center x
         41.86,        // center y
         0.008,        // radius x
         0.5     }     // radius y
    };

    // Initialize all the diagrams
    for(uint i=0; i<arrname.size();i++){ // Loop F7, F11
        vector<vector<TH2D>> temp2d;
        for(uint j=0; j<arrname.at(0).size();j++){ // Loop F3X,F5X,...
            vector<TH2D> temp1d;
            for(uint k=0;k<=corrcount; k++){ // Loop Correction number
                string arr = arrname.at(i).at(j) + to_string(k); // Generate array name
                string arrn = arrtitle.at(i).at(j) + " Corr: " + to_string(k);
                double ymax = 150;
                double ymin = -100;
                if(j == 4) { // Beta diagrams need smaller bins
                    ymax = 0.7;
                    ymin = 0.6 -0.2*i; // outgoing have lower velocities
                }
                temp1d.push_back(TH2D(arr.c_str(),arrn.c_str(),2000,2.2,3.2,
                                      500,ymin,ymax));
                temp1d.back().SetOption("colz");
                temp1d.back().GetXaxis()->SetTitle("A/Q");
                temp1d.back().GetYaxis()->SetTitle(arrname.at(i).at(j).c_str());
            }
            temp2d.push_back(temp1d);
        }
        culpritdiag.push_back(temp2d);
    }

    printf("Successfully generated Histograms for higher order...\n");
    //Attempting first real correction with F5 x-position

    p1.F7absF5X = 2.644;
    p1.F7linF5X = -1.448E-5;
    p1.F7linF5A = -0.0002524;
    p1.F7linF3X = 0.000103;
    p1.F7absF5X0 = cutval[0][0];

    p1.F11absF9X0 = cutval[1][0];
    p1.F11absF9X = 2.647;
    p1.F11linF9X = -5.688E-5;
    p1.F11linF9A = 0.0002825;
    p1.F11linF11A = 7.872E-5;

    vector<vector <double>> fillvals(2,vector<double>(9,0)); // Fill dependent variable and

    // Progress Bar setup
    int eventno=0; // counting variable
    Long64_t totevents = tree.NumEntries();
    const int downscale = 500; // every n-th event

    while(tree.singleloop()){
        fillvals.at(0).at(0) = tree.F3X;
        fillvals.at(0).at(1) = tree.F3A;
        fillvals.at(0).at(2) = tree.F5X;
        fillvals.at(0).at(3) = tree.F5A;
        fillvals.at(1).at(0) = tree.F9X;
        fillvals.at(1).at(1) = tree.F9A;
        fillvals.at(1).at(2) = tree.F11X;
        fillvals.at(1).at(3) = tree.F11A;

        // Fill corrected Values
        fillvals.at(0).at(6) = tree.BigRIPSBeam_aoq[beam.at(0)] + p1.F7absF5X0 -
                               (p1.F7absF5X +tree.F5X*p1.F7linF5X);
        fillvals.at(0).at(7) = fillvals.at(0).at(6) + tree.F5A*p1.F7linF5A;
        fillvals.at(0).at(8) = fillvals.at(0).at(7) + tree.F3X*p1.F7linF3X;

        fillvals.at(1).at(6) = tree.BigRIPSBeam_aoq[beam.at(1)] + p1.F11absF9X0-
                               (p1.F11absF9X+tree.F9X*p1.F11linF9X);
        fillvals.at(1).at(7) = fillvals.at(1).at(6) - tree.F9A*p1.F11linF9A;
        fillvals.at(1).at(8) = fillvals.at(1).at(7) - tree.F11A*p1.F11linF11A;

        // Fill pre and post events
        for(uint ii=0; ii<beam.size();ii++){
            if((pow(1./cutval.at(ii).at(2)*(tree.BigRIPSBeam_aoq[beam.at(ii)] -
                                           cutval.at(ii).at(0)),2) +
                pow(1/cutval.at(ii).at(3)*(tree.BigRIPSBeam_zet[beam.at(ii)] -
                                          cutval.at(ii).at(1)),2)) <1){
                // Applying the elliptic cut for 85Ge
                fillvals.at(ii).at(4) = tree.BigRIPSBeam_beta[beam.at(ii)];
                fillvals.at(ii).at(5) = tree.BigRIPSBeam_aoq[beam.at(ii)];

                for(uint k=0; k<culpritdiag.at(0).size(); k++){
                    for(uint j=0; j<=corrcount; j++){
                        culpritdiag.at(ii).at(k).at(j).Fill(fillvals.at(ii).at(5+j),
                                                           fillvals.at(ii).at(k));
                    }
                }
            }
        }
        eventno++;
        if(!(eventno%downscale)) progressbar(eventno,totevents);
    }
    printf("\nSuccessfully looped for higher order...\n");
    // Get linear fits of the projected means
    vector<vector<vector<TProfile *>>> projections;
    const double cutfrac = 0.8; // fraction to consider for fit
    auto corrlinfit = new TF1("Linear Fit", linfit, cutval[0][1]-cutfrac*cutval[0][3],
                               cutval[0][1]+cutfrac*cutval[0][3], 2);
    corrlinfit->SetParNames("absolute", "linear");
    for (auto &elem : culpritdiag) {
        vector<vector<TProfile*>> proftemp2d;
        for(auto &elem2: elem){
            vector<TProfile*> proftemp1d;
            for(uint i_corr=0; i_corr<=corrcount;i_corr++){
                proftemp1d.push_back(elem2.at(i_corr).ProfileY());
                proftemp1d.back()->Fit("Linear Fit");
            }
            proftemp2d.push_back(proftemp1d);
        }
        projections.push_back(proftemp2d);
    }

    // Setting folder structure, this should correspond to the order of corrections
    const vector<vector<string>> folders{
        {"Corrections/Pre/Raw",
         "Corrections/Pre/F5Xcorr",
         "Corrections/Pre/F5XF5Acorr",
         "Corrections/Pre/F5XF5AF3X"},
        {"Corrections/Post/Raw",
         "Corrections/Post/F9Xcorr",
         "Corrections/Post/F9Acorr",
         "Corrections/Post/F11Acorr"}};

    for (auto &i_fold: folders) for(auto &j: i_fold) output.mkdir(j.c_str());

    for(uint i=0; i<projections.size();i++){  // Pre, Post
        for(uint k=0; k<projections.at(0).at(0).size(); k++){ // Corr 0,1,2,

            string tempfolder = folders.at(i).at(k) + "/Profiles";
            output.mkdir(tempfolder.c_str());

            for(uint j=0; j<projections.at(0).size();j++){ // F3X, F5X, ...
                output.cd(tempfolder.c_str());
                projections.at(i).at(j).at(k)->Write();

                output.cd(folders.at(i).at(k).c_str());
                culpritdiag.at(i).at(j).at(k).Write();
            }
        }
    }

    output.cd("");
    printf("Finished with higher order corrections!\n");
}

void dalicalib(treereader &tree, TFile &output){
    // This Method aims to calibrate the 187 detectors of DALI
    printf("Now beginning the Calibration of the NaI crystals... \n");

    vector<string> keys{"DALINaI", "DALINaI.fADC", "DALINaI.id"};
    tree.setloopkeys(keys);

    TH2D gammadetectors("dalispectra", "Spectrum of each Gamma Detector",
                        186,0,186,4096,0,4096);
    gammadetectors.SetOption("colz");
    gammadetectors.GetXaxis()->SetTitle("Detector Number");
    gammadetectors.GetYaxis()->SetTitle("ADC Channel");

    int numdet =0;

    // Progress Bar setup
    int currevt=0; // counting variable
    Long64_t totevents = tree.NumEntries();
    const int downscale = 500; // every n-th event

    while(tree.singleloop()){
        numdet = tree.DALINaI_;
        for(int i=0; i<numdet;i++){
            if(tree.DALINaI_fADC[i]) gammadetectors.Fill(tree.DALINaI_id[i],
                tree.DALINaI_fADC[i]);
        }
        currevt++;
        if(!(currevt%downscale)) progressbar(currevt,totevents);
    }

    output.mkdir("DALI");
    output.cd("DALI");
    gammadetectors.Write();
    output.cd("");
    printf("\nFinished DALI Calibration");
}

void makepid(treereader &tree, TFile &output, const vector<bool> &goodevents){
    // 5 beams, 2incoming, 3outgoing
    vector <string> keys{"BigRIPSBeam.aoq", "BigRIPSBeam.zet", "F3X", "F3A",
                         "F5X", "F5A", "F9X", "F9A", "F11X", "F11A"};
    tree.setloopkeys(keys);

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

    for(auto &elem: PID){
        for(auto &uelem: elem){
            uelem.SetOption("colz");
            uelem.GetXaxis()->SetTitle("A/Q");
            uelem.GetYaxis()->SetTitle("Z");
        }
    }

    vector<double> cutval{
            2.7084, // center x
            41.0,  // center y
            0.009, // radius x
            0.5    // radius y
    };

    // Progress Bar setup
    uint eventcounter =0;
    Long64_t totevents = tree.NumEntries();
    const int downscale = 500; // every n-th event

    vector<vector<double>> valinc;     // Store temporary beam values
    while(tree.singleloop()){
        if(goodevents.at(eventcounter)){
            double beamaoqcorr = tree.BigRIPSBeam_aoq[0] + p1.F7absF5X0 -
                                 (p1.F7absF5X+tree.F5X*p1.F7linF5X) +
                                 tree.F5A*p1.F7linF5A +
                                 tree.F3X*p1.F7linF3X;
            double beamaoqcorr2 = tree.BigRIPSBeam_aoq[4] +  p1.F11absF9X0-
                                  (p1.F11absF9X+tree.F9X*p1.F11linF9X) -
                                  tree.F9A*p1.F11linF9A - tree.F11A*p1.F11linF11A;

            //Loop over all elements in the tree
            valinc.push_back({tree.BigRIPSBeam_aoq[0],tree.BigRIPSBeam_aoq[1]});
            valinc.push_back({tree.BigRIPSBeam_zet[0],tree.BigRIPSBeam_zet[1]});
            if(closeness(valinc.at(0)) && closeness(valinc.at(1))){
                // Cut Particles that have variating aoq or zet
                PID.at(0).at(0).Fill(tree.BigRIPSBeam_aoq[0],
                                     tree.BigRIPSBeam_zet[0]);
                PID.at(1).at(0).Fill(beamaoqcorr, tree.BigRIPSBeam_zet[0]);
            }

            valinc.push_back({tree.BigRIPSBeam_aoq[2],tree.BigRIPSBeam_aoq[3],
                              tree.BigRIPSBeam_aoq[4]});
            valinc.push_back({tree.BigRIPSBeam_zet[2],tree.BigRIPSBeam_zet[3],
                              tree.BigRIPSBeam_zet[4]});
            if(closeness(valinc.at(2)) && closeness(valinc.at(3))){
                PID.at(0).at(1).Fill(tree.BigRIPSBeam_aoq[4],
                                     tree.BigRIPSBeam_zet[4]);
                PID.at(1).at(1).Fill(beamaoqcorr2, tree.BigRIPSBeam_zet[4]);
            }
            valinc.clear();

            // We now fill the cutted data (cut by ellipsoid)
            if((pow(1./cutval.at(2)*(tree.BigRIPSBeam_aoq[0]-cutval.at(0)),2) +
                pow(1/cutval.at(3)*(tree.BigRIPSBeam_zet[0]-cutval.at(1)),2))<1){
                PID.at(0).at(2).Fill(tree.BigRIPSBeam_aoq[0],
                                     tree.BigRIPSBeam_zet[0]);
                PID.at(1).at(2).Fill(beamaoqcorr, tree.BigRIPSBeam_zet[0]);
                PID.at(0).at(3).Fill(tree.BigRIPSBeam_aoq[4],
                                     tree.BigRIPSBeam_zet[4]);
                PID.at(1).at(3).Fill(beamaoqcorr2, tree.BigRIPSBeam_zet[4]);
            }
        }
        eventcounter++;
        if(!(eventcounter%downscale)) progressbar(eventcounter,totevents);
    }
    printf("\nFinished making PIDs!\n");
    vector<string> folders{"PID/Uncorrected", "PID/Corrected"};
    for (auto i :folders) output.mkdir(i.c_str());

    for(uint i=0; i<folders.size(); i++){
        output.cd(folders.at(i).c_str());
        for(auto &elem: PID.at(i)) elem.Write();
    }
    output.cd("");
}

void makehistograms(const string input){
    TTree* inputtree;
    TFile f(input.c_str());
    f.GetObject("tree", inputtree);

    treereader alt2dtree(inputtree); // Opens the input file...

    // Generating an outputfile that matches names with the input file
    string output = input.substr(0,22) + "hist.root";

    TFile outputfile(output.c_str(), "RECREATE");
    if(!outputfile.IsOpen()) __throw_invalid_argument("Output file not valid");

    vector<bool> options{
        true,  // Plastics
        true,  // ppac
        true,  // Ionisationchamber
        true,  // highordercorrection
        true,  // makepid
        true   // DALIcalib
    };

    // Store events that cannot be used
    vector<bool> goodevents(alt2dtree.NumEntries(), true);

    // Then we read in the tree (lazy)
    //TTreeReader mytreereader("tree", &inputfile);

    // Rebuild F7 Trigger combined charge threshhold
    if (options.at(0)) plastics(alt2dtree, outputfile, goodevents);

    // Cut the PPAC's
    if (options.at(1)) ppacs(alt2dtree, outputfile, goodevents);

    if (options.at(2)) ionisationchamber(alt2dtree, outputfile, goodevents);
    printf("Finished with PPAC Consistency checks\n");
    // Get Corrections
    if (options.at(3)) highordercorrection(alt2dtree, outputfile);

    // Get Z vs. A/Q
    if (options.at(4)) makepid(alt2dtree, outputfile, goodevents);
    printf("Made PID histograms in %s\n", output.c_str());
    //Get ADC Spectra for DALI
    if(options.at(5)) dalicalib(alt2dtree, outputfile);

    cout << "Runs has " <<accumulate(goodevents.begin(),goodevents.end(),0)
         << " good Elements" << endl;

    //outputfile.Close();
}
