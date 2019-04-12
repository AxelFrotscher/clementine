//#include <Rtypes.h>
#include "MakeAllTree_78Ni.hh"
#include "progress.h"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtRecoPID.hh"
#include "TVectorD.h"
#include "TFile.h"
#include "TArtPPAC.hh"
#include "TArtCalibPPAC.hh"
#include "TArtEventInfo.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtFocalPlane.hh"
#include <iostream>
#include <TH2F.h>
#include <TF1.h>
#include "TGraph.h"
#include "TArtMINOSParameters.hh"
#include "TArtCalibMINOS.hh"
#include "TMath.h"
#include "TMinuit.h"
#include <Math/Vector3D.h>
#include "zdssetting.h"
#include "minos.h"

R__LOAD_LIBRARY(libanacore.so)
R__LOAD_LIBRARY(libminos.so)

using std::cout, std::endl, std::string, std::vector,
      std::__throw_invalid_argument, std::ifstream, std::stringstream, std::max,
      std::__throw_bad_function_call;

void generatetree(const string infile, const string output) {
    //  signal(SIGINT,stop_interrupt); // CTRL + C , interrupt
    cout << "Now in Estore -> " << infile << endl;

    // Create StoreManager both for calibration "TArtCalib..." and 
    // treatment "TArtReco..."
    auto *sman = TArtStoreManager::Instance();
    // Create EventStore to control the loop and get the EventInfo
    TArtEventStore estore;
    //estore.SetInterrupt(&stoploop);
    if (!estore.Open(infile.c_str()))
        __throw_invalid_argument(("Could not open" +
                                  infile).c_str());

    uint setting = getset(estore); // 0 -> 2014, 1 -> 2015
    // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes 
    // parameters from ".xml" files
    TArtBigRIPSParameters para;
    vector<vector<string>> parameterfiles{
        {"../config/db2014/BigRIPSPPAC.xml", "../config/db2014/BigRIPSPlastic.xml",
         "../config/db2014/BigRIPSIC.xml",   "../config/db2014/FocalPlane.xml"},
        {"../config/db/BigRIPSPPAC.xml",     "../config/db/BigRIPSPlastic.xml",
         "../config/db/BigRIPSIC.xml",       "../config/db/FocalPlane.xml"}};

    for (auto i : parameterfiles.at(setting)) {
        if (!para.LoadParameter(&i[0])) __throw_invalid_argument(&i[0]);
    }

    para.SetFocusPosOffset(8, 138.5);
    // Create CalibPID to get and calibrate raw data ( CalibPID -> [CalibPPAC , 
    // CalibIC, CalibPlastic , CalibFocalPlane] 
    auto brcalib = new TArtCalibPID();
    auto ppaccalib = brcalib->GetCalibPPAC();
    auto cfpl = brcalib->GetCalibFocalPlane();

    // Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z, ... 
    // (RecoPID -> [ RecoTOF , RecoRIPS , RecoBeam] )
    TArtRecoPID recopid;

    // Definition of observables we want to reconstruct

    printf("Defining BigRIPS parameters\n");
    const vector<int> focalplanes{3, 5, 7, 8, 9, 11};

    vector<vector<string>> mfil{  // matrixfiles
        {"../config/matrix/mat1.mat",                   "D3"},  // F3 - F5 => D3
        {"../config/matrix/mat2.mat",                   "D5"},  // F5 - F7 => D5
        {"../config/matrix/F8F9_LargeAccAchr.mat",      "D7"},  // F8 - F9 => D7
        {"../config/matrix2014/F9F11_LargeAccAchr.mat", "D8"},  // F9-  F11=> D8
        {"../config/matrix/F9F11_LargeAccAchr.mat",     "D8"}}; // F9 - F11=> D8

    vector<TArtRIPS *> rips{
            recopid.DefineNewRIPS(3, 5, &mfil[0][0][0], &mfil[0][1][0]),
            recopid.DefineNewRIPS(5, 7, &mfil[1][0][0], &mfil[1][1][0]),
            recopid.DefineNewRIPS(8, 9, &mfil[2][0][0], &mfil[2][1][0]),
            recopid.DefineNewRIPS(9, 11, &mfil[3 + setting][0][0],
                                  &mfil[3 + setting][1][0])};

    // Reconstruction of TOF DefineNewTOF(first plane,second plane, time offset)
    vector<vector<string>> fplname{{"F3pl", "F7pl"},
                                   {"F8pl", "F11pl-1"}};
    vector<vector<double>> tofoff{ //300.25 F3-F7 init -159.45 F8-F11 init
        {300.83, -161},
        {304.17 + 0.17,    // good Offset Value for F3-F7,  empty-target run  300.85
         -161.64 - 0.08}}; // good Offset Value for F8-F11, empty-target run -160.45

    class setting set;
    if(set.getsetname() == "70Fe") tofoff.at(0) = vector<double>{300.91,-158.73};
    if(set.getsetname() == "78Ni") tofoff.at(0) = vector<double>{300.23,-159.19};


    vector<TArtTOF *> tof{
        recopid.DefineNewTOF(&fplname[0][0][0], &fplname[0][1][0],
                             tofoff[setting][0], 5),
        recopid.DefineNewTOF(&fplname[1][0][0], &fplname[1][1][0],
                             tofoff[setting][1], 9)};

    // Reconstruction of IC observables for ID
    vector<TArtBeam *> beam{  // br = BigRIPS, zd = ZeroDegree
        recopid.DefineNewBeam(rips[0], rips[1], tof[0], (char *) "F7IC"),  //br_37
        recopid.DefineNewBeam(rips[1], tof[0], (char *) "F7IC"),           //br_57
        recopid.DefineNewBeam(rips[2], tof[1], (char *) "F11IC"),          //zd_89
        recopid.DefineNewBeam(rips[3], tof[1], (char *) "F11IC"),          //zd_911
        recopid.DefineNewBeam(rips[2], rips[3], tof[1], (char *) "F11IC")};//zd_811

    // Create Minos, TArtMINOSPar Class too incompetent to handle regular construction
    auto setup = new TArtMINOSParameters("MINOSParameters", "MINOSParameters");
    vector<string> minosparameters{"../config/db2014/MINOS.xml",
                                   "../config/db/MINOS.xml"};
    setup->LoadParameters(&minosparameters.at(setting)[0]);

    auto minoscalib = new TArtCalibMINOS; // regular construction yields double free

    // Create TTree and output file for it
    TFile fout(output.c_str(), "RECREATE");
    if(!fout.IsOpen()) __throw_invalid_argument("Could not open output file!\n");

    auto tree = new TTree("tree","tree");

    // Define data nodes which are supposed to be dumped to tree 
    // EventInfo is important for the fBit information to know the trigger!
    vector<string> datanodes{"EventInfo", "BigRIPSPPAC", "BigRIPSPlastic",
                             "BigRIPSIC","BigRIPSFocalPlane","BigRIPSRIPS",
                             "BigRIPSTOF", "BigRIPSBeam", "EventInfo_FBIT"};

    // Trigger bits are broken in anaroot manual recover (kudos to n.hupin)
    uint EventInfo_FBIT =0;
    TClonesArray * fbitarray;

    for(auto i : datanodes){
        TClonesArray * array;
        if(i=="EventInfo_FBIT") {
            fbitarray = (TClonesArray *)sman->FindDataContainer("EventInfo");
            //tree->Branch("EventInfo_FBIT", &EventInfo_FBIT, "EventInfo_FBIT/I");
        }
        else{
            array = (TClonesArray *)sman->FindDataContainer(&i[0]);
            tree->Branch(array->GetName(), &array);
        }
        printf("%s", array->GetName());
    }

    if(!fbitarray) __throw_bad_function_call();
    //Making new branches
    vector<double> Xpad, Ypad, Qpad;
    tree->Branch("MinosClustX", &Xpad);
    tree->Branch("MinosClustY", &Ypad);
    tree->Branch("MinosClustQ", &Qpad);

    double VDrift = 0.03786,       // in cm/(1E-6*s)
           DelayTrigger = 2055,    // in ns both Values from
           Stoptime = 9980;        //       ConfigMINOSDirft_fixed.txt
    tree->Branch("VDrift", &VDrift, "MINOS.VDrift/D");
    tree->Branch("DelayTrig", &DelayTrigger, "DelayTrigger/D");

    int filled; // Number of tracks measured
    tree->Branch("Trackamount", &filled);

    vector<vector<double>> minoscalibvalues;
    tree->Branch("Minoscalibvalues", &minoscalibvalues);

    // Coordinates of minos pads hit {{x1,y1},{x2,y2},...}
    vector<vector<double>> minostrackxy;
    tree->Branch("minostrackxy", &minostrackxy);

    vector<vector<double>> minostime;
    tree->Branch("minostime", &minostime);

    ///Make MINOS drift time analysis
    string runname = infile.substr(34,9);
    TH1I timedrift(runname.c_str(), "TPC time drift", 1536,0,15360);
    timedrift.GetXaxis()->SetTitle("t [ns]");
    timedrift.GetYaxis()->SetTitle("trigger number");

    //Parameters for MINOS analysis
    double minosthresh, TimeBinElec, Tshaping;
    vector<string> minosdrift{"../config/db2014/ConfigMINOSDrift.txt",
                              "../config/db/ConfigMINOSDrift_fixed.txt"};

    ifstream minosconfigfile(minosdrift.at(setting));
    assert(minosconfigfile.is_open());

    string trash = " ";
    getline(minosconfigfile, trash);
    minosconfigfile >> minosthresh >> TimeBinElec >> Tshaping;
    getline(minosconfigfile, trash);

    while(getline(minosconfigfile, trash)){
        stringstream ss(trash);
        ss >> trash >> DelayTrigger >> Stoptime >> VDrift;
        if(infile.find(trash) != string::npos){
            cout << endl << "Found parameters for MINOS file "<< trash  << "  "
                 << infile << endl;
            break;
        }
    }
    minosconfigfile.close();

    tree->Branch("TimeBinElec", &TimeBinElec);
    tree->Branch("Tshaping", &Tshaping);

    //BigRIPS
    const double magnum = -999; // Magic Number for unset events
    const int planesperppac = 4;
    vector<vector<double>> fF8PPAC(3, vector<double>(planesperppac, magnum));
 
    // 6 focal planes (3,5,7,8,9,11) and 4 value (X,A,Y,B)
    vector<vector<double>> focal{6, vector<double>(4)}; // F3X, F3A

    vector<int> fGoodPPACFocus(12, 0);
    vector<int> fGoodPPACFocusOr(12, 0);

    // F8 PPAC positions relative to F8 0: x 1: y
    const vector<vector<double>> f8z{{-1310.7,-1302.1},{-1273.3,-1281.9},
                                     {-810.7,-802.1},{-773.3,-781.9}};

    vector<vector<string>> focalname{{"F3X","F3A"}, {"F5X","F5A"}, 
                                     {"F7X","F7A"}, {"F8X","F8A"},
                                     {"F9X","F9A"}, {"F11X","F11A"}};

    for(uint i=0; i<focal.size(); i++){
        tree->Branch(focalname[i][0].c_str(), &focal.at(i).at(0));
        tree->Branch(focalname[i][1].c_str(), &focal.at(i).at(1));
    }

    TVectorD *vec;

    tree->Branch("fgoodppacfocus",fGoodPPACFocus.data(),"fGoodPPACFocus[12]/I");
    tree->Branch("fgoodppacfocusor",fGoodPPACFocusOr.data(),
                 "fGoodPPACFocusOr[12]/I");

    const vector<string> ppacname{
        "F8PPAC-1A", "F8PPAC-1B", "F8PPAC-2A", "F8PPAC-2B",
        "F9PPAC-1A", "F9PPAC-1B", "F9PPAC-2A", "F9PPAC-2B", 
        "F11PPAC-1A","F11PPAC-1B","F11PPAC-2A","F11PPAC-2B"};

    // Progress Bar setup
    int neve = 0; // counting variable
    uint totevents = 5000000; //50000
    const int downscale = 1000; // every n-th event

    progressbar progress(totevents, 0);

    cout << endl << "Finished Initialization. Starting transcription... " << endl;
    while(estore.GetNextEvent() && (neve < totevents)){
        //Making the BigRIPS tree calibration
        brcalib->ClearData();
        brcalib->ReconstructData();

        //Reconstructiong the PID
        recopid.ClearData();
        recopid.ReconstructData();

        //Reconstructing the scattering angle:
        //Beam spot on target reconstruction
        EventInfo_FBIT = 0;

        for(auto &elem: fF8PPAC) fill(elem.begin(), elem.end(), magnum);
        fill(fGoodPPACFocus.begin(), fGoodPPACFocus.end(), 0);
        fill(fGoodPPACFocusOr.begin(), fGoodPPACFocusOr.end(), 0);

        // Reading out the PPAC'S
        vector<TArtPPAC*> fppac; // Storage Vector for all PPAC's
        vector<bool> fired;

        for(auto i: ppacname){ 
            fppac.push_back(ppaccalib->FindPPAC(&i[0])); // Get all PPAC planes
            fired.push_back(fppac.back() && fppac.back()->IsFiredX() &&
                            fppac.back()->IsFiredY());
        }

        // Determine well focused events
        const vector<uint> ppacuse{8,9,11}; // Which PPAC's to use (match w/ ppacname!)

        for(uint i=0; i<ppacuse.size(); i++){ // Loop over every PPAC
            uint ppacno = ppacuse.at(i);
            const int j = planesperppac;
            // If every Plane has a signal, set fGoodPPACFocus(Or)
            if(fppac.at(j*i+0) && fppac.at(j*i+1) && 
               fppac.at(j*i+2) && fppac.at(j*i+3)){
                fGoodPPACFocus.at(ppacno) = 
                    fppac.at(j*i+0)->IsFiredX() && fppac.at(j*i+1)->IsFiredX() && 
                    fppac.at(j*i+2)->IsFiredX() && fppac.at(j*i+3)->IsFiredX();
                fGoodPPACFocusOr.at(ppacno) = 
                    (fppac.at(j*i+0)->IsFiredX() || fppac.at(j*i+1)->IsFiredX()) &&
                    (fppac.at(j*i+2)->IsFiredX() || fppac.at(j*i+3)->IsFiredX());
            }
        }

        // 880mm F8-1B to target

        for(uint i=0; i<focalplanes.size(); i++){ // Loop all focal planes
            if(cfpl->FindFocalPlane(focalplanes.at(i))){
                vec = cfpl->FindFocalPlane(focalplanes.at(i))->GetOptVector();
                for(uint j=0; j<planesperppac; j++){
                    focal.at(i).at(j) = (*vec)(j);
                }
            }
        }

        // Manual filling of the trigger
        auto rawevent = (TArtRawEventObject *)sman->FindDataContainer("RawEvent");
        for(uint i=0; i<rawevent->GetNumSeg(); i++){
            auto seg = rawevent->GetSegment(i);
            if(seg->GetFP()==63 && seg->GetDetector()==10){
                for(uint j=0; j<seg->GetNumData();j++){
                    EventInfo_FBIT = seg->GetData(j)->GetVal();
                }
            }
        }

        ((TArtEventInfo *)fbitarray->At(0))->SetTriggerBit(EventInfo_FBIT);
        //cout << EventInfo_FBIT << endl;

        ////////////////////////////////// Making MINOS ////////////////////////
        minoscalib->ClearData();
        minoscalib->ReconstructData();

        /// MINOS 1. Filling vectors with (x,y) data ///////////////////////////
        filled = 0; // Number of tracks
        Xpad.clear();
        Ypad.clear();
        Qpad.clear();
        minoscalibvalues.clear();
        minostrackxy.clear();
        minostime.clear();
        auto *minos = new TArtCalibMINOSData;
        for(int i=0; i<minoscalib->GetNumCalibMINOS(); i++){
            minos = minoscalib->GetCalibMINOS(i);
            minoscalibvalues.push_back(*minos->GetCalibValueArray());
            minostime.push_back(*minos->GetCalibTimeArray());

            double maxcharge =0;
            double x_mm = minos->GetX();
            double y_mm = minos->GetY();
            minostrackxy.push_back(vector<double>{x_mm, y_mm});

            if(!(abs(x_mm)<0.01 && abs(y_mm)< 0.01)){ // Cut not connected channels
                for(int j=0; j<minos->GetNData(); j++)
                    maxcharge = max(maxcharge, minos->GetCalibValue(j));
                if(maxcharge>minosthresh){
                    Xpad.push_back(x_mm);
                    Ypad.push_back(y_mm);
                    Qpad.push_back(maxcharge);
                    filled++;
                }
            }
        }
        if(!minostrackxy.size()) minostrackxy.push_back({0,0});

        tree->Fill();
        neve++;

        progress.increaseevent();
        if(!(neve%downscale)) progress.draw();

        // Add Graphical Feedback

        vector<bool> clusterringbool;
        vector<int> clusternbr, clusterpads;
        vector<double> Xpadnew, Ypadnew, Qpadnew, Zpadnew;
        int trackNbr = 0;

        /// Do Hough-transform for data
        if (filled) {
            for (int iteration = 0; (Xpad.size() > 9 && iteration < 20); iteration++) {
                int filter_result = minosana::Obertelli_filter(Xpad, Ypad,
                        Qpad, Xpadnew, Ypadnew, Qpadnew, clusterringbool);
                if (filter_result > 10 && clusterringbool.back()) trackNbr++;
            }
        }

        if(trackNbr <1 || trackNbr > 4) continue;


        /// Fill minos tpcdrift histogram
        for(int i=0; i<minoscalibvalues.size(); i++){ /// Loop over all pads

            double x_mm = minostrackxy.at(i).at(0);
            double y_mm = minostrackxy.at(i).at(1);
            bool fitbool = false;

            for (int j = 0; j < Xpadnew.size(); j++) {
                if (abs(Xpadnew[j] - x_mm) < 0.01 && abs(Ypadnew[j] - y_mm) < 0.01) {
                    fitbool = true;
                    break;
                }
            }
            // Check if new channel is of interest
            if (!fitbool) continue;

            /// Generate Q(t) diagram
            //cout << "Pad NO: " << i <<  " Fitbool: " << fitbool << endl;
            TH1F hfit("hfit","hfit", 512, 0, 512);
            for(int j=0; j<minoscalibvalues.at(i).size(); j++){
                if(minoscalibvalues.at(i).at(j) >0)
                    hfit.SetBinContent(hfit.FindBin(minostime.at(i).at(j)),
                            minoscalibvalues.at(i).at(j)+250);
            }

            if(hfit.GetSumOfWeights() == 0){  continue;} // no empty diagrams

            hfit.GetXaxis()->SetRange(0, 510);
            double hfit_max = hfit.GetMaximum();
            double hfit_max_T = hfit.GetMaximumBin();

            // Find T_min and T_max limits for non-0 signals
            double T_min = hfit.FindLastBinAbove(250);
            double T_max = hfit.FindFirstBinAbove(0) - 1;

            // Take only 1.5*shaping time before max if other signals before
            if (hfit_max_T - 3.5 * (Tshaping / TimeBinElec) > T_min)
                T_min = hfit_max_T - 2 * Tshaping / TimeBinElec;
            if ((hfit_max_T + 10) < T_max || T_max == -1) T_max = hfit_max_T + 10;

            T_min = std::max(T_min, 0.);
            T_max = std::min(T_max, 510.);

            // Set fit parameters
            TF1 fit_function("fit_function",
                [](double *x, double *p){ /// Check for boundaries of x
                    if(x[0] < p[1] || x[0] > 512) return 250.;
                    else return p[0] * exp(-3.*(x[0]-p[1])/p[2])*sin((x[0]-p[1])/p[2]) *
                                pow((x[0]-p[1])/p[2], 3) + 250;}, 0, 511, 3);

            fit_function.SetParameters( hfit_max - 250,
                                        hfit_max_T - Tshaping / TimeBinElec,
                                        Tshaping / TimeBinElec);
            fit_function.SetParNames("Amplitude", "trigger time", "shaping time");
            fit_function.SetParLimits(0, 0, 1E5);
            fit_function.SetParLimits(1, -20, 512);
            fit_function.SetParLimits(2, 0, 512);

            // --> parameter n (no store no draw) crucial for multithread <--
            int fit2DStatus = hfit.Fit(&fit_function, "RQN", "", T_min, T_max);

            if(!fit2DStatus && fit_function.GetParameter(1) > .15
                            && fit_function.GetParameter(1) < 512)
                timedrift.Fill(fit_function.GetParameter(1)*TimeBinElec);
        }

    }
    progress.reset();
    cout << "Writing the tree..."<<endl;

    fout.Write();
    fout.Close();
    printf("Writing TTree complete!\n");

    TFile tpcdrift(Form("output/MINOS/%s.root", runname.c_str()), "Update");
    if (!tpcdrift.GetListOfKeys()->Contains(runname.c_str()))
        timedrift.Write();
    tpcdrift.Close();
    //std::abort(); // "Fixing" dirty cleanup from MINOS
}

uint getset(TArtEventStore &estore){
    const char *run = estore.GetRunInfo()->GetRunName()->Data();
    if(!strcmp(run,(const char *)"psp14")) return 0;
    else if(!strcmp(run,(const char *)"psp15")) return 1;
    else __throw_invalid_argument(Form("Unknown Run %s\n", run));
}