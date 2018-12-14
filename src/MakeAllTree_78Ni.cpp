#include <Rtypes.h>
#include "MakeAllTree_78Ni.hh"
#include "progress.h"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtRecoPID.hh"
//#include "TArtCalibDALI.hh"
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
//#include "TArtDALIParameters.hh"
#include "TArtMINOSMap.hh"
#include "TArtMINOSPara.hh"
#include "TArtMINOSParameters.hh"

#include "TArtCalibMINOS.hh"
#include "TArtAnalyzedMINOS.hh"
#include "TArtTrackMINOS.hh"
#include "TArtVertexMINOS.hh"

#include "TMath.h"

R__LOAD_LIBRARY(libanacore.so)
R__LOAD_LIBRARY(libminos.so)

using namespace std;

void generatetree(const string infile, const string output){
    //  signal(SIGINT,stop_interrupt); // CTRL + C , interrupt
    cout << "Now in Estore -> " << infile << endl;

    // Create StoreManager both for calibration "TArtCalib..." and 
    // treatment "TArtReco..."
    auto * sman = TArtStoreManager::Instance();
    // Create EventStore to control the loop and get the EventInfo
    TArtEventStore estore;
    //estore.SetInterrupt(&stoploop);
    if(!estore.Open(infile.c_str())) __throw_invalid_argument(("Could not open"+
                                                               infile).c_str());

    // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes 
    // parameters from ".xml" files
    TArtBigRIPSParameters para;
    vector<string> parameterfiles{
        "config/db/BigRIPSPPAC.xml", "config/db/BigRIPSPlastic.xml",
        "config/db/BigRIPSIC.xml",   "config/db/FocalPlane.xml"};

    for(auto i : parameterfiles){
        if(!para.LoadParameter(&i[0])) __throw_invalid_argument(&i[0]);
    }

    para.SetFocusPosOffset(8,138.5);
    // Create CalibPID to get and calibrate raw data ( CalibPID -> [CalibPPAC , 
    // CalibIC, CalibPlastic , CalibFocalPlane] 
    auto brcalib      = new TArtCalibPID();
    auto ppaccalib    = brcalib->GetCalibPPAC();
    auto cfpl         = brcalib->GetCalibFocalPlane();
    //auto *plasticcalib = brcalib->GetCalibPlastic();
    //auto *iccalib      = brcalib->GetCalibIC();

    // Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z, ... 
    // (RecoPID -> [ RecoTOF , RecoRIPS , RecoBeam] )
    TArtRecoPID recopid;

    // Definition of observables we want to reconstruct

    printf("Defining BigRIPS parameters\n");
    const vector<int> focalplanes{3, 5, 7, 8, 9, 11};

    vector<vector<string>> mfil{  // matrixfiles
            {"config/matrix/mat1.mat",               "D3"},    // F3 - F5 => D3
            {"config/matrix/mat2.mat",               "D5"},    // F5 - F7 => D5
            {"config/matrix/F8F9_LargeAccAchr.mat",  "D7"},    // F8 - F9 => D7
            {"config/matrix/F9F11_LargeAccAchr.mat", "D8"}};   // F9 - F11=> D8

    vector<TArtRIPS *> rips{
            recopid.DefineNewRIPS(3, 5,  &mfil[0][0][0], &mfil[0][1][0]),
            recopid.DefineNewRIPS(5, 7,  &mfil[1][0][0], &mfil[1][1][0]),
            recopid.DefineNewRIPS(8, 9,  &mfil[2][0][0], &mfil[2][1][0]),
            recopid.DefineNewRIPS(9, 11, &mfil[3][0][0], &mfil[3][1][0])};

    // Reconstruction of TOF DefineNewTOF(first plane,second plane, time offset)
    vector<vector<string>> fplname{{"F3pl", "F7pl"},
                                   {"F8pl", "F11pl-1"}};
    vector<double> tofoff{ //300.25 F3-F7 init -159.45 F8-F11 init
        304.17+0.17,   // good Offset Value for F3-F7,  empty-target run  300.85
        -161.64-0.08}; // good Offset Value for F8-F11, empty-target run -160.45

    vector<TArtTOF *> tof{
        recopid.DefineNewTOF(&fplname[0][0][0], &fplname[0][1][0], tofoff[0], 5),
        recopid.DefineNewTOF(&fplname[1][0][0], &fplname[1][1][0], tofoff[1], 9)};

    // Reconstruction of IC observables for ID
    vector<TArtBeam *> beam{  // br = BigRIPS, zd = ZeroDegree
        recopid.DefineNewBeam(rips[0], rips[1], tof[0], (char *) "F7IC"),  //br_37
        recopid.DefineNewBeam(rips[1], tof[0], (char *) "F7IC"),           //br_57
        recopid.DefineNewBeam(rips[2], tof[1], (char *) "F11IC"),          //zd_89
        recopid.DefineNewBeam(rips[3], tof[1], (char *) "F11IC"),          //zd_911
        recopid.DefineNewBeam(rips[2], rips[3], tof[1], (char *) "F11IC")};//zd_811

    // Create DALIParameters to get ".xml"
    //auto dpara = TArtDALIParameters::Instance();
    //dpara->LoadParameter((char *)"config/db/DALI.xml");

    // Create CalibDALI to get and calibrate raw data
    //auto dalicalib = new TArtCalibDALI();
    //dalicalib->LoadFile((char *)"config/db/DALI");

    // Create Minos
    TArtMINOSParameters setup("MINOSParameters", "MINOSParameters");
    setup.LoadParameters((char*)"config/db/MINOS.xml");

    TArtCalibMINOS minoscalib;
    //TArtAnalyzedMINOS minosanalyzed;
    //TArtTrackMINOS minostrack;
    //TArtVertexMINOS minosvertex;

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

    // Create output values for minos by hand -.-
    TMinosClust fitdata;
    TMinosResult dataresult;

    tree->Branch("MinosClustX", &fitdata.x_mm);
    tree->Branch("MinosClustY", &fitdata.y_mm);
    tree->Branch("MinosClustZ", &fitdata.z_mm);
    tree->Branch("MinosClustT", &fitdata.t_ns);
    tree->Branch("MinosClustQ", &fitdata.Chargemax);
    tree->Branch("MinosClustNCl", &fitdata.n_Cluster);
    tree->Branch("MinosClustNP", &fitdata.n_Pads);
    tree->Branch("MinosClustChi2", &fitdata.Chi2);

    tree->Branch("MinosResultX", &dataresult.x_mm);
    tree->Branch("MinosResultY", &dataresult.y_mm);
    tree->Branch("MinosResultZ", &dataresult.z_mm);
    tree->Branch("MinosResultQ", &dataresult.Chargemax);
    tree->Branch("MinosResultNCl", &dataresult.n_Cluster);
    tree->Branch("MinosResultNP", &dataresult.n_Pads);
    tree->Branch("MinosResultChi2", &dataresult.Chi2);

    int trackNbr =0, trackNbr_FINAL = 0, evtOrig = 0;
    tree->Branch("trackNbr",&trackNbr, "trackNbr/I");
    tree->Branch("trackNbr_final", &trackNbr_FINAL, "trackNbr_final/I");
    tree->Branch("evtOrig", &evtOrig, "evtOrig/I");

    double z_vertex =0, x_vertex=0, y_vertex =0, r_vertex =0, phi_vertex=0,
           thetaz1=0, thetaz2=0,
           VDrift,          // in cm/(1E-6*s)
           DelayTrigger,    // in ns
           Stoptime;
    tree->Branch("VDrift", &VDrift, "MINOS.VDrift/D");
    tree->Branch("DelayTrig", &DelayTrigger, "DelayTrigger/D");
    tree->Branch("x_vertex", &x_vertex, "x_vertex/D");
    tree->Branch("y_vertex", &y_vertex, "y_vertex/D");
    tree->Branch("z_vertex", &z_vertex, "z_vertex/D");
    tree->Branch("r_vertex", &r_vertex, "r_vertex/D");
    tree->Branch("phi_vertex", &phi_vertex, "phi_vertex/D");
    tree->Branch("thetaz1", &thetaz1, "thetaz1/D");
    tree->Branch("thetaz2", &thetaz1, "thetaz2/D");

    //Parameters for MINOS analysis
    double minosthresh, TimeBinElec, Tshaping;
    ifstream minosconfigfile;
    minosconfigfile.open("config/db/ConfigMINOSDrift_fixed.txt");
    if(!minosconfigfile.is_open())__throw_invalid_argument("Could not open "
                                        "config/db/ConfigMINOSDrift_fixed.txt");
    string trash;
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

    //BigRIPS
    const double magnum = -999; // Magic Number for unset events
    const int planesperppac = 4;
    vector<vector<double>> fF8PPAC(3, vector<double>(planesperppac, magnum));

    Double_t fXTar = magnum, fYTar = magnum;
 
    // 6 focal planes (3,5,7,8,9,11) and 4 value (X,A,Y,B)
    vector<vector<double>> focal{6, vector<double>(4)}; // F3X, F3A

    vector<int> fGoodPPACFocus(12, 0);
    vector<int> fGoodPPACFocusOr(12, 0);

    // F8 PPAC positions relative to F8 0: x 1: y
    const vector<vector<double>> f8z{{-1310.7,-1302.1},{-1273.3,-1281.9},
                                     {-810.7,-802.1},{-773.3,-781.9}};

    /*tree->Branch("xtar",&fXTar,"fXTar/D");
    tree->Branch("ytar",&fYTar,"fYTar/D");*/

    vector<vector<string>> focalname{{"F3X","F3A"}, {"F5X","F5A"}, 
                                     {"F7X","F7A"}, {"F8X","F8A"},
                                     {"F9X","F9A"}, {"F11X","F11A"}};

    for(uint i=0; i<focal.size(); i++){
        tree->Branch(focalname[i][0].c_str(), &focal.at(i).at(0));
        tree->Branch(focalname[i][1].c_str(), &focal.at(i).at(1));
    }

    //TArtFocalPlane *tfpl;
    TVectorD *vec;

    tree->Branch("fgoodppacfocus",fGoodPPACFocus.data(),"fGoodPPACFocus[12]/I");
    tree->Branch("fgoodppacfocusor",fGoodPPACFocusOr.data(),
                 "fGoodPPACFocusOr[12]/I");

    //DALI
    /*Int_t dalimultwotime        = 0;
    Int_t dalimult              = 0;
    Int_t dalitimetruemult      = 0;
    Int_t dalimultthres         = 0;
    Int_t dalitimetruemultthres = 0;

    tree->Branch("dalimultwotime",&dalimultwotime,"dalimultwotime/I");
    tree->Branch("dalimult",&dalimult,"dalimult/I");
    tree->Branch("dalitimetruemult",&dalitimetruemult,"dalitimetruemult/I");
    tree->Branch("dalimultthres",&dalimultthres,"dalimultthres/I");

    tree->Branch("dalitimetruemultthres",&dalitimetruemultthres,
                 "dalitimetruemultthres/I");*/

    const vector<string> ppacname{
        "F8PPAC-1A", "F8PPAC-1B", "F8PPAC-2A", "F8PPAC-2B",
        "F9PPAC-1A", "F9PPAC-1B", "F9PPAC-2A", "F9PPAC-2B", 
        "F11PPAC-1A","F11PPAC-1B","F11PPAC-2A","F11PPAC-2B"};

    // Progress Bar setup
    int neve = 0; // counting variable
    uint totevents = 50000; //50000
    const int downscale = 50; // every n-th event

    progressbar progress(totevents, 0);

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

        /*
        //Double_t modif = - 0.0011*(F11X-8) - 0.0000003*(F11X-8)*(F11X-8)*(F11X-8)
        //                 - 0.000023*(F9X-10) - 0.0000005*(F9X)*(F9X)
        //                 + 0.000000005*(F9X+10)*(F9X+10)*(F9X+10);

        //AoQ_Z_BR->Fill(beam_br_57->GetAoQ(),beam_br_57->GetZet());
        //AoQ_Z_ZD->Fill(beam_zd_911->GetAoQ()+modif,beam_zd_911->GetZet());*/

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Making MINOS
        minoscalib.ClearData();
        //minosanalyzed.ClearData();
        //minostrack.ClearData();
        //minosvertex.ClearData();
        minoscalib.ReconstructData();
        //minosanalyzed.ReconstructData();
        //minostrack.ReconstructData();
        //minosvertex.ReconstructVertex();
        // Convert the reconstructed vertex in z (mm)
        //z_vertex = (minosvertex.GetZv()*TimeBinElec - DelayTrigger)*1e-3*VDrift*10;
        //time bin*(ns/us)*vdrift(cm/us)*(mm/cm)

        // MINOS 1. Filling vectors with (x,y) data ////////////////////////////
        vector<double> Xpad, Ypad, Qpad, Xpadnew,Ypadnew, Qpadnew, Zpadnew;
        int filled = 0; // Number of tracks
        fitdata.reset();
        TArtCalibMINOSData *minos = new TArtCalibMINOSData;
        for(int i=0; i<minoscalib.GetNumCalibMINOS(); i++){
            minos = minoscalib.GetCalibMINOS(i);
            double maxcharge =0;
            double x_mm = minos->GetX();
            double y_mm = minos->GetY();
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

        // MINOS 2. Modify with hough-transformation
        int padsleft = Xpad.size();
        vector<bool> clusterringbool;
        vector<int> clusternbr, clusterpads;
        trackNbr = trackNbr_FINAL = 0;

        if(filled){
            for(int iteration =0; (Xpad.size()>=10 && iteration <20); iteration++){
                int filter_result = Obertelli_filter(Xpad, Ypad, Qpad, Xpadnew,
                                                     Ypadnew, Qpadnew,
                                                     clusterringbool);
            for(int i=0; i<filter_result; i++){
                clusternbr.push_back(iteration);
                clusterpads.push_back(filter_result);
            }
            if(filter_result > 10 && clusterringbool.back()) trackNbr++;
            }

        }

        for(int i=0; i<Xpadnew.size(); i++){
            fitdata.add(Xpadnew[i],Ypadnew[i], -1E4,-1E4,Qpadnew[i],
                        clusternbr[i], clusterpads[i],0);
            Zpadnew.push_back(-1E4);
        }

        // Bonus for Tracknumber between 1 and 4
        TH1F hfit("hfit", "hfit", 512,0,512);
        auto fit_function = new TF1("fit_function",conv_fit, 0,511,3);
        if(trackNbr > 0 && trackNbr < 4){
            if(!filled) cerr << "Tracknumber:" << trackNbr << " but no evts." << endl;


            // MINOS 3: Fitting the taken pads for Qmax and Ttrig information //
            padsleft -= Xpadnew.size();
            for(int i=0; i<minoscalib.GetNumCalibMINOS(); i++){
                minos = minoscalib.GetCalibMINOS(i);
                double x_mm = minos->GetX();
                double y_mm = minos->GetY();
                hfit.Reset(); // reset fitting histogram
                bool fitbool = false;
                int indexfill = 0;

                for(int j=0; j < Xpadnew.size(); j++){
                    if(abs(Xpadnew[j]-x_mm) <0.01 && abs(Ypadnew[j]-y_mm)<0.01){
                        fitbool = true;
                        indexfill = j;
                        break;
                    }
                }
                // Check if new channel is of interest
                if(fitbool){
                    for(int k=0; k<minos->GetNData(); k++){
                        if(minos->GetCalibValue(k) >= 0)
                            hfit.SetBinContent(hfit.FindBin(minos->GetCalibValue(k)),
                                                 minos->GetCalibValue(k) + 250);
                    }
                    // Fitting the hfit histogram of last ch. if not empty
                    if(hfit.GetSumOfWeights()>0){
                        hfit.GetXaxis()->SetRange(0,510);
                        double hfit_max = hfit.GetMaximum();
                        double hfit_max_T = hfit.GetMaximumBin();

                        // Find T_min and T_max limits for non-0 signals
                        double T_min = hfit.FindLastBinAbove(250);
                        double T_max = hfit.FindFirstBinAbove(0)-1;

                        // Take only 1.5*shaping time before max if other signals before
                        if(hfit_max_T-3.5*(Tshaping/TimeBinElec) > T_min)
                            T_min = hfit_max_T-2*Tshaping/TimeBinElec;
                        if((hfit_max_T +10)< T_max || T_max == -1) T_max = hfit_max_T+10;

                        T_min = std::max(T_min, 0.);
                        T_max = std::min(T_max, 510.);

                        // Set fit parameters
                        fit_function->SetParameters(hfit_max-250,
                                                   hfit_max_T- Tshaping/TimeBinElec,
                                                   Tshaping/TimeBinElec);
                        fit_function->SetParLimits(0,0,1E5);
                        fit_function->SetParLimits(1,-20,512);
                        fit_function->SetParLimits(2,0,512);

                        int fit2DStatus = hfit.Fit(fit_function, "Q", "", T_min,T_max);
                        double fit_function_max =0, fit_function_Tpad =0, Chi2 =0;

                        if(!fit2DStatus){
                            auto fit_result = hfit.GetFunction("fit_function");
                            Chi2 = fit_result->GetChisquare();
                            fit_function_max = fit_function->GetMaximum();
                            fit_function_Tpad = fit_function->GetParameter(1);
                        }

                        // attribute q_pad and z_mm value
                        double q_pad =0, z_mm=0, t_pad=0;
                        if(fit2DStatus || fit_function_max <=20 ||
                           fit_function_max > 1E5 || fit_function_Tpad < .15 ||
                           fit_function_Tpad > 512 || fit_function->GetParameter(2) < .15 ||
                           fit_function->GetParameter(2) > 512){
                            // No correct fit here :(
                            q_pad = hfit_max -250;
                            z_mm = -1E4;
                        }
                        else{ // Add to the variables the fit parameter
                            t_pad = fit_function_Tpad;
                            z_mm = (t_pad*TimeBinElec-DelayTrigger)*VDrift;
                            q_pad = fit_function_max-250;
                        }

                        fitdata.replace(Xpadnew[indexfill], Ypadnew[indexfill],
                                        t_pad*TimeBinElec, z_mm, q_pad,
                                        clusternbr[indexfill], clusterpads[indexfill],
                                        Chi2, indexfill);
                        Zpadnew[indexfill] = z_mm;
                        Qpadnew[indexfill] = q_pad;

                    } // End if histogram not empty
                } // end if fitbool true
                else continue;
            } // End of all entries

            // 4: MINOS Filtering the tracks off possible noise with Hough3D ///
            int padsleft2 = Xpadnew.size();
            vector<double> xin, yin, zin, qin, xout, yout, zout, qout;
            int cluster_temp = 0, ringsum =0, npoint1=0, npoint2=0, array_final = 0;
            double zmax =0 ;
            vector<int> ringtouch(18,0);
            auto gryz_1 = new TGraph;
            auto grxz_1 = new TGraph;
            auto gryz_2 = new TGraph;
            auto grxz_2 = new TGraph;
            for(int i=0; i<padsleft2;i++){
                if(xin.size() && ((cluster_temp !=int(clusternbr[i]) && !i) ||
                   i==(Xpadnew.size()-1))){
                    Hough_filter(xin,yin,zin,qin,xout,yout,zout,qout);

                    for(int k=0; k<xout.size(); k++){
                        zmax = max(zmax,zout[k]);
                        ringtouch.at(int((sqrt(pow(xout[k],2.) +
                                               pow(yout[k],2.))-45.2)/2.1))++;
                    }
                    for(auto &j: ringtouch) if(j>0) ringsum++;
                    if(zmax>290) ringsum=16;
                    if(xout.size()>10 && ringsum >= 15) {
                        cout << xout.size() <<" "<< ringsum << " Padsleft2" << endl;
                        trackNbr_FINAL++;
                        //cluster1 = cluster_temp; delete if no use
                        for (int l = 0; l < xout.size(); l++) {
                            grxz_1->SetPoint(npoint1, zout[l], xout[l]);
                            gryz_1->SetPoint(npoint1, zout[l], yout[l]);
                            dataresult.add(xout[l], yout[l], zout[l],
                                           qout[l], trackNbr_FINAL, xout.size(), zmax);
                            array_final++;
                            npoint2++;
                        }
                    }
                    xin.clear();
                    yin.clear();
                    zin.clear();
                    qin.clear();
                    xout.clear();
                    yout.clear();
                    zout.clear();
                    qout.clear();
                    //npoint_temp =0; useless?
                    ringsum = 0;
                    zmax =0;
                    for(auto &i:ringtouch) i=0;
                }

                cluster_temp = clusternbr[i];
                if(!(clusterpads[i]>=10 && clusterringbool[i] == true &&
                     Zpadnew[i]> -1E4 && Zpadnew[i]<=320)) continue;
                else{
                    xin.push_back(Xpadnew[i]);
                    yin.push_back(Ypadnew[i]);
                    zin.push_back(Zpadnew[i]);
                    qin.push_back(Qpadnew[i]);
                    //npoint_temp++;
                }
            } // end of loop on pads
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Making DALI
        /*dalicalib->ClearData();
        //dalicalib->SetPlTime(plasticcalib->FindPlastic("F8pl")->GetTime());
        dalicalib->SetVertex(z_vertex);   // WE DID IT... Your turn to make it
                                        // compile :)
        //Add above to remove F8plastic tof.
        dalicalib->ReconstructData();
        dalimultwotime = dalicalib->GetMultWithoutT();
        dalimult = dalicalib->GetMult();
        dalitimetruemult = dalicalib->GetTimeTrueMult();
        dalimultthres = dalicalib->GetMultThres();
        dalitimetruemultthres = dalicalib->GetTimeTrueMultThres(); */

        tree->Fill();
        neve++;

        progress.increaseevent();
        if(!(neve%downscale)) progress.draw();

        // Add Graphical Feedback
    }
    progress.reset();
    cout << "Writing the tree..."<<endl;

    fout.Write();
    fout.Close();

    printf("Writing TTree complete!\n");

}

int Obertelli_filter(vector<double> &x,vector<double> &y,vector<double> &q,
                     vector<double> &x_out,vector<double> &y_out,
                     vector<double> &q_out, vector<bool> &ringbool) {
    double bint1=2.;
    double bint2=2.;
    int maxt = 360;
    int mint = 0;
    int nt1 = (maxt-mint)/bint1;
    int nt2 = (maxt-mint)/bint1;

    double PI = TMath::Pi();
    double Rint = 45.2;
    double Rext = 45.2 + 18*2.1;
    int filter_result = 0;

    TH2F hp_xy("hp_xy", "hp_xy", nt1, mint, maxt, nt2, mint, maxt);
    TH2F hpDiag_xy("hpDiag_xy", "hpDiag_xy", nt1, mint, maxt, nt2, mint, maxt);
//	TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-85,85,100,-85,85);
//	TH2F *hcnew_xy = new TH2F("hcnew_xy","Cluster in xy plane AFTER Obertelli"
//                                        "transform",100,-85,85,100,-85,85);

    double max_xy;
//	TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);

    vector<double> xTemp, yTemp, qTemp;

    double theta1, theta2 =0, xt1, yt1, xt2, yt2,
           line0=0., line1=0.,
           delta=0., AA=0., BB=0., CC=0.,
           maxtheta1=0., maxtheta2=0., xmax1=0., ymax1=0., xmax2=0., ymax2=0.,
           par0=0., par1=0.,
           r_mm=0.;
    int ringsum=0;
    bool maxfound = false;

    for(unsigned int i=0;i<x.size();i++){
        xTemp.push_back(x.at(i));
        yTemp.push_back(y.at(i));
        qTemp.push_back(q.at(i));

        //Fill coordinate space histograms for plots
//		hc_xy->Fill(x->at(i),y->at(i),q->at(i));

        //Loop of indices
        for(int j=0; j<nt1; j++){
            theta1 = (j+0.5)*bint1 + mint;
            xt1 = Rint * TMath::Cos(theta1*PI/180.);
            yt1 = Rint * TMath::Sin(theta1*PI/180.);
            line1 = (yt1 - y.at(i))/(xt1 - x.at(i));
            line0 = yt1 - xt1 * line1;
            AA = 1 + line1*line1;
            BB = 2*line0*line1;
            CC = line0*line0 - Rext*Rext;

            delta = BB*BB - 4*AA*CC;
            if(delta >= 0){
                xt2 = -(BB + sqrt(delta))/(2*AA);
                yt2 = line0 + line1*xt2;
                if(xt2 <= 0)	        theta2 = 180 - asin(yt2/Rext)*180/PI;
                else{
                    if     (yt2 >  0)   theta2 = asin(yt2/Rext)*180/PI;
                    else if(yt2 <= 0)	theta2 = 360 + asin(yt2/Rext)*180/PI;
                }

                //if(yt2>0){theta2 = 180./PI*acos(xt2/Rext);}
                //else{theta2=360. - 180./PI*acos(xt2/Rext);}

                if((xt1*x.at(i) + yt1*y.at(i))>=0 &&
                   (xt2*x.at(i) + yt2*y.at(i))>=0 && (xt1*xt2+yt1*yt2)>=0){
                    hp_xy.Fill(theta1,theta2);
                    if(abs(theta1-theta2) <= 10) hpDiag_xy.Fill(theta1,theta2);
                }
                else{
                    if(delta != 0){
                        xt2 = (-BB + sqrt(delta))/(2*AA);
                        yt2 = line0 + line1*xt2;
                        if(xt2 <= 0)	      theta2 = 180 - asin(yt2/Rext)*180/PI;
                        else if(xt2 > 0){
                            if(yt2 > 0)	      theta2 = asin(yt2/Rext)*180/PI;
                            else if(yt2 <=0 ) theta2 = 360 + asin(yt2/Rext)*180/PI;
                        }
                        //if(yt2>0){theta2 = 180./PI*acos(xt2/Rext);}
                        //else{theta2=360. - 180./PI*acos(xt2/Rext);}
                        if( (xt1*x.at(i) + yt1*y.at(i)) >= 0 &&
                            (xt2*x.at(i) + yt2*y.at(i)) >= 0 &&
                            (xt1*xt2+yt1*yt2) >= 0){
                            hp_xy.Fill(theta1,theta2);
                            if(abs(theta1-theta2) <= 10) hpDiag_xy.Fill(theta1,theta2);
                        }
                    }
                }
            }
        }
    }
    x.clear();
    y.clear();
    q.clear();

    if(hpDiag_xy.GetMaximum() >= 10) max_xy = hpDiag_xy.GetMaximum();
//		cout << "Max taken in diag... withh value=" << max_xy << endl;
    else max_xy = hp_xy.GetMaximum();

    for(int ii=0; ii<nt1; ii++){
        if(maxfound) break;
        for(int jj=0; jj<nt2; jj++){
            if(hp_xy.GetBinContent(ii+1, jj+1) == max_xy){
                maxtheta1 = (ii+0.5)*bint1 + mint;
                maxtheta2 = (jj+0.5)*bint2 + mint;
                maxfound = true;
                //cout << "xy: theta max are " << maxtheta1 << " , " << maxtheta2 << endl;

            }
            if(maxfound) break;
        }
    }

    xmax1 = Rint * TMath::Cos(maxtheta1*PI/180.);
    ymax1 = Rint * TMath::Sin(maxtheta1*PI/180.);
    xmax2 = Rext * TMath::Cos(maxtheta2*PI/180.);
    ymax2 = Rext * TMath::Sin(maxtheta2*PI/180.);

    // xy PEAK
    par1 = (ymax2-ymax1)/(xmax2-xmax1);
    par0 = (ymax1 - xmax1*par1);

    /*cout<<"xmax1 "<<xmax1<<" ymax1 "<<ymax1<<" xmax2 "<<xmax2<<" ymax2 "<<ymax2<<endl;
    line_xy->SetParameter(0,par0);
    line_xy->SetParameter(1,par1);
    hc_xy->GetListOfFunctions()->Add(line_xy);
	line_xy->SetLineWidth(1);*/

    //Selection of x,y points IN the maxmean+/-1 found in Obertelli transform of xy plane
    for(unsigned int i=0; i<xTemp.size(); i++){
        if( (abs(par1*xTemp[i]-yTemp[i]+par0)/sqrt(1+par1*par1))<= 6 &&
            ((xmax1*xTemp[i] + ymax1*yTemp[i]) >= 0) &&
            ((xmax2*xTemp[i] + ymax2*yTemp[i]) >= 0) &&
            ((xmax1*xmax2 + ymax1*ymax2) >= 0)){
			/*cout << "Taken points= " << xTemp[i] << " , " << yTemp[i] << " , " << zTemp[i] << endl;
			hcnew_xy->Fill(xTemp[i],yTemp[i],qTemp[i]);*/
            x_out.push_back(xTemp[i]);
            y_out.push_back(yTemp[i]);
            q_out.push_back(qTemp[i]);
            filter_result++;
            r_mm = sqrt(xTemp[i]*xTemp[i]+yTemp[i]*yTemp[i]);
            if(r_mm < (45.2+5*2.1)) ringsum++;
        }
        else{
            x.push_back(xTemp[i]);
            y.push_back(yTemp[i]);
            q.push_back(qTemp[i]);
        }
    }

    for(int ip=0; ip<filter_result; ip++){
        if(ringsum>2) ringbool.push_back(true);
        else ringbool.push_back(false);
    }
/*
	c1->Divide(3,1);
	// Coordinate space
	c1->cd(1);
	hc_xy->Draw("colz");
	// Hough space
	c1->cd(2);
	hp_xy->Draw("colz");
	// Coordinate space : New plot
	c1->cd(3);
	hcnew_xy->Draw("colz");

	c1->Update();
*/

    return filter_result;
}

double conv_fit(double *x, double *p) {
    // Check for boundaries of x
    if(x[0]<p[1] || x[0] > 512) return 250;
    else return p[0] * exp(-3.*(x[0]-p[1])/p[2]) * sin((x[0]-p[1])/p[2]) *
                pow((x[0]-p[1])/p[2], 3) + 250;
}


void Hough_filter(vector<double> &x,vector<double> &y,vector<double> &z,
                  vector<double> &q,vector<double> &x_out,vector<double> &y_out,
                  vector<double> &z_out,vector<double> &q_out) {
    int nt_xy = 180, nt_xz = 180, nt_yz = 180,
        nr_xy = 45,  nr_xz = 300, nr_yz = 300;
    double bint_xy = 2., bint_xz = 2., bint_yz = 2.,
           binr_xy = 3., binr_xz = 3., binr_yz = 3.;
    int nt = nt_xy,nr = nr_xy;
    double PI = TMath::Pi();

    double rho_xy, rho_xz, rho_yz;
    double theta_xy, theta_xz, theta_yz;

    TH2F hp_xy("hp_xy", "hp_xy", nt_xy, 0, 180, nr_xy, -1*nr_xy, nr_xy);
    TH2F hp_xz("hp_xz", "hp_xz", nt_xz, 0, 180, nr_xz, -1*nr_xz, nr_xz);
    TH2F hp_yz("hp_yz", "hp_yz", nt_yz, 0, 180, nr_yz, -1*nr_yz, nr_yz);

    /*TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-85,85,100,-85,85);
    TH2F *hc_xz = new TH2F("hc_xz","Track in xz plane",250,-50,450,100,-85,85);
    TH2F *hc_yz = new TH2F("hc_yz","Track in yz plane",250,-50,450,100,-85,85);

   	TH2F *hcnew_xy = new TH2F("hcnew_xy","Track in xy plane AFTER Hough transform",100,-85,85,100,-85,85);
   	TH2F *hcnew_xz = new TH2F("hcnew_xz","Track in xz plane AFTER Hough transform",250,-50,450,100,-85,85);
   	TH2F *hcnew_yz = new TH2F("hcnew_yz","Track in yz plane AFTER Hough transform",250,-50,450,100,-85,85);*/

    //	int npeaks_xy, npeaks_xz, npeaks_yz;
    vector<double> thetapeaks_xy, rpeaks_xy, thetapeaks_xz, rpeaks_xz,
                   thetapeaks_yz, rpeaks_yz;
    double max_xy, max_xz, max_yz, rmean_xy=0, thetamean_xy=0, rmean_xz=0,
           thetamean_xz=0, rmean_yz=0, thetamean_yz=0;
    /*TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);
    TF1* line_xz = new TF1("line_xz","[0] + [1]*x",-50,450);
    TF1* line_yz = new TF1("line_yz","[0] + [1]*x",-50,450);*/

    double r0_xy=0., r0_xz=0., r0_yz=0., rmin_xy=0., rmin_xz=0., rmin_yz=0.,
           rmax_xy=0., rmax_xz=0., rmax_yz=0.,
           tmin=0., tmax=0.,
           rinf=0., rsup=0.;

    if(nt<nt_xz)nt=nt_xz;
    if(nr<nr_xz)nr=nr_xz;
    if(nt<nt_yz)nt=nt_yz;
    if(nr<nr_yz)nr=nr_yz;

    //cout<<"In filter !!!"<< endl << "size: "<<x.size()<<endl;

    for(unsigned int i=0; i < x.size(); i++){
        //Fill coordinate space histograms for plots
		/*hc_xy.Fill(x.at(i),y.at(i),q.at(i));
		hc_xz.Fill(z.at(i),x.at(i),q.at(i));
		hc_yz.Fill(z.at(i),y.at(i),q.at(i));*/

        //Loop of indices and fill Histograms
        for(int j=0; j < nt; j++){
            //xy
            theta_xy = j * 180./nt_xy;
            rho_xy = x.at(i) * TMath::Cos(theta_xy*PI/180.) +
                     y.at(i) * TMath::Sin(theta_xy*PI/180.);
            if(abs(theta_xy) < 180. && abs(rho_xy) < nr_xy){
                //if(i%40==0) cout<<"i="<<i<<" xy "<<rho_xy<<" "<<theta_xy<<endl;
                hp_xy.Fill(theta_xy, rho_xy);
            }

            //xz
            theta_xz = j * 180./nt_xz;
            rho_xz = z.at(i) * TMath::Cos(theta_xz*PI/180.) +
                     x.at(i) * TMath::Sin(theta_xz*PI/180.);
            if(abs(theta_xz) < 180. && abs(rho_xz) < nr_xz){
                //if(i%40==0) cout<<"i="<<i<<" xz "<<rho_xz<<" "<<theta_xz<<endl;
                hp_xz.Fill(theta_xz, rho_xz);
            }

            //yz
            theta_yz = j * 180./nt_yz;
            rho_yz = z.at(i) * TMath::Cos(theta_yz*PI/180.) +
                     y.at(i) * TMath::Sin(theta_yz*PI/180.);
            if(abs(theta_yz) < 180. && abs(rho_yz) < nr_yz) {
                //if(i%40==0) cout<<"i="<<i<<" yz "<<rho_yz<<" "<<theta_yz<<endl;
                hp_yz.Fill(theta_yz, rho_yz);
            }
        }
    }

    max_xy = hp_xy.GetMaximum();
    max_xz = hp_xz.GetMaximum();
    max_yz = hp_yz.GetMaximum();

    for(int ii=0; ii<nt; ii++){
        for(int jj=0; jj<nr; jj++){
            if(hp_xy.GetBinContent(ii+1, jj+1) == max_xy && jj<nr_xy){
                thetapeaks_xy.push_back((ii+0.5) * nt_xy/nt);
                rpeaks_xy.push_back((jj+0.5) * 2 - nr_xy);
                rmean_xy     += rpeaks_xy.back();
                thetamean_xy += thetapeaks_xy.back();
                // cout << "xy: " << thetapeaks_xy.back() << " , "
                //      << rpeaks_xy.back() << endl;
            }
            if(hp_xz.GetBinContent(ii+1, jj+1) == max_xz){
                thetapeaks_xz.push_back((ii+0.5) * nt_xz/nt);
                rpeaks_xz.push_back((jj+0.5)*2 - nr_xz);
                rmean_xz     += rpeaks_xz.back();
                thetamean_xz += thetapeaks_xz.back();
                // cout << "xz: " << thetapeaks_xz.back() << " , "
                //      << rpeaks_xz.back() << endl;
            }
            if(hp_yz.GetBinContent(ii+1, jj+1) == max_yz){
                thetapeaks_yz.push_back((ii+0.5)*nt_yz/nt);
                rpeaks_yz.push_back((jj+0.5)*2 - nr_yz);
                rmean_yz     += rpeaks_yz.back();
                thetamean_yz += thetapeaks_yz.back();
                // cout << "yz: " << thetapeaks_yz.back() << " , "
                //      << rpeaks_yz.back() << endl;
            }
        }
    }

    /*cout << "Number of max found :::     IN xy = " << rpeaks_xy.size()
         << " ,     IN xz = " << rpeaks_xz.size() << " ,     IN yz = "
         << rpeaks_yz.size() << endl;

    // xy PEAK
    rmean_xy = rmean_xy / rpeaks_xy.size();
    thetamean_xy = thetamean_xy / thetapeaks_xy.size();
    line_xy.SetParameter(0,rmean_xy / (TMath::Sin(thetamean_xy*PI/180)));
    line_xy.SetParameter(1,( -(TMath::Cos(thetamean_xy*PI/180)) /
                               (TMath::Sin(thetamean_xy*PI/180)) ));
    hc_xy.GetListOfFunctions()->Add(line_xy);

    // xz PEAK
    rmean_xz = rmean_xz / rpeaks_xz.size();
    thetamean_xz = thetamean_xz / thetapeaks_xz.size();
    line_xz.SetParameter(0,rmean_xz / (TMath::Sin(thetamean_xz*PI/180)));
    line_xz.SetParameter(1,( -(TMath::Cos(thetamean_xz*PI/180)) /
                                 (TMath::Sin(thetamean_xz*PI/180)) ));
    hc_xz.GetListOfFunctions()->Add(line_xz);

    // yz PEAK
    rmean_yz = rmean_yz / rpeaks_yz.size();
    thetamean_yz = thetamean_yz / thetapeaks_yz.size();
    line_yz.SetParameter(0,rmean_yz / (TMath::Sin(thetamean_yz*PI/180)));
    line_yz.SetParameter(1,( -(TMath::Cos(thetamean_yz*PI/180)) /
                               (TMath::Sin(thetamean_yz*PI/180)) ));
    hc_yz.GetListOfFunctions()->Add(line_yz); */

    rmean_xy = rpeaks_xy[0];
    thetamean_xy = thetapeaks_xy[0];
    rmean_xz = rpeaks_xz[0];
    thetamean_xz = thetapeaks_xz[0];
    rmean_yz = rpeaks_yz[0];
    thetamean_yz = thetapeaks_yz[0];

    /*line_xy.SetLineWidth(1);
    line_xz.SetLineWidth(1);
    line_yz.SetLineWidth(1);*/

    // Selection of x,y,z points COMMON to the 3 maxmean+/-1 found in Hough
    // spaces for xy, xz and yz spaces
    for(unsigned int i=0; i<x.size(); i++){
        r0_xy = x.at(i) * TMath::Cos(thetamean_xy * PI/180.) +
                y.at(i) * TMath::Sin(thetamean_xy * PI/180.);
        tmin = thetamean_xy - bint_xy;
        tmax = thetamean_xy + bint_xy;
        if((tmin) < 0)   tmin = tmin + 180.;
        if((tmax) > 180) tmax = tmax - 180.;
        rmin_xy = x.at(i) * TMath::Cos(tmin * PI/180.) +
                  y.at(i) * TMath::Sin(tmin * PI/180.);
        rmax_xy = x.at(i) * TMath::Cos(tmax * PI/180.) +
                  y.at(i) * TMath::Sin(tmax * PI/180.);

        rinf = min( rmean_xy - binr_xy, rmean_xy + binr_xy);
        rsup = max( rmean_xy - binr_xy, rmean_xy + binr_xy);
        if((r0_xy >= rinf || rmin_xy >= rinf || rmax_xy >= rinf) &&
           (r0_xy <= rsup || rmin_xy <= rsup || rmax_xy <= rsup)){
            r0_xz = z.at(i) * TMath::Cos(thetamean_xz * PI/180.) +
                    x.at(i) * TMath::Sin(thetamean_xz * PI/180.);
            tmin = thetamean_xz - bint_xz;
            tmax = thetamean_xz + bint_xz;
            if((tmin) < 0)   tmin = tmin + 180.;
            if((tmax) > 180) tmax = tmax - 180.;
            rmin_xz = z.at(i) * TMath::Cos(tmin * PI/180.) +
                      x.at(i) * TMath::Sin(tmin * PI/180.);
            rmax_xz = z.at(i) * TMath::Cos(tmax * PI/180.) +
                      x.at(i) * TMath::Sin(tmax * PI/180.);

            rinf = min( rmean_xz - binr_xz, rmean_xz + binr_xz);
            rsup = max( rmean_xz - binr_xz, rmean_xz + binr_xz);

            if((r0_xz >= rinf || rmin_xz >= rinf || rmax_xz >= rinf) &&
               (r0_xz <= rsup || rmin_xz <= rsup || rmax_xz <= rsup)){
                r0_yz = z.at(i) * TMath::Cos(thetamean_yz * PI/180.)+
                        y.at(i) * TMath::Sin(thetamean_yz * PI/180.);
                tmin = thetamean_yz - bint_yz;
                tmax = thetamean_yz + bint_yz;
                if((tmin) < 0)   tmin = tmin + 180.;
                if((tmax) > 180) tmax = tmax - 180.;
                rmin_yz = z.at(i) * TMath::Cos(tmin * PI/180.) +
                          y.at(i) * TMath::Sin(tmin * PI/180.);
                rmax_yz = z.at(i) * TMath::Cos(tmax * PI/180.) +
                          y.at(i) * TMath::Sin(tmax * PI/180.);

                rinf = min( rmean_yz - binr_yz, rmean_yz + binr_yz);
                rsup = max( rmean_yz - binr_yz, rmean_yz + binr_yz);

                if((r0_yz >= rinf || rmin_yz >= rinf || rmax_yz >= rinf) &&
                   (r0_yz <= rsup || rmin_yz <= rsup || rmax_yz <= rsup)){
                    // cout << "Taken points= " << x->at(i) << " , " << y->at(i)
                    //      << " , " << z->at(i) << endl;
                    // hcnew_xy->Fill(x->at(i),y->at(i),q->at(i));
                    // hcnew_xz->Fill(z->at(i),x->at(i),q->at(i));
                    // hcnew_yz->Fill(z->at(i),y->at(i),q->at(i));
                    x_out.push_back(x.at(i));
                    y_out.push_back(y.at(i));
                    z_out.push_back(z.at(i));
                    q_out.push_back(q.at(i));
                }
            }
        }
    }
    /*

     c1->Divide(3,3);
     // Coordinate space
     c1->cd(1);
     hc_xy->Draw("colz");
     c1->cd(2);
     hc_xz->Draw("colz");
     c1->cd(3);
     hc_yz->Draw("colz");

     // Hough space
     c1->cd(4);
     hp_xy->Draw("colz");
     c1->cd(5);
     hp_xz->Draw("colz");
     c1->cd(6);
     hp_yz->Draw("colz");

     // Coordinate space : New plots
     c1->cd(7);
     hcnew_xy->Draw("colz");
     c1->cd(8);
     hcnew_xz->Draw("colz");
     c1->cd(9);
     hcnew_yz->Draw("colz");

     c1->Update();
     */
}

