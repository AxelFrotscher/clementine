#include <Rtypes.h>
#include "MakeAllTree_78Ni.hh"

#include "progress.h"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtRecoPID.hh"
#include "TArtCalibDALI.hh"
#include "TVectorD.h"
#include "TFile.h"
#include "TArtPPAC.hh"
#include "TArtCalibPPAC.hh"
#include "TArtEventInfo.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtFocalPlane.hh"
#include <iostream>
#include "TArtDALIParameters.hh"
#include "TArtMINOSMap.hh"
#include "TArtMINOSPara.hh"
#include "TArtMINOSParameters.hh"

#include "TArtCalibMINOS.hh"
#include "TArtAnalyzedMINOS.hh"
#include "TArtTrackMINOS.hh"
#include "TArtVertexMINOS.hh"


R__LOAD_LIBRARY(libanacore.so)

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
    auto dpara = TArtDALIParameters::Instance();
    dpara->LoadParameter((char *)"config/db/DALI.xml");

    // Create CalibDALI to get and calibrate raw data
    auto dalicalib = new TArtCalibDALI();
    dalicalib->LoadFile((char *)"config/db/DALI");

    // Create Minos
    TArtMINOSParameters setup("MINOSParameters", "MINOSParameters");
    setup.LoadParameters((char*)"config/db/MINOS.xml");

    TArtCalibMINOS minoscalib;
    TArtAnalyzedMINOS minosanalyzed;
    TArtTrackMINOS minostrack;
    TArtVertexMINOS minosvertex;

    // Create TTree and output file for it
    TFile fout(output.c_str(), "RECREATE");
    if(!fout.IsOpen()) __throw_invalid_argument("Could not open output file!\n");

    auto tree = new TTree("tree","tree");

    // Define data nodes which are supposed to be dumped to tree 
    // EventInfo is important for the fBit information to know the trigger!
    vector<string> datanodes{"EventInfo", "BigRIPSPPAC", "BigRIPSPlastic",
                             "BigRIPSIC","BigRIPSFocalPlane","BigRIPSRIPS",
                             "BigRIPSTOF", "BigRIPSBeam", "DALINaI",
                             "EventInfo_FBIT", "TrackMINOS"};

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
    /*TClonesArray fidata, dataresult;
    fidata.SetClass("TMinosClust");
    dataresult.SetClass("TMinosResult");
    tree->Branch("fitdata", &fidata);
    tree->Branch("dataresult", &dataresult); */

    int track_no, track_no_final, evtOrig, padsleft;
    tree->Branch("track_no",&track_no, "track_no/I");
    tree->Branch("track_no_final", &track_no_final, "track_no_final/I");
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

    tree->Branch("xtar",&fXTar,"fXTar/D");
    tree->Branch("ytar",&fYTar,"fYTar/D");

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
    Int_t dalimultwotime        = 0;
    Int_t dalimult              = 0;
    Int_t dalitimetruemult      = 0;
    Int_t dalimultthres         = 0;
    Int_t dalitimetruemultthres = 0;

    tree->Branch("dalimultwotime",&dalimultwotime,"dalimultwotime/I");
    tree->Branch("dalimult",&dalimult,"dalimult/I");
    tree->Branch("dalitimetruemult",&dalitimetruemult,"dalitimetruemult/I");
    tree->Branch("dalimultthres",&dalimultthres,"dalimultthres/I");

    tree->Branch("dalitimetruemultthres",&dalitimetruemultthres,
                 "dalitimetruemultthres/I");

    const vector<string> ppacname{
        "F8PPAC-1A", "F8PPAC-1B", "F8PPAC-2A", "F8PPAC-2B",
        "F9PPAC-1A", "F9PPAC-1B", "F9PPAC-2A", "F9PPAC-2B", 
        "F11PPAC-1A","F11PPAC-1B","F11PPAC-2A","F11PPAC-2B"};

    // Progress Bar setup
    int neve = 0; // counting variable
    uint totevents = 5000000;
    const int downscale = 5000; // every n-th event

    progressbar progress(totevents, 0);

    while(estore.GetNextEvent() && (neve < totevents)){ //&& neve < 100000
        //Making the BigRIPS tree calibration
        brcalib->ClearData();
        brcalib->ReconstructData();

        //Reconstructiong the PID
        recopid.ClearData();
        recopid.ReconstructData();

        //Reconstructing the scattering angle:
        //Beam spot on target reconstruction
        fXTar = fYTar = magnum;
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
        minosvertex.ClearData();
        minoscalib.ReconstructData();
        //minosanalyzed.ReconstructData();
        //minostrack.ReconstructData();
        minosvertex.ReconstructVertex();
        // Convert the reconstructed vertex in z (mm)
        z_vertex = (minosvertex.GetZv()*TimeBinElec - DelayTrigger)*1e-3*VDrift*10;
        //time bin*(ns/us)*vdrift(cm/us)*(mm/cm)


        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Making DALI
        dalicalib->ClearData();
        //dalicalib->SetPlTime(plasticcalib->FindPlastic("F8pl")->GetTime());
        dalicalib->SetVertex(z_vertex);   // WE DID IT... Your turn to make it
                                        // compile :)
        //Add above to remove F8plastic tof.
        dalicalib->ReconstructData();
        dalimultwotime = dalicalib->GetMultWithoutT();
        dalimult = dalicalib->GetMult();
        dalitimetruemult = dalicalib->GetTimeTrueMult();
        dalimultthres = dalicalib->GetMultThres();
        dalitimetruemultthres = dalicalib->GetTimeTrueMultThres();

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