#include "MakeAllTree_78Ni.hh"
R__LOAD_LIBRARY(libanacore.so)

using namespace std;

// function to exit loop at keyboard interrupt. 
bool stoploop = false;
void stop_interrupt(){
    printf("keyboard interrupt\n");
    stoploop = true;
}

void generatetree(const string infile, const string output){
    //  signal(SIGINT,stop_interrupt); // CTRL + C , interrupt
    cout << "Now in Estore -> " << infile << endl;

    // Create StoreManager both for calibration "TArtCalib..." and 
    // treatment "TArtReco..."
    auto * sman = TArtStoreManager::Instance();
    // Create EventStore to control the loop and get the EventInfo
    TArtEventStore estore;
    estore.SetInterrupt(&stoploop); 
    estore.Open(infile.c_str());

    // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes 
    // parameters from ".xml" files
    TArtBigRIPSParameters para;
    vector<string> parameterfiles{
        "config/db/BigRIPSPPAC.xml", "config/db/BigRIPSPlastic.xml",
        "config/db/BigRIPSIC.xml",   "config/db/FocalPlane.xml"};

    for(auto i : parameterfiles) para.LoadParameter(&i[0]);

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
            {"config/matrix/mat1.mat",               "D3"},                // F3 - F5 => D3
            {"config/matrix/mat2.mat",               "D5"},                // F5 - F7 => D5
            {"config/matrix/F8F9_LargeAccAchr.mat",  "D7"},    // F8 - F9 => D7
            {"config/matrix/F9F11_LargeAccAchr.mat", "D8"}}; // F9 - F11  => D8

    vector<TArtRIPS *> rips{
            recopid.DefineNewRIPS(3, 5, &mfil[0][0][0], &mfil[0][1][0]),
            recopid.DefineNewRIPS(5, 7, &mfil[1][0][0], &mfil[1][1][0]),
            recopid.DefineNewRIPS(8, 9, &mfil[2][0][0], &mfil[2][1][0]),
            recopid.DefineNewRIPS(9, 11, &mfil[3][0][0], &mfil[3][1][0])};

    // Reconstruction of TOF DefineNewTOF(first plane,second plane, time offset)
    vector<vector<string>> fplname{{"F3pl", "F7pl"},
                                   {"F8pl", "F11pl-1"}};
    vector<double> tofoff{ //300.25 F3-F7 init -159.45 F8-F11 init
            304.50,   // good Offset Value for F3-F7,  empty-target run  300.85
            -162.49}; // good Offset Value for F8-F11, empty-target run -160.45

    vector<TArtTOF *> tof{
            recopid.DefineNewTOF(&fplname[0][0][0], &fplname[0][1][0], tofoff[0], 5),
            recopid.DefineNewTOF(&fplname[1][0][0], &fplname[1][1][0], tofoff[1], 9)};

    // Reconstruction of IC observables for ID
    vector<TArtBeam *> beam{  // br = BigRIPS, zd = ZeroDegree
            recopid.DefineNewBeam(rips[0], rips[1], tof[0], (char *) "F7IC"),   //br_37
            recopid.DefineNewBeam(rips[1], tof[0], (char *) "F7IC"),   //br_57
            recopid.DefineNewBeam(rips[2], tof[1], (char *) "F11IC"),  //zd_89
            recopid.DefineNewBeam(rips[3], tof[1], (char *) "F11IC"),  //zd_911
            recopid.DefineNewBeam(rips[2], rips[3], tof[1], (char *) "F11IC")}; //zd_811

    // To my knowledge, only [0] and [4] usable (tof-Focalplane mismatch)

    // Create DALIParameters to get ".xml"
    auto dpara = TArtDALIParameters::Instance();
    dpara->LoadParameter((char *)"config/db/DALI.xml");

    // Create CalibDALI to get and calibrate raw data
    auto dalicalib = new TArtCalibDALI();

    // Create TTree and output file for it
    TFile fout(output.c_str(), "RECREATE");
    auto tree = new TTree("tree","tree");

    // Define data nodes which are supposed to be dumped to tree 
    // EventInfo is important for the fBit information to know the trigger!
    vector<string> datanodes{"EventInfo", "BigRIPSPPAC", "BigRIPSPlastic",
                                  "BigRIPSIC","BigRIPSFocalPlane","BigRIPSRIPS",
                                  "BigRIPSTOF", "BigRIPSBeam", "DALINaI"};

    for(auto i : datanodes){
        auto *array = (TClonesArray *)sman->FindDataContainer(&i[0]);
        tree->Branch(array->GetName(), &array);
        printf("%s", array->GetName());
    }
    cout << output << endl;
    //Making new branches

    //BigRIPS
    const double magnum = -999; // Magic Number for unset events
    const int planesperppac = 4;
    vector<vector<double>> fF8PPAC(3, vector<double>(planesperppac, magnum));
    vector<vector<double>> fF8POS(3, vector<double>(2,magnum));

    Double_t fXTar = magnum, fYTar = magnum;
 
    // 6 focal planes (3,5,7,8,9,11) and 4 value (X,A,Y,B)
    vector<vector<double>> focal{6, vector<double>(4)}; // F3X, F3A

    vector<int> fGoodPPACFocus(12, 0);
    vector<int> fGoodPPACFocusOr(12, 0);

    tree->Branch("xtar",&fXTar,"fXTar/D");
    tree->Branch("ytar",&fYTar,"fYTar/D");
    /*tree->Branch("f8ppacx",fF8PPAC.at(0),"fF8PPAC[0][6]/D");
    tree->Branch("f8ppacy",fF8PPAC.at(1),"fF8PPAC[1][6]/D");
    tree->Branch("f8posx",fF8Pos.at(0),"fF8Pos[0][3]/D");
    tree->Branch("f8posy",fF8Pos.at(1),"fF8Pos[1][3]/D");*/

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

    int neve = 0;

    const vector<string> ppacname{
        "F8PPAC-1A", "F8PPAC-1B", "F8PPAC-2A", "F8PPAC-2B",
        "F9PPAC-1A", "F9PPAC-1B", "F9PPAC-2A", "F9PPAC-2B", 
        "F11PPAC-1A","F11PPAC-1B","F11PPAC-2A","F11PPAC-2B"};

    while(estore.GetNextEvent() && (neve < 2000000)){ //&& neve < 100000
        if(!(neve%10000)) printf("Event %i\n", neve);

        //Making the BigRIPS tree calibration
        brcalib->ClearData();
        brcalib->ReconstructData();

        //Reconstructiong the PID
        recopid.ClearData();
        recopid.ReconstructData();

        //Reconstructing the scattering angle:
        //Beam spot on target reconstruction
        fXTar = fYTar = magnum;

        for(auto &elem: fF8PPAC) fill(elem.begin(), elem.end(), magnum);
        for(auto &elem: fF8POS)  fill(elem.begin(), elem.end(), magnum);
        fill(fGoodPPACFocus.begin(), fGoodPPACFocus.end(), 0);
        fill(fGoodPPACFocusOr.begin(), fGoodPPACFocusOr.end(), 0);

        // Reading out the PPAC'S
        bool  posRec = true; // First PPAC has 4 fired planes

        vector<TArtPPAC*> fppac; // Storage Vector for all PPAC's
        vector<bool> fired;

        for(auto i: ppacname){ 
            fppac.push_back(ppaccalib->FindPPAC(&i[0])); // Get all PPAC planes
            fired.push_back(fppac.back() && fppac.back()->IsFiredX() &&
                            fppac.back()->IsFiredY());
        }

        // Determine well focussed events
        const vector<uint> ppacuse{8,9,11}; // Which PPAC's to use (match w/ ppacname!)

        for(uint i=0; i<ppacuse.size(); i++){ // Loop over every PPAC
            uint ppacno = ppacuse.at(i);
            int j = planesperppac;
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

        for(uint i=0; i<planesperppac; i++){ // Check first PPAC for hits
            if(fired.at(i)){
                fF8PPAC.at(0).at(i) = fppac.at(i)->GetX();
                fF8PPAC.at(1).at(i) = fppac.at(i)->GetY();
            }
            posRec *= fired.at(i); // check for hits on all planes
        }

        if(posRec){
            double mX = (fF8POS.at(0).at(1)-fF8POS.at(0).at(0))/  // dX/dZ
                        (fF8POS.at(2).at(1)-fF8POS.at(2).at(0));
            double mY = (fF8POS.at(1).at(1)-fF8POS.at(1).at(0))/  // dY/dZ
                        (fF8POS.at(2).at(1)-fF8POS.at(2).at(0));

            fXTar = fF8POS.at(0).at(1) + 880 * mX;
            fYTar = fF8POS.at(1).at(1) + 880 * mY;

            for(uint i=0; i<focalplanes.size(); i++){ // Loop all focal planes
                if(cfpl->FindFocalPlane(focalplanes.at(i))){
                    vec = cfpl->FindFocalPlane(focalplanes.at(i))->GetOptVector();
                    for(uint j=0; j<planesperppac; j++){
                        focal.at(i).at(j) = (*vec)(j);
                    }
                }
            }
        }

        //Double_t modif = - 0.0011*(F11X-8) - 0.0000003*(F11X-8)*(F11X-8)*(F11X-8)
        //                 - 0.000023*(F9X-10) - 0.0000005*(F9X)*(F9X)
        //                 + 0.000000005*(F9X+10)*(F9X+10)*(F9X+10);

        //AoQ_Z_BR->Fill(beam_br_57->GetAoQ(),beam_br_57->GetZet());
        //AoQ_Z_ZD->Fill(beam_zd_911->GetAoQ()+modif,beam_zd_911->GetZet());

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Making MINOS
        /* minoscalib->ClearData();
        minosanalyzed->ClearData();
        minostrack->ClearData();
        minosvertex->ClearData();
        minoscalib->ReconstructData();
        minosanalyzed->ReconstructData(); 
        minostrack->ReconstructData();
        minosvertex->ReconstructVertex();
        // Convert the reconstructed vertex in z (mm)
        z_vertex = (minosvertex->GetZv()*TimeBinElec - DelayTrig)*1e-3*VDrift*10; 
        //time bin*(ns/us)*vdrift(cm/us)*(mm/cm)
        */

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //Making DALI
        dalicalib->ClearData();
        //dalicalib->SetPlTime(plasticcalib->FindPlastic("F8pl")->GetTime());
        //dalicalib->SetVertex(z_vertex);   // WE DID IT... Your turn to make it 
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
    }
    cout << "Writing the tree..."<<endl;

    fout.Write();
    fout.Close();

    printf("Writing TTree complete!\n");
}