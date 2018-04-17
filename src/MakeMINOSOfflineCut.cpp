//----------! Developed by C. Santamaria / CEA Saclay !----------
//----------!      Version date :: 2014/11/12         !----------
// Code to decode the ridf format file for MINOS
// Calibration, Analysis and Tracking

// > make
// > ./MakeMINOSOfflineCut ../ridf/FileName.ridf              -> Creates a root 
//file: ../rootfiles/FileName.root
// > ./MakeMINOSOfflineCut ../ridf/FileName.ridf test.root
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "TSystem.h"
#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtDALIParameters.hh"
#include "TArtMINOSParameters.hh"
#include "TArtCalibPID.hh"
#include "TArtCalibDALI.hh"
#include "TArtCalibMINOS.hh"
#include "TArtCalibMINOSData.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtFocalPlane.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtPPAC.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtPPAC.hh"
#include "TArtBeam.hh"
#include "TArtTOF.hh"
#include "TArtRIPS.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include <vector>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph.h"
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <time.h>
#include "/Users/npaul/rootbuild/include/Math/Vector3D.h"
#include "/Users/npaul/MINOS_analysis/liboffline/TMinosClust.h"
#include "/Users/npaul/MINOS_analysis/liboffline/TMinosResult.h"
#include "/Users/npaul/MINOS_analysis/liboffline/Tracking.h"

#include "TCutG.h"

using namespace std;
using namespace ROOT::Math;

//char *ROOTFILEDIR = "./";
char *ROOTFILEDIR = "../rootfiles/";

//Definitions to use TMinuit
Tracking *Tracking_functions;
TClonesArray data_result;
TMinosResult *minosdata_result;

void SumDistance1(int &, double *, double & sum, double * par,  int);
void SumDistance2(int &, double *, double & sum, double * par, int);

// function to exit loop at keyboard interrupt.
bool stoploop = false;
void stop_interrupt(){
    printf("keyboard interrupt\n");
    stoploop = true;
}

int main(int argc, char** argv) {
    
    // Variables to calculate the elapsed time of the process
    time_t start,stop;
    time(&start);
    
    // Filename
    char* ridffile;
    char* rootfile;
    if(argc < 2)
    {
        cerr << "Missing RIDF file argument" << endl;
    }
    ridffile = argv[1];
    
    cout << " *** RIDF file: " << ridffile << endl;
    
    ofstream out;
	out.open("../counts/F7dsTotalCounts_xsruns.txt",ofstream::app);
    
    TArtStoreManager *sman = TArtStoreManager::Instance();
    
    TArtEventStore *estore = new TArtEventStore();
    estore->SetInterrupt(&stoploop);
    estore->Open(ridffile);
    TArtRawEventObject *rawevent = estore->GetRawEventObject();
    
    TArtCalibCoin*myInfo = new TArtCalibCoin();  //test AC
    myInfo->LoadData(); 
    
    
    // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes 
parameters from ".xml" files
    
//------------------------------------------------------------------------------
--------------------
    TArtBigRIPSParameters *para = TArtBigRIPSParameters::Instance();
    para->LoadParameter("db/PPAC/BigRIPSPPAC.xml.89As");
    para->LoadParameter("db/Plastic/BigRIPSPlastic.xml.89As");
    para->LoadParameter("db/IC/BigRIPSIC.xml.89As");
    para->LoadParameter("db/FocalPlane.xml");
    para->SetFocusPosOffset(8,138.5);
    
    // Create CalibPID to get and calibrate raw data ( CalibPID ->
    //[CalibPPAC , CalibIC, CalibPlastic , CalibFocalPlane]
    TArtCalibPID *brcalib     = new TArtCalibPID();
    TArtCalibPPAC *ppaccalib  = brcalib->GetCalibPPAC();
    TArtCalibPlastic *plasticcalib = brcalib->GetCalibPlastic();
    
    TArtCalibPID *cpid  = new TArtCalibPID();
    TArtCalibFocalPlane *cfpl= cpid->GetCalibFocalPlane();
    
    
    
    // Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z, ... 
(RecoPID ->
    //[ RecoTOF , RecoRIPS , RecoBeam] )
    TArtRecoPID *recopid = new TArtRecoPID();
    
    //para->PrintListOfPPACPara();
    //return;
    
    // Definition of observables we want to reconstruct
    

     ////////////////////////////////////////////// ONLINE !!!!! 
//////////////////////////////////////////////

    //Cr66 setting
    std::cout << "Defining bigrips parameters" << std::endl;
    TArtRIPS *rips3to5 = 
recopid->DefineNewRIPS(3,5,"../run/matrix/mat1.mat","D3"); // F3 - F5
    TArtRIPS *rips5to7 = 
recopid->DefineNewRIPS(5,7,"../run/matrix/mat2.mat","D5"); // F5 - F7
    TArtRIPS *rips8to9 = 
recopid->DefineNewRIPS(8,9,"../run/matrix/F8F9_LargeAccAchr.mat","D7"); // F8 - 
F9
    TArtRIPS *rips9to11 = 
recopid->DefineNewRIPS(9,11,"../run/matrix/F9F11_LargeAccAchr.mat","D8"); // F9 
- F11
    // Reconstruction of TOF DefineNewTOF(fisrt plane, second plane, time 
offset)
    /*
    // TOF offset for setting 1: 84Zn
    TArtTOF *tof3to7  = recopid->DefineNewTOF("F3pl","F7pl",304.7,5); // F3 - F7
    TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",-162.9,9); // F8 
- F11
    */
	/*   
    // TOF offset for setting 2: 110Zr
    //TArtTOF *tof3to7  = recopid->DefineNewTOF("F3pl","F7pl",304.3,5); // F3 - 
F7, transmission runs
    TArtTOF *tof3to7  = recopid->DefineNewTOF("F3pl","F7pl",304.0,5); // F3 - F7
    //TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",-162.4,9); // 
F8 - F11, old TOF offset from online anlysis
    //TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",-161.65,9); // 
F8 - F11 full target transmission run ***110Zr setting
    TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",-161.89,9); // F8 
- F11 physics run **110Zr setting
    //TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",-163.95,9); // 
F8 - F11 for empty target run **110Zr setting
    */
    // TOF offset for setting 4: 89As->88Se run
    TArtTOF *tof3to7  = recopid->DefineNewTOF("F3pl","F7pl",304.1,5); // F3 - F7
    TArtTOF *tof8to11 = recopid->DefineNewTOF("F8pl","F11pl-1",-162.35,9); // F8 
- F11 physics run 
    
    // Reconstruction of IC observables for ID
    TArtBeam *beam_br_35 = recopid->DefineNewBeam(rips3to5,tof3to7,"F7IC");
    TArtBeam *beam_br_57 = recopid->DefineNewBeam(rips5to7,tof3to7,"F7IC");
    TArtBeam *beam_br_37 = 
recopid->DefineNewBeam(rips3to5,rips5to7,tof3to7,"F7IC");
    TArtBeam *beam_zd_89 = recopid->DefineNewBeam(rips8to9,tof8to11,"F11IC");
    TArtBeam *beam_zd_911 = recopid->DefineNewBeam(rips9to11,tof8to11,"F11IC");
    TArtBeam *beam_zd_811 = 
recopid->DefineNewBeam(rips8to9,rips9to11,tof8to11,"F11IC");
     
    //to get trigger pattern
    TArtEventInfo *evtinfo=new TArtEventInfo();
    
    //  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //The Focalplane:
    Double_t F3X=-9999; Double_t F3A=-9999; Double_t F3Y=-9999; Double_t 
F3B=-9999;
    Double_t F5X=-9999; Double_t F5A=-9999; Double_t F5Y=-9999; Double_t 
F5B=-9999;
    Double_t F7X=-9999; Double_t F7A=-9999; Double_t F7Y=-9999; Double_t 
F7B=-9999;
    Double_t F8X=-9999; Double_t F8A=-9999; Double_t F8Y=-9999; Double_t 
F8B=-9999;
    Double_t F9X=-9999; Double_t F9A=-9999; Double_t F9Y=-9999; Double_t 
F9B=-9999;
    Double_t F11X=-9999; Double_t F11A=-9999; Double_t F11Y=-9999; Double_t 
F11B=-9999;
    
     //The plastic Q values:

  Double_t F3PLA_QL_raw=-9999;  Double_t F3PLA_QR_raw=-9999;  
  Double_t F7PLA_QL_raw=-9999;  Double_t F7PLA_QR_raw=-9999;  
  Double_t F8PLA_QL_raw=-9999;  Double_t F8PLA_QR_raw=-9999;  
  Double_t F11PLA_QL_raw=-9999; Double_t F11PLA_QR_raw=-9999; 

  Double_t F3PLA_Q_ave=-9999; 
  Double_t F7PLA_Q_ave=-9999;   
  Double_t F8PLA_Q_ave=-9999; 
  Double_t F11PLA_Q_ave=-9999;
  
   //The plastic time:

  Double_t F3PLA_TL_raw=-9999;  Double_t F3PLA_TR_raw=-9999;  
  Double_t F7PLA_TL_raw=-9999;  Double_t F7PLA_TR_raw=-9999;  
  Double_t F8PLA_TL_raw=-9999;  Double_t F8PLA_TR_raw=-9999;  
  Double_t F11PLA_TL_raw=-9999; Double_t F11PLA_TR_raw=-9999; 

  Double_t F3PLA_TL=-9999;  Double_t F3PLA_TR=-9999;  Double_t F3PLA_T=-9999;  
  Double_t F7PLA_TL=-9999;  Double_t F7PLA_TR=-9999;  Double_t F7PLA_T=-9999;  
  Double_t F8PLA_TL=-9999;  Double_t F8PLA_TR=-9999;  Double_t F8PLA_T=-9999;  
  Double_t F11PLA_TL=-9999; Double_t F11PLA_TR=-9999; Double_t F11PLA_T=-9999; 
  

        
  
    // Create DALIParameters to get ".xml"
    //------------------------------------
    
    TArtDALIParameters *dpara = TArtDALIParameters::Instance();
    dpara->LoadParameter("db/DALI.xml");
    
    // Create CalibDALI to get and calibrate raw data
    //-----------------------------------------------
    TArtCalibDALI *dalicalib = new TArtCalibDALI();
    // Create MINOSParameters to get ".xml"
    //------------------------------------
    TArtMINOSParameters *setup = new 
TArtMINOSParameters("MINOSParameters","MINOSParameters");
    setup->LoadParameters("./db/MINOS.xml");
    //setup->PrintListOfMINOSPara();
    
    TArtCalibMINOS *CalibMINOS = new TArtCalibMINOS();
    //Load function library
    Tracking_functions = new Tracking();
    
    
    char* infile;
    if(argv[2]==NULL) {
        char*  pch;
        char* pch2;
        pch = strtok(ridffile, "/");
        while (pch != NULL)
        {
            infile = pch;
            pch = strtok (NULL, "/");
        }
        cerr << infile << endl;
        pch2 = strtok(infile, ".");
        
        char OutFile[200]="";
        strcat(OutFile, ROOTFILEDIR);
        strcat(OutFile, pch2);
	    //strcat(OutFile, ".root");
        
        rootfile=OutFile;
    }
    else rootfile=argv[2];
    
    cout << endl;
    cout << " *** ROOT file: " << rootfile << endl;
    
    TFile *fout = new 
TFile(Form("%s_xs_noPIDcut_nocleaningcuts.root",rootfile),"RECREATE");
    TTree * tree = new TTree("tree","ridf tree");
    //tree->Branch("rawdata",&rawevent);
    
    //Here lets define the histograms that we want to fill

    TH2F*pid_37=new TH2F("pid_37","pid_37",500,2.55,2.77,500,39,45);
    TH2F*pid_811=new TH2F("pid_811","pid_811",500,2.55,2.8,500,38,45);
    TH2F*brho_cutcheck=new 
TH2F("brho_cutcheck","brhocutcheck",1000,0.92,1.1,1000,4.5,5.1);
    TH2F*f3plastic=new TH2F("f3plastic","f3plastic",100,-3,-1,100,-0.5,0.5);
    TH2F*f7plastic=new TH2F("f7plastic","f7plastic",100,-3.5,0,100,-1.5,1);
    TH2F*f8plastic=new TH2F("f8plastic","f8plastic",100,-3,3.5,100,-0.8,0.5);
    TH2F*f11plastic=new TH2F("f11plastic","f11plastic",100,-5,6,100,-3,4);
    
    //TH2F*pid_in=new TH2F("pid_in","pid_in",500,2.5,3,500,31,41);
    //TH2F*pid_out=new TH2F("pid_out","pid_out",500,2.5,3,500,31,41);
    
    
    TClonesArray fitdata;
    fitdata.SetClass("TMinosClust");
    tree->Branch("fitdata",&fitdata);
    int trackNbr;
    //tree->Branch("trackNbr",&trackNbr,"trackNbr/I");
    int trackNbr_FINAL;
    tree->Branch("trackNbr_FINAL",&trackNbr_FINAL,"trackNbr_FINAL/I");
    int evtOrig;
    tree->Branch("evtOrig",&evtOrig,"evtOrig/I");
    int padsleft;
    //tree->Branch("padsleft",&padsleft,"padsleft/I");
    data_result.SetClass("TMinosResult");
    tree->Branch("data_result",&data_result);
    double z_vertex=0., x_vertex=0., y_vertex=0., r_vertex=0., phi_vertex = 0.;
    double aoq_br=0., z_br=0., aoq_zd=0., z_zd=0.;
    // Variables only filled when trackNbr_FINAL==2
    double VDrift; //in cm/µs
    double DelayTrig; // in ns
    double StopT;
    int EventInfo_FBIT;
    double delta[4];
    double brho[6];
    tree->Branch("EventInfo_FBIT",&EventInfo_FBIT,"EventInfo_FBIT/I");
    tree->Branch("VDrift",&VDrift,"VDrift/D");
    tree->Branch("brho",brho,"brho[6]/D");    
    tree->Branch("DelayTrig",&DelayTrig,"DelayTrig/D");
    tree->Branch("aoq_br",&aoq_br,"aoq_br/D");
    tree->Branch("z_br",&z_br,"z_br/D");
    tree->Branch("aoq_zd",&aoq_zd,"aoq_zd/D");
    tree->Branch("z_zd",&z_zd,"z_zd/D");
    tree->Branch("x_vertex",&x_vertex,"x_vertex/D");
    tree->Branch("y_vertex",&y_vertex,"y_vertex/D");
    tree->Branch("z_vertex",&z_vertex,"z_vertex/D");
    tree->Branch("r_vertex",&r_vertex,"r_vertex/D");
    tree->Branch("phi_vertex",&phi_vertex,"phi_vertex/D"); // angle between two 
tracks in 3D in degrees
    double thetaz1=0., thetaz2=0.;
    tree->Branch("thetaz1",&thetaz1,"thetaz1/D"); // angle between 1st track and 
z axis in 3D in degrees
    tree->Branch("thetaz2",&thetaz2,"thetaz2/D"); // angle between 2nd track and 
z axis in 3D in degrees
    double parFit_1[4], parFit_2[4], parFit_2_F8[4], parFit_2_DSSSD[4], 
parFit_1r[4], parFit_2r[4];
    tree->Branch("parFit_1",parFit_1,"parFit_1[4]/D");
    tree->Branch("parFit_2",parFit_2,"parFit_2[4]/D");
    tree->Branch("parFit_1r",parFit_1r,"parFit_1r[4]/D");
    tree->Branch("parFit_2r",parFit_2r,"parFit_2r[4]/D");
    double errFit_1[4], errFit_2[4];
    tree->Branch("errFit_1",errFit_1,"errFit_1[4]/D");
    tree->Branch("errFit_2",errFit_2,"errFit_2[4]/D");
    tree->Branch("parFit_2_F8",parFit_2_F8,"parFit_2_F8[4]/D");
    tree->Branch("parFit_2_DSSSD",parFit_2_DSSSD,"parFit_2_DSSSD[4]/D");
    double x_DSSSD, y_DSSSD, z_DSSSD;
    tree->Branch("x_DSSSD",&x_DSSSD,"x_DSSSD/D");
    tree->Branch("y_DSSSD",&y_DSSSD,"y_DSSSD/D");
    tree->Branch("z_DSSSD",&z_DSSSD,"z_DSSSD/D");
    double x_F8PPAC, y_F8PPAC; //position @ z_DSSSD
    tree->Branch("x_F8PPAC",&x_F8PPAC,"x_F8PPAC/D");
    tree->Branch("y_F8PPAC",&y_F8PPAC,"y_F8PPAC/D");
    double z_vertexPPAC=0., x_vertexPPAC=0., y_vertexPPAC=0.;
    tree->Branch("x_vertexPPAC",&x_vertexPPAC,"x_vertexPPAC/D");
    tree->Branch("y_vertexPPAC",&y_vertexPPAC,"y_vertexPPAC/D");
    tree->Branch("z_vertexPPAC",&z_vertexPPAC,"z_vertexPPAC/D");
    double z_vertexDSSSD=0., x_vertexDSSSD=0., y_vertexDSSSD=0.;
    tree->Branch("x_vertexDSSSD",&x_vertexDSSSD,"x_vertexDSSSD/D");
    tree->Branch("y_vertexDSSSD",&y_vertexDSSSD,"y_vertexDSSSD/D");
    tree->Branch("z_vertexDSSSD",&z_vertexDSSSD,"z_vertexDSSSD/D");
    tree->Branch("delta",delta,"delta[4]/D");    
    tree->Branch("F3X",&F3X,"F3X/D");
    tree->Branch("F3A",&F3A,"F3A/D");
    tree->Branch("F3Y",&F3Y,"F3Y/D");
    tree->Branch("F3B",&F3B,"F3B/D");
    
    tree->Branch("F5X",&F5X,"F5X/D");
    tree->Branch("F5A",&F5A,"F5A/D");
    tree->Branch("F5Y",&F5Y,"F5Y/D");
    tree->Branch("F5B",&F5B,"F5B/D");
    
    tree->Branch("F7X",&F7X,"F7X/D");
    tree->Branch("F7A",&F7A,"F7A/D");
    tree->Branch("F7Y",&F7Y,"F7Y/D");
    tree->Branch("F7B",&F7B,"F7B/D");
    
    tree->Branch("F8X",&F8X,"F8X/D");
    tree->Branch("F8A",&F8A,"F8A/D");
    tree->Branch("F8Y",&F8Y,"F8Y/D");
    tree->Branch("F8B",&F8B,"F8B/D");
    
    tree->Branch("F9X",&F9X,"F9X/D");
    tree->Branch("F9A",&F9A,"F9A/D");
    tree->Branch("F9Y",&F9Y,"F9Y/D");
    tree->Branch("F9B",&F9B,"F9B/D");
    
    tree->Branch("F11X",&F11X,"F11X/D");
    tree->Branch("F11A",&F11A,"F11A/D");
    tree->Branch("F11Y",&F11Y,"F11Y/D");
    tree->Branch("F11B",&F11B,"F11B/D");
    
    tree->Branch("F8PLA_Q_ave",&F8PLA_Q_ave,"F8PLA_Q_ave/D");
    tree->Branch("F8PLA_T",&F8PLA_T,"F8PLA_T");
    tree->Branch("F11PLA_T",&F8PLA_T,"F11PLA_T");
    
    
   //EventInfo is important for the fBit information to know the trigger!
    TClonesArray *info_array = (TClonesArray 
*)sman->FindDataContainer("EventInfo");
    std::cout<<info_array->GetName()<<std::endl;
    Int_t EventInfo_fBit; //=(Int_t)((TArtEventInfo 
*)info_array->At(0))->GetTriggerBit(); //test AC
    tree->Branch(info_array->GetName(),&info_array);
   
   
    /*
    // define data nodes which are supposed to be dumped to tree
    //EventInfo is important for the fBit information to know the trigger!
    TClonesArray * info_array = (TClonesArray 
*)sman->FindDataContainer("EventInfo");
    std::cout<<info_array->GetName()<<std::endl;
    Int_t EventInfo_fBit = (Int_t)((TArtEventInfo 
*)info_array->At(0))->GetTriggerBit();
    tree->Branch(info_array->GetName(),&info_array);
    */
    
    
    
    //Dali data
    TClonesArray * dali_array=
    (TClonesArray *)sman->FindDataContainer("DALINaI");
    tree->Branch(dali_array->GetName(),&dali_array);
    
    //PID reconstructed data:
    TClonesArray *beam_array =
    (TClonesArray *)sman->FindDataContainer("BigRIPSBeam");
    std::cout<<beam_array->GetName()<<std::endl;
    tree->Branch(beam_array->GetName(),&beam_array);
    
    
    //PPAC array
    TClonesArray * ppac_array =
    (TClonesArray *)sman->FindDataContainer("BigRIPSPPAC");
    std::cout<<ppac_array->GetName()<<std::endl;
    tree->Branch(ppac_array->GetName(),&ppac_array);
    
    //Plastic array
    TClonesArray * pla_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSPlastic");
    tree->Branch(pla_array->GetName(),&pla_array);
    
    //ionization chamber array
    TClonesArray * ic_array = 
    (TClonesArray *)sman->FindDataContainer("BigRIPSIC");
    std::cout<<ic_array->GetName()<<std::endl;
    tree->Branch(ic_array->GetName(),&ic_array);
    
  // define data nodes which are supposed to be dumped to tree 

     
    //para->PrintListOfPPACPara();

    //%%%%%%%%%%%%%%%%%%%%%%
    //DALI
    Int_t dalimultwotime = 0;
    Int_t dalimult = 0;
    Int_t dalitimetruemult = 0;
    Int_t dalimultthres = 0;
    Int_t dalitimetruemultthres = 0;
    
    tree->Branch("dalimultwotime",&dalimultwotime,"dalimultwotime/I");
    tree->Branch("dalimult",&dalimult,"dalimult/I");
    tree->Branch("dalitimetruemult",&dalitimetruemult,"dalitimetruemult/I");
    tree->Branch("dalimultthres",&dalimultthres,"dalimultthres/I");
    
tree->Branch("dalitimetruemultthres",&dalitimetruemultthres,
"dalitimetruemultthres/I");
    
    
    //--------------------------------------------------------------------
    //PID Cuts

    // For 113Tc(p,2p)112Mo
    //TFile *BRcut = new TFile("../cut/brcut_113Tc.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_112Mo.root","READ");
    
    // For 112Mo inelastic
    //TFile *BRcut = new TFile("../cut/brcut_112Mo.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_112Mo.root","READ");
    
    // For 110Zr inelastic
    //TFile *BRcut = new TFile("../cut/brcut_110Zr.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_110Zr.root","READ");

    // For 110Zr
    TFile *BRcut = new TFile("../cut/brcut_Z42_A112.root","READ");
    TFile *ZDcut = new TFile("../cut/zdcut_Z40_A110.root","READ");
    
    // For 108Zr inelastic
    //TFile *BRcut = new TFile("../cut/brcut_108Zr.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_108Zr.root","READ");
    
    // For 108Zr
    //TFile *BRcut = new TFile("../cut/brcut_109nb.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_108Zr.root","READ");
    
    // For 112Mo(p,3p)110Zr
    //TFile *BRcut = new TFile("../cut/brcut_112Mo.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_110zr.root","READ");
    
    // For 113Mo(p,3pn)110Zr
    //TFile *BRcut = new TFile("../cut/brcut_113mo.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_110zr.root","READ");

	// For 112Nb(p,2pn)110Zr
    //TFile *BRcut = new TFile("../cut/brcut_112nb.root","READ");
    //TFile *ZDcut = new TFile("../cut/zdcut_110zr.root","READ");

    // For 88Ge
    // TFile *BRcut = new TFile("../cut/88Ge/brcut_89As.root","READ");
    // TFile *ZDcut = new TFile("../cut/88Ge/zdcut_88Ge.root","READ");
    // For 84Se
    //TFile *BRcut = new TFile("../cut/94Se/brcut_95Br.root","READ");
    //TFile *ZDcut = new TFile("../cut/94Se/zdcut_94Se.root","READ");

    // For 96Kr
    //TFile *BRcut = new TFile("../cut/96Kr/brcut_97Rb.root","READ");
    //TFile *ZDcut = new TFile("../cut/96Kr/zdcut_96Kr.root","READ");

    // For 90Se
    //TFile *BRcut = new TFile("../cut/90Se/brcut_91Br.root","READ");
    //TFile *ZDcut = new TFile("../cut/90Se/zdcut_90Se.root","READ");

    // For 98Kr
    //TFile *BRcut = new TFile("../cut/98Kr/brcut_99Rb.root","READ");
    //TFile *ZDcut = new TFile("../cut/98Kr/zdcut_98Kr.root","READ");

    // For 100Kr
    //TFile *BRcut = new TFile("../cut/100Kr/brcut_101Rb.root","READ");
    //TFile *ZDcut = new TFile("../cut/100Kr/zdcut_100Kr.root","READ");

    // For 100Sr (for checki from (p,pn) )
    //TFile *BRcut = new TFile("../cut/100Sr/brcut_101Sr.root","READ");
    //TFile *ZDcut = new TFile("../cut/100Sr/zdcut_100Sr.root","READ");
    
    //OTHER CUTS TO CLEAN UP SPECTRUM
    
     TFile *Brhocut = new TFile("../cut/brhozdcut_setting3.root","READ");
     TFile *F3plasticcut = new 
TFile("../cut/f3plasticcut_setting3.root","READ");
     TFile *F7plasticcut = new 
TFile("../cut/f7plasticcut_setting3.root","READ");
     TFile *F8plasticcut = new 
TFile("../cut/f8plasticcut_setting3.root","READ");
     TFile *F11plasticcut = new 
TFile("../cut/f11plasticcut_setting3.root","READ");


    //TCutG *brcut1;
    //TCutG *brcut2;
    //TCutG *brcut3;
    TCutG *brcut;
    TCutG *zdcut;
    TCutG *brhocut;
    TCutG *f3plasticcut;
    TCutG *f7plasticcut;
    TCutG *f8plasticcut;
    TCutG *f11plasticcut;
    //BRcut1->GetObject("CUTG",brcut1);
    //BRcut2->GetObject("CUTG",brcut2);
    //BRcut3->GetObject("CUTG",brcut3);
    BRcut->GetObject("CUTG",brcut);
    ZDcut->GetObject("CUTG",zdcut);
    Brhocut->GetObject("CUTG",brhocut);
    F3plasticcut->GetObject("CUTG",f3plasticcut);
    F7plasticcut->GetObject("CUTG",f7plasticcut);
    F8plasticcut->GetObject("CUTG",f8plasticcut);
    F11plasticcut->GetObject("CUTG",f11plasticcut);
    
    
    // Parameters for the MINOS ANALYSIS
    double MINOSthresh;
    double TimeBinElec; //in ns
    double Tshaping; // in ns
    
    ifstream ConfigFile;
    ConfigFile.open("../configs/ConfigMINOSDrift_fixed.txt");
    string HeaderFile;
    string dummystring;
    bool ThisFileConfig = false;
    getline(ConfigFile,dummystring);
    ConfigFile >> MINOSthresh >> TimeBinElec >> Tshaping;
    cout << MINOSthresh << " " << TimeBinElec << " " << Tshaping << endl;
    getline(ConfigFile,dummystring, '\n'); //takes present line
    getline(ConfigFile,dummystring, '\n'); //discard comment line
    while(ConfigFile.is_open())
    {
        ConfigFile >> HeaderFile >> DelayTrig >> StopT >> VDrift;
        if(HeaderFile == infile) {
            cout << HeaderFile << endl;
            break;
        }
    }
    
    ConfigFile.close();
    
    cout << endl;
    cout << " *** MINOS Configuration Parameters *** " << endl;
    cout << " Electronics        :::   MINOSthresh = " << MINOSthresh << " 
(bins)  ;  TimeBinElec = " << TimeBinElec << " (ns)  ;    Tshaping = " << 
Tshaping << " (ns)  ;  DelayTrig = " << DelayTrig << " (ns)  ;   VDrift = " << 
VDrift << " (cm/micros)" << endl;
    cout << endl;
    
    double PI = TMath::Pi();
    int InTrig = 0;
    int neve = 0;
    
    //    vector<TCanvas*> Filter_canvas;
    double ChargeBin=0.,maxCharge=0.;
    int filled = 0;
    vector<double> Xpad, Ypad, Qpad, XpadNew, YpadNew, QpadNew, ZpadNew;
    vector<int> clusterringbool;
    vector<int> clusternbr;
    vector<int> clusterpads;
    int Iteration=0;
    int filter_result=0;
    //int padsleft=0;
    int indexfill=0;
    bool fitbool = false;
    int fit2DStatus = 0;
    double Chi2=0.;
    double x_mm,y_mm,z_mm,q_pad,t_pad;
    
    //DSSSD Variables
    vector<double> XpadDSSSD, YpadDSSSD, QXpadDSSSD, QYpadDSSSD;
    double XbaryDSSSD, YbaryDSSSD;
    double ThDSSSD = 1000.;

    
    TF1 *fit_function = new 
TF1("fit_function",Tracking_functions,&Tracking::conv_fit, 0, 511, 
3,"Tracking","conv_fit"); // only 3 param. because baseline is fixed
    double hfit_max, hfit_max_T, T_min, T_max;
    TH1F *hfit = new TH1F("hfit","hfit",512,0,512);
    
    
    //2nd step variables
    int npoint_temp=0, cluster_temp=0;
    int padsleft2=0;
    int cluster1=0, cluster2=0;
    int ringsum=0;
    int ringtouch[18]={0};
    double zmax=0.;
    int allevt_2pfiltered=0, allevt_2pvertex=0;
    int array_final=0;
    //    vector<TCanvas*> Hough_canvas;
    bool padsnbr1=false, padsnbr2=false;
    double Qmean1=0., Qmean2=0.;
    TGraph * gryz_1;
    TGraph * grxz_1;
    TGraph * gryz_2;
    TGraph * grxz_2;
    vector<double> xin, yin, zin, qin, xout, yout, zout, qout;
    vector<int> cluster;
    
    
    
    TMinuit *min ;
    //while(estore->GetNextEvent() && neve<100000){
    while(estore->GetNextEvent() ){
        if(neve%1000==0) cout << "Event " << neve << endl;
        //Clear & Reset variables
        fitdata.Clear();
        data_result.Clear();
        Xpad.clear();
        Ypad.clear();
        Qpad.clear();
        XpadNew.clear();
        YpadNew.clear();
        QpadNew.clear();
        clusterringbool.clear();
        clusternbr.clear();
        clusterpads.clear();
        //	Filter_canvas.clear();
        hfit->Reset();
        XpadDSSSD.clear();
        YpadDSSSD.clear();
        QXpadDSSSD.clear();
        QYpadDSSSD.clear();
        
        
        filled=0;
        indexfill = 0;
        ChargeBin = 0.;
        maxCharge = 0.;
        Iteration=0;
        filter_result=0;
        fit2DStatus=0;
        EventInfo_FBIT=0;
        trackNbr=0;
        trackNbr_FINAL=0;
        padsleft=0;
        x_mm = 0.; y_mm = 0.; z_mm = 0.; q_pad = 0.; t_pad = 0.;
        fitbool = false;
        array_final=0;
        ringsum=0;
        cluster1=0; cluster2=0;
        z_vertex=-10000.; x_vertex=-10000.; y_vertex=-10000.; r_vertex=-10000.;
        z_vertexPPAC=-10000.; x_vertexPPAC=-10000.; y_vertexPPAC=-10000.;
        z_vertexDSSSD=-10000.; x_vertexDSSSD=-10000.; y_vertexDSSSD=-10000.;
        zmax=0.;
        padsnbr1=false;
        padsnbr2=false;
        Qmean1=0.; Qmean2=0.;
        padsleft2=0;
        thetaz1=-10.; thetaz2=-10.;
        phi_vertex=0.;
        XbaryDSSSD = 0.; YbaryDSSSD = 0.;
        x_DSSSD = -10000.; y_DSSSD = -10000.; z_DSSSD = -10000.;
        x_F8PPAC = -10000.; y_F8PPAC = -10000.;
        
        int npoint1=0,npoint2=0;
        gryz_1 = new TGraph();
        grxz_1 = new TGraph();
        gryz_2 = new TGraph();
        grxz_2 = new TGraph();
        
        for(int iit=0;iit<4;iit++) {
            parFit_1[iit] = -1000.;
            errFit_1[iit] = -1000.;
            parFit_2[iit] = -1000.;
            errFit_2[iit] = -1000.;
            parFit_2_F8[iit] = -1000.;
            parFit_2_DSSSD[iit] = -1000.;
        }
        
        evtOrig = neve;
        
        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
        //Making the BigRIPS tree calibration and find cuts and triggers
      
        
        brcalib->ClearData();
        brcalib->ReconstructData();
        
        cpid->ClearData();
        cpid->ReconstructData();
        cfpl = cpid->GetCalibFocalPlane();
        
        //Reconstructing the PID
        recopid->ClearData();
        recopid->ReconstructData();
        
        myInfo->LoadData(); //this loads the trigger bit information 
        EventInfo_FBIT = (Int_t)((TArtEventInfo 
*)info_array->At(0))->GetTriggerBit();
        //cout<<EventInfo_FBIT<<" "<<((TArtEventInfo 
*)info_array->At(0))->GetTriggerBit()<<" "<<((TArtEventInfo 
*)info_array->At(0))->GetTimeStamp()<<endl;
        
        //The Focalplane:
        F3X=-9999; F3A=-9999; F3Y=-9999; F3B=-9999;
        F5X=-9999; F5A=-9999; F5Y=-9999; F5B=-9999;
        F7X=-9999; F7A=-9999; F7Y=-9999; F7B=-9999;
        F8X=-9999; F8A=-9999; F8Y=-9999; F8B=-9999;
        F9X=-9999; F9A=-9999; F9Y=-9999; F9B=-9999;
        F11X=-9999; F11A=-9999; F11Y=-9999; F11B=-9999;
        
        TArtFocalPlane *tfpl;
        TVectorD *vec;
        tfpl = cfpl->FindFocalPlane(3);
        if(tfpl){vec=tfpl->GetOptVector(); F3X=(*vec)(0); F3A=(*vec)(1); 
F3Y=(*vec)(2); F3B=(*vec)(3);}
        tfpl = cfpl->FindFocalPlane(5);
        if(tfpl){vec=tfpl->GetOptVector(); F5X=(*vec)(0); F5A=(*vec)(1); 
F5Y=(*vec)(2); F5B=(*vec)(3);}
        tfpl = cfpl->FindFocalPlane(7);
        if(tfpl){vec=tfpl->GetOptVector(); F7X=(*vec)(0); F7A=(*vec)(1); 
F7Y=(*vec)(2); F7B=(*vec)(3);}
        tfpl = cfpl->FindFocalPlane(8);
        if(tfpl){vec=tfpl->GetOptVector(); F8X=(*vec)(0); F8A=(*vec)(1); 
F8Y=(*vec)(2); F8B=(*vec)(3);}
        tfpl = cfpl->FindFocalPlane(9);
        if(tfpl){vec=tfpl->GetOptVector(); F9X=(*vec)(0); F9A=(*vec)(1); 
F9Y=(*vec)(2); F9B=(*vec)(3);}
        tfpl = cfpl->FindFocalPlane(11);
        if(tfpl){vec=tfpl->GetOptVector(); F11X=(*vec)(0); F11A=(*vec)(1); 
F11Y=(*vec)(2); F11B=(*vec)(3);}
        
            //The plastic Q values:

  Double_t F3PLA_QL_raw=-9999;  Double_t F3PLA_QR_raw=-9999;  
  Double_t F7PLA_QL_raw=-9999;  Double_t F7PLA_QR_raw=-9999;  
  Double_t F8PLA_QL_raw=-9999;  Double_t F8PLA_QR_raw=-9999;  
  Double_t F11PLA_QL_raw=-9999; Double_t F11PLA_QR_raw=-9999; 

  Double_t F3PLA_Q_ave=-9999; 
  Double_t F7PLA_Q_ave=-9999;   
  Double_t F8PLA_Q_ave=-9999; 
  Double_t F11PLA_Q_ave=-9999;
  
   //The plastic time:

  Double_t F3PLA_TL_raw=-9999;  Double_t F3PLA_TR_raw=-9999;  
  Double_t F7PLA_TL_raw=-9999;  Double_t F7PLA_TR_raw=-9999;  
  Double_t F8PLA_TL_raw=-9999;  Double_t F8PLA_TR_raw=-9999;  
  Double_t F11PLA_TL_raw=-9999; Double_t F11PLA_TR_raw=-9999; 

  Double_t F3PLA_TL=-9999;  Double_t F3PLA_TR=-9999;  Double_t F3PLA_T=-9999;  
  Double_t F7PLA_TL=-9999;  Double_t F7PLA_TR=-9999;  Double_t F7PLA_T=-9999;  
  Double_t F8PLA_TL=-9999;  Double_t F8PLA_TR=-9999;  Double_t F8PLA_T=-9999;  
  Double_t F11PLA_TL=-9999; Double_t F11PLA_TR=-9999; Double_t F11PLA_T=-9999; 
  
   TArtPlastic *tpla;
        tpla = plasticcalib->FindPlastic("F3pl");
        if(tpla){
            F3PLA_TL_raw = tpla->GetTLRaw(); F3PLA_TR_raw = tpla->GetTRRaw(); 
            F3PLA_TL = tpla->GetTimeL(); F3PLA_TR = tpla->GetTimeR(); 
            F3PLA_T = tpla->GetTime();
            
            //cout<<tpla->GetTLRaw()<<endl;
            //cout<<"F3PLA_TL_raw  "<<F3PLA_TL_raw<<endl;
                        
            F3PLA_QL_raw = tpla->GetQLRaw(); F3PLA_QR_raw = tpla->GetQRRaw();
            F3PLA_Q_ave = tpla->GetQAveRaw();
        }

        tpla = plasticcalib->FindPlastic("F7pl");
        if(tpla){
            F7PLA_TL_raw = tpla->GetTLRaw(); F7PLA_TR_raw = tpla->GetTRRaw(); 
            F7PLA_TL = tpla->GetTimeL(); F7PLA_TR = tpla->GetTimeR(); 
            F7PLA_T = tpla->GetTime();
      
            F7PLA_QL_raw = tpla->GetQLRaw(); F7PLA_QR_raw = tpla->GetQRRaw(); 
            F7PLA_Q_ave = tpla->GetQAveRaw();
        }
    
        tpla = plasticcalib->FindPlastic("F8pl");
        if(tpla){
            F8PLA_TL_raw = tpla->GetTLRaw(); F8PLA_TR_raw = tpla->GetTRRaw(); 
            F8PLA_TL = tpla->GetTimeL(); F8PLA_TR = tpla->GetTimeR(); 
            F8PLA_T = tpla->GetTime();

            F8PLA_QL_raw = tpla->GetQLRaw(); F8PLA_QR_raw = tpla->GetQRRaw(); 
            F8PLA_Q_ave = tpla->GetQAveRaw();
        }

        tpla = plasticcalib->FindPlastic("F11pl-1");
        if(tpla){
            F11PLA_TL_raw = tpla->GetTLRaw(); F11PLA_TR_raw = tpla->GetTRRaw(); 
            F11PLA_TL = tpla->GetTimeL(); F11PLA_TR = tpla->GetTimeR(); 
            F11PLA_T = tpla->GetTime();

            F11PLA_QL_raw = tpla->GetQLRaw(); F11PLA_QR_raw = tpla->GetQRRaw(); 
            F11PLA_Q_ave = tpla->GetQAveRaw();
            
            //cout<<"F11PLA_TL_raw  "<<F11PLA_TL_raw<<endl;
        }
        

        
        //Starting with the analysis conditions
        //bool brCutBool = true;
        //bool zdCutBool = true;
        
    
        bool brCutBool = true;
        bool zdCutBool = true;
        bool brhoBool = true;
        bool f3plasticcutBool = true;
        bool f7plasticcutBool = true;
        bool f8plasticcutBool = true;
        bool f11plasticcutBool = true;
   
        double betazd;
                
        brho[0]=beam_br_35->GetBrho(); 
        brho[1]=beam_br_57->GetBrho();
        brho[2]=beam_br_37->GetBrho(); 
        brho[3]=beam_zd_89->GetBrho();
        brho[4]=beam_zd_911->GetBrho();
        brho[5]=beam_zd_811->GetBrho();
        
           	
    	//define the PID parameters to cut on
    	betazd=beam_zd_811->GetBeta();
    	//aoq_br=beam_br_37->GetAoQ();
       	//z_br=beam_br_37->GetZet();
       

       	//for 110Zr setting
       	
//aoq_zd=beam_zd_811->GetAoQ()+0.000065*F9X+0.00000025*(F9X-20)*(F9X-20)+0.0001*
F9A-0.00003*F9A*F9A-0.0001*F11X-0.00006*F11A-0.000001*(F11A+6)*(F11A+6);
        //z_zd=beam_zd_811->GetZet()+(betazd-0.4974)*12;
        
        //for 89As setting
        z_br=beam_br_37->GetZet()-0.01+1.36932*(beam_br_37->GetBeta()-0.645);
        
aoq_br=beam_br_37->GetAoQ()+0.000027*(F7PLA_Q_ave-650)-0.00005*(F3PLA_Q_ave-440)
+0.0002*F3X+0.00001*F3A-0.000003*F3A*F3A-0.00003*F5A+0.000003*F5A*F5A+0.00002*
F5X-0.0000005*F5X*F5X+0.00006*F7X+0.00001*F7X*F7X-0.00011*F7A;
        aoq_zd = 
beam_zd_811->GetAoQ()+0.001-0.0003*(F8PLA_Q_ave-240)+0.000035*(F11PLA_Q_ave-650)
+0.0006*F8X-0.00001*F8X*F8X-0.000025*(F9X-80)-0.0000004*(F9X-80)*(F9X-80)-0.0002
*F9A+0.0003*F11X-0.00001*F11X*F11X-0.00085*(F11A+10);
        z_zd=beam_zd_811->GetZet()+0.03-3.38688*(betazd-0.55);
        
         
        if( brcut->IsInside(aoq_br,z_br) ) brCutBool = true;
       
        if ( 
f3plasticcut->IsInside((F3PLA_TR-F3PLA_TL),log(F3PLA_QL_raw/F3PLA_QR_raw))) 
f3plasticcutBool = true;
  
        if ( 
f7plasticcut->IsInside((F7PLA_TR-F7PLA_TL),log(F7PLA_QL_raw/F7PLA_QR_raw))) 
f7plasticcutBool = true;

        if ( 
f8plasticcut->IsInside((F8PLA_TR-F8PLA_TL),log(F8PLA_QL_raw/F8PLA_QR_raw))) 
f8plasticcutBool = true; 

        if ( 
f11plasticcut->IsInside((F11PLA_TR-F11PLA_TL),log(F11PLA_QL_raw/F11PLA_QR_raw))) 
f11plasticcutBool = true; 
 
        if ( brhocut->IsInside((brho[3]/brho[4]),brho[3])) brhoBool = true;      
    
       
           
       if( zdcut->IsInside(aoq_zd,z_zd) ) zdCutBool = true;
       //if( zdcut->IsInside(aoqc,zetc) ) zdCutBool = true;
        
        //Trigger register information
         if( brCutBool == true && zdCutBool == true && EventInfo_fBit !=6) 
InTrig++;
          
       cout<<"db6"<<endl; 
        delta[0] = rips3to5->GetDelta();
	    delta[1] = rips5to7->GetDelta();
	    delta[2] = rips8to9->GetDelta();
	    delta[3] = rips9to11->GetDelta();
    
      
    //cout << EventInfo_fBit << endl ;
        
        //Need to check the settings for PSP
        //if(EventInfo_fBit==13||EventInfo_fBit==15) F7F11DALITrigger = true;
        
        if(brCutBool == false || zdCutBool == false  || brhoBool == false || 
f3plasticcutBool == false || f7plasticcutBool == false || f8plasticcutBool == 
false || f11plasticcutBool == false ){
            neve++;
            continue;
        }
        
    //histograms to check and make sure that the cuts are applied correctly    
    pid_37->Fill(aoq_br,z_br); 
    pid_811->Fill(aoq_zd,z_zd); 
    brho_cutcheck->Fill((brho[3]/brho[4]),brho[3]);
    f3plastic->Fill((F3PLA_TR-F3PLA_TL),log(F3PLA_QL_raw/F3PLA_QR_raw));
    f7plastic->Fill((F7PLA_TR-F7PLA_TL),log(F7PLA_QL_raw/F7PLA_QR_raw));
    f8plastic->Fill((F8PLA_TR-F8PLA_TL),log(F8PLA_QL_raw/F8PLA_QR_raw));
    f11plastic->Fill((F11PLA_TR-F11PLA_TL),log(F11PLA_QL_raw/F11PLA_QR_raw));
    
        cout<<"db4"<<endl; 
    
    
        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
        //Making DALI
        dalicalib->ClearData();
        //dalicalib->SetPlTime(plasticcalib->FindPlastic("F8pl")->GetTime());
        //Add above to remove F8plastic tof.
        dalicalib->ReconstructData();
        
        dalimultwotime = dalicalib->GetMultWithoutT();
        dalimult = dalicalib->GetMult();
        dalitimetruemult = dalicalib->GetTimeTrueMult();
        dalimultthres = dalicalib->GetMultThres();
        dalitimetruemultthres = dalicalib->GetTimeTrueMultThres();
        
                      
         //cout<<"EventInfo_fBit "<<EventInfo_fBit<<endl;
        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
        
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
        //Making MINOS OFFLINE Reconstruction
       
        CalibMINOS->ClearData();
        CalibMINOS->ReconstructData();

        TMinosClust *minosfitdata;
        TArtCalibMINOSData *minos = new TArtCalibMINOSData();
        
        
        //1:MINOS://///////////////// Filling vectors with (x,y) information for 
TPC & DSSSD ///////////////////:MINOS://
        
        for(Int_t i=0;i<CalibMINOS->GetNumCalibMINOS();i++) {
            minos = CalibMINOS->GetCalibMINOS(i);
            ChargeBin = 0.;
            maxCharge = 0.;
            x_mm = minos->GetX();
            y_mm = minos->GetY();
            
            if(minos->GetDetID() >= 1) continue;
            
            // look at electronic signals from TPC
            else if(minos->GetDetID() == 0) {
                if( !(abs(x_mm)<0.0001 && abs(y_mm)<0.0001) ) { // NON connected 
MINOS (TPC) channels...
                    for(Int_t j=0; j<minos->GetNData(); j++) {
                        if(minos->GetCalibValue(j)>maxCharge) maxCharge = 
minos->GetCalibValue(j);
                    }
                    if(maxCharge>=MINOSthresh) /*continue;*/{
                        Xpad.push_back(x_mm);
                        Ypad.push_back(y_mm);
                        Qpad.push_back(maxCharge);
                        filled++;
                    }
                }
            }
            
            
        }
        
        
        //2:MINOS://////////////////////// Modified (xy) Hough transform (only 
for TPC) ////////////////////////:MINOS://
        if(filled>0) {
            padsleft = Xpad.size();
            
            while(Xpad.size()>=10 && Iteration<20) {
                filter_result = 0;
                Iteration++;
                //  			Filter_canvas.push_back(new 
TCanvas(Form("Event%d_cluster%d", i, Iteration), Form("Event%d_cluster%d", i, 
Iteration)));
                //			filter_result = Hough_modified(Filter_canvas.back(), 
&Xpad, &Ypad, &Qpad, &XpadNew, &YpadNew, &QpadNew, &clusterringbool);
                filter_result = Tracking_functions->Hough_modified(&Xpad, &Ypad, 
&Qpad, &XpadNew, &YpadNew, &QpadNew, &clusterringbool);
                //cout << "Obertelli filter n°" << Iteration << " ::: size of 
cluster=" << filter_result << " , " << XpadNew.size() << " , size of old pads=" 
<< Xpad.size() << endl;
                //			Filter_canvas.back()->Write();
                for(int ik=0; ik<filter_result; ik++) {
                    clusternbr.push_back(Iteration);
                    clusterpads.push_back(filter_result);
                }
                if(filter_result>10 && clusterringbool.back()==1) trackNbr++;
            }
        }
        
        for(unsigned int il=0; il<XpadNew.size(); il++)
        {
            minosfitdata = (TMinosClust*)fitdata.ConstructedAt(il);
            minosfitdata->Set(XpadNew[il], YpadNew[il],-10000, -10000, 
QpadNew[il], clusternbr[il],clusterpads[il], 0.);
            ZpadNew.push_back(-10000.);
        }
        
        
        //////////////////////////////////////// Get beam direction from F8PPAC 
////////////////////////////////////////
        TArtPPAC *tppac;
        Double_t F8PPAC1A[3], F8PPAC1B[3], F8PPAC2A[3], F8PPAC2B[3];
        Double_t F8PPAC1[3], F8PPAC2[3];
        tppac = ppaccalib->FindPPAC("F8PPAC-1A");
        if(tppac){F8PPAC1A[0] = tppac->GetX(); F8PPAC1A[1] = 
tppac->GetY();F8PPAC1A[2]=-1310.7;}
        tppac = ppaccalib->FindPPAC("F8PPAC-1B");
        if(tppac){F8PPAC1B[0] = tppac->GetX(); F8PPAC1B[1] = 
tppac->GetY();F8PPAC1B[2]=-1273.3;}
        tppac = ppaccalib->FindPPAC("F8PPAC-2A");
        if(tppac){F8PPAC2A[0] = tppac->GetX(); F8PPAC2A[1] = 
tppac->GetY();F8PPAC2A[2]=-802.1;}
        tppac = ppaccalib->FindPPAC("F8PPAC-2B");
        if(tppac){F8PPAC2B[0] = tppac->GetX(); F8PPAC2B[1] = 
tppac->GetY();F8PPAC2B[2]=-781.9;}
        double z_F8=-144; 
        if(tppac){
            F8PPAC1[2]=(F8PPAC1A[2]+F8PPAC1B[2])/2+z_F8;//this is the time
            F8PPAC2[2]=(F8PPAC2A[2]+F8PPAC2B[2])/2+z_F8;
            bool flag1=false, flag2=false;
            bool flag1A[2] = {false}, flag1B[2] = {false}, flag2A[2] = {false}, 
flag2B[2] = {false};
            for(int i1=0; i1<2; i1++) {
                if(F8PPAC1A[i1]<100 && F8PPAC1A[i1]>-100) flag1A[i1]=true; 
//condition that signals are physical
                if(F8PPAC1B[i1]<100 && F8PPAC1B[i1]>-100) flag1B[i1]=true;
                if(F8PPAC2A[i1]<100 && F8PPAC2A[i1]>-100) flag2A[i1]=true;
                if(F8PPAC2B[i1]<100 && F8PPAC2B[i1]>-100) flag2B[i1]=true;
            }
            
            if(flag1A[0]==true && flag1A[1]==true && flag1B[0]==true && 
flag1B[1]==true && flag2A[0]==false && flag2A[1]==false && flag2B[0]==false && 
flag2B[1]==false) {
                F8PPAC1[0]=F8PPAC1A[0]; //case where only first PPAC is hit
                F8PPAC1[1]=F8PPAC1A[1];
                F8PPAC2[0]=F8PPAC1B[0];
                F8PPAC2[1]=F8PPAC1B[1];
            }
            if(flag2A[0]==true && flag2A[1]==true && flag2B[0]==true && 
flag2B[1]==true && flag1A[0]==false && flag1A[1]==false && flag1B[0]==false && 
flag1B[1]==false) {
                F8PPAC1[0]=F8PPAC2A[0];  //case where only second PPAC is hit
                F8PPAC1[1]=F8PPAC2A[1];
                F8PPAC2[0]=F8PPAC2B[0];
                F8PPAC2[1]=F8PPAC2B[1];
            }
            
            
            if(flag1A[0]==true && flag1A[1]==true && flag1B[0]==true && 
flag1B[1]==true) {
                F8PPAC1[0]=(F8PPAC1A[0]+F8PPAC1B[0])/2; //take average of ppac1 
x signal
                F8PPAC1[1]=(F8PPAC1A[1]+F8PPAC1B[1])/2; //take average of ppac1 
y signal
            }
            else if(flag1A[0]==true && flag1A[1]==true) {
                F8PPAC1[0]=F8PPAC1A[0];
                F8PPAC1[1]=F8PPAC1A[1];
            }
            else if(flag1B[0]==true && flag1B[1]==true) {
                F8PPAC1[0]=F8PPAC1B[0];
                F8PPAC1[1]=F8PPAC1B[1];
            }
            else flag1=true;
            ////
            if(flag2A[0]==true && flag2A[1]==true && flag2B[0]==true && 
flag2B[1]==true) {
                F8PPAC2[0]=(F8PPAC2A[0]+F8PPAC2B[0])/2;
                F8PPAC2[1]=(F8PPAC2A[1]+F8PPAC2B[1])/2;
            }
            else if(flag2A[0]==true && flag2A[1]==true) {
                F8PPAC2[0]=F8PPAC2A[0];
                F8PPAC2[1]=F8PPAC2A[1];
            }
            else if(flag2B[0]==true && flag2B[1]==true) {
                F8PPAC2[0]=F8PPAC2B[0];
                F8PPAC2[1]=F8PPAC2B[1];
            }
            else flag2=true;
            
            if(flag1==true && flag2 == true) continue;
            
            if(flag2A[0]==false && flag2A[1]==false && flag2B[0]==false && 
flag2B[1]==false || flag1A[0]==false && flag1A[1]==false && flag1B[0]==false && 
flag1B[1]==false) continue;
            
            
Tracking_functions->ParFor_Vertex(F8PPAC1,F8PPAC2,parFit_2_F8);//Calc line pars. 
used in Vertex function
        
         
         //Adjust for experimental offset seen in vertex position
         double MINOSxoffset=-3.5; double MINOSyoffset=2.0;
         parFit_2_F8[0]= parFit_2_F8[0]+MINOSxoffset;
         parFit_2_F8[2]= parFit_2_F8[2]+MINOSyoffset;
         
        }
        
        //cout << parFit_2_F8[0] << " ; " << parFit_2_F8[1] << " ; " << 
parFit_2_F8[2] << " ; " << parFit_2_F8[3] << endl;

        
        
        ////////////////////        following analysis ONLY if #tracks = 1 || 2 
|| 3 || 4  for TPC       ///////////////////
        
        if(trackNbr>0 && trackNbr<5)
        {
            if(filled==0) cerr << "Error !!!" << endl;
            
            //3:MINOS://////////// Fitting the taken pads for Qmax and Ttrig 
information ////////////:MINOS://
            
            padsleft = padsleft - XpadNew.size();
            
            for(Int_t i=0;i<CalibMINOS->GetNumCalibMINOS();i++) {
                minos = CalibMINOS->GetCalibMINOS(i);
                hfit->Reset();
                fitbool = false;
                
                if(minos->GetDetID() != 0) continue;
                
                x_mm = minos->GetX();
                y_mm = minos->GetY();
                
                for(unsigned int jj=0; jj<XpadNew.size(); jj++) {
                    if( abs(XpadNew[jj]-x_mm)<0.0001 && 
abs(YpadNew[jj]-y_mm)<0.0001) {
                        fitbool = true;
                        indexfill=jj;
                        break;
                    }
                }
                
                // Check if New channel is of interest
                //(if so, we read the Q(t) signal, and we should fill the 
vectors w/ t & q information after fitting E(T))
                if( fitbool==true ) {
                    for(Int_t j=0; j<minos->GetNData(); j++) {
                        if(minos->GetCalibValue(j)>=0){
                            
hfit->SetBinContent(hfit->FindBin(minos->GetCalibTime(j)), 
minos->GetCalibValue(j)+250);
                            //cerr << "    * t=" << minos->GetCalibTime(j) << " 
* q=" << minos->GetCalibValue(j) << endl;
                        }
                    }
                    
                    // Fitting the hfit histogram of last channel if not empty
                    if(hfit->GetSumOfWeights()>0) {
                        hfit->GetXaxis()->SetRange(0,510);
                        hfit_max = hfit->GetMaximum();
                        hfit_max_T = hfit->GetMaximumBin();
                        T_min=-1;
                        T_max=-1;
                        
                        // Find the T_min & T_max limits of the signal non zero
                        for(int h=hfit_max_T;h>0;h--) {
                            if(T_min == -1 && (hfit->GetBinContent(h))<=250 ) {
                                T_min = h;
                                break;
                            }
                        }
                        for(int h=hfit_max_T;h<510;h++) {
                            if(T_max == -1 && (hfit->GetBinContent(h))==0 ) {
                                T_max = h;
                                break;
                            }
                        }
                        //Take only 1.5*Tshaping before the max if other signals 
before...
                        if((hfit_max_T-3.5*(Tshaping/TimeBinElec)) > T_min) 
T_min = hfit_max_T-2*Tshaping/TimeBinElec;
                        if((hfit_max_T+10) < T_max || T_max==-1) T_max = 
hfit_max_T+10.;
                        
                        T_min = max(T_min,0.);
                        if(T_max>510) T_max = 510;
                        
                        
                        // Set fit parameters
                        fit_function->SetParameter(0, hfit_max-250.);
                        fit_function->SetParameter(1,hfit_max_T - 
Tshaping/TimeBinElec);
                        fit_function->SetParameter(2, Tshaping/TimeBinElec);
                        fit_function->SetParLimits(0,0,100000);
                        fit_function->SetParLimits(1,-20,512);
                        fit_function->SetParLimits(2,0,512);
                        
                        // Fit of the signal within the range defined: [T_min, 
T_max]
                        fit2DStatus = 
hfit->Fit(fit_function,"Q","",T_min,T_max);
                        //gStyle->SetOptFit(1111);
                        
                        double fit_function_max = 0., fit_function_Tpad = 0.;
                        
                        if(fit2DStatus==0) {
                            TF1 *fit_result = hfit->GetFunction("fit_function");
                            Chi2 = fit_result->GetChisquare();
                            fit_function_max = fit_function->GetMaximum();
                            fit_function_Tpad = fit_function->GetParameter(1);
                        }
                	    
                        //attribute q_pad and z_mm value
                        if(fit2DStatus!=0 || fit_function_max<=20. || 
fit_function_max>=100000. || fit_function_Tpad<=0.15 || fit_function_Tpad>=513. 
|| fit_function->GetParameter(2)<=0.15 || fit_function->GetParameter(2)>=513.) {
                            //cout << "NOT CORRECTLY FITTED !!!!!!! chi2 = " << 
vectout->chi2.back() << endl;
                            q_pad = hfit_max-250.;
                            z_mm = -10000;
                        }
                        else {
                            // Add to the variables the fit parameters
                            t_pad = fit_function_Tpad;//trigger time
                            z_mm = ((t_pad*TimeBinElec-DelayTrig)*VDrift); 
//time bin*(ns/us)*vdrift(mm/ns) =>mm
                            q_pad = fit_function_max-250.;  // Max charge/event
                        }
                        
                        // Comment out for saving histograms of Q(t)
                        // TH1F *hfit_clone = 
(TH1F*)hfit->Clone(Form("E(T)_%d_%f_%f",neve,x_mm,y_mm));
                        //hfit_clone->Write();
                        //cout << " Histogram results: " << i << " --- @ (" << 
lastx_mm << ", " << lasty_mm << "); q_pad = " << q_pad << endl;
                        minosfitdata = 
(TMinosClust*)fitdata.ConstructedAt(indexfill);
                        minosfitdata->Set(XpadNew[indexfill], 
YpadNew[indexfill], t_pad*TimeBinElec, z_mm , q_pad, clusternbr[indexfill], 
clusterpads[indexfill], Chi2);
                        //Fill the z and q information for next steps (3D Hough 
filter & 3D fit weighted by charge)
                        ZpadNew[indexfill] = z_mm;
                        QpadNew[indexfill] = q_pad;
                        
                    }//end if histogram not empty
                    
                    //cerr << "Event " << neve << "::: x=" << XpadNew[indexfill] 
<< ", y=" << YpadNew[indexfill] << ", Qmax=" << q_pad << ", t=" << t_pad << ", 
z=" << z_mm << endl;
                    
                }// end if fitbool==true
                else continue;
            }//END of entries in tclonesarray for the event
            
            
            
            
            //4:MINOS://////////// Filtering the tracks off possible noise with 
Hough3D (3*2D planes) ////////////:MINOS://
            
            padsleft2 = XpadNew.size();
            for(unsigned int i=0;i<(XpadNew.size());i++)
            {
                
                if( (cluster_temp==int(clusternbr[i])) && i==(XpadNew.size() - 
1) && clusterpads[i]>=10 && clusterringbool[i]==1 && ZpadNew[i]>-10000 && 
ZpadNew[i]<=320)
            	{
					//cout<<minosdata->x_mm<<" "<<minosdata->z_mm<<" 
"<<minosdata->Phi<<endl;
					//cout<<"trackNbr "<<trackNbr<<endl;
					xin.push_back(XpadNew[i]);
					yin.push_back(YpadNew[i]);
					zin.push_back(ZpadNew[i]);
					qin.push_back(QpadNew[i]);
					npoint_temp++;
                }
                
                //cerr << "Event " << neve << ", xpadnew=" << i << endl;
                if(xin.size()>0 && ((cluster_temp!=int(clusternbr[i]) && i!=0) 
|| i==(XpadNew.size() - 1)))
                {
                    //   			Hough_canvas.push_back(new 
TCanvas(Form("Event%d_cluster%d", entry, cluster_temp), 
Form("Event%d_cluster%d", entry, cluster_temp)));
                    //				Hough_3D(Hough_canvas.back(), &xin, &yin, 
&zin, &qin, &xout, &yout, &zout, &qout);
                    Tracking_functions->Hough_3D(&xin, &yin, &zin, &qin, &xout, 
&yout, &zout, &qout);
                    //				Hough_canvas.back()->Write();
                    for(unsigned int ij=0; ij<xout.size();ij++)
                    {
                    	if(zout[ij]>zmax) zmax = zout[ij];
                    	
ringtouch[int((sqrt(xout[ij]*xout[ij]+yout[ij]*yout[ij])-45.2)/2.1)]++;
                    }
                    for(int ko=0; ko<18; ko++)
                    {
                        if(ringtouch[ko]>0) ringsum++;
                    }
                	if(zmax>290) ringsum=16;
                	if(xout.size()>10 && ringsum>=15)
                	{
                    	trackNbr_FINAL++;
                    	if(trackNbr_FINAL==1)
                    	{
                        	cluster1 = cluster_temp;
                        	for(unsigned int ij=0; ij<xout.size(); ij++)
                        	{
                                grxz_1->SetPoint(npoint1,zout[ij],xout[ij]);
                                gryz_1->SetPoint(npoint1,zout[ij],yout[ij]);
                                minosdata_result = 
(TMinosResult*)data_result.ConstructedAt(array_final);
                                minosdata_result->Set(xout[ij], yout[ij], 
zout[ij], qout[ij], 1, xout.size(), zmax);
                                array_final++;
                                npoint1++;
                        	}
						}
                    	else if(trackNbr_FINAL==2)
                    	{
                        	cluster2 = cluster_temp;
                        	for(unsigned int ij=0; ij<xout.size(); ij++)
                        	{
                                grxz_2->SetPoint(npoint2,zout[ij],xout[ij]);
                                gryz_2->SetPoint(npoint2,zout[ij],yout[ij]);
                                minosdata_result = 
(TMinosResult*)data_result.ConstructedAt(array_final);
                                minosdata_result->Set(xout[ij], yout[ij], 
zout[ij], qout[ij], 2, xout.size(), zmax);
                                array_final++;
                                npoint2++;
                        	}
                    	}
                	}
                    //cout<<"Evt "<<entry<<" data point "<<i<<" track Nr = 
"<<trackNbr_FINAL<<" xout.size ="<<xout.size()<<" data cluster 
"<<minosdata->n_Cluster<<endl;
                    
                	xin.clear();
                	yin.clear();
                	zin.clear();
                	qin.clear();
                	xout.clear();
               		yout.clear();
                	zout.clear();
                	qout.clear();
                	npoint_temp=0;
                	ringsum=0;// angle between 1st track and z axis in 3D in 
degrees
                	zmax=0.;
                	for(int ko=0; ko<18; ko++) ringtouch[ko] = 0;
                    
            	}
                
            	cluster_temp = clusternbr[i];
                
            	if(!(clusterpads[i]>=10 && clusterringbool[i]==1 && 
ZpadNew[i]>-10000 && ZpadNew[i]<=320)) continue;
            	else
            	{
					//cout<<minosdata->x_mm<<" "<<minosdata->z_mm<<" 
"<<minosdata->Phi<<endl;
					//cout<<"trackNbr "<<trackNbr<<endl;
					xin.push_back(XpadNew[i]);
					yin.push_back(YpadNew[i]);
					zin.push_back(ZpadNew[i]);
					qin.push_back(QpadNew[i]);
					npoint_temp++;
            	}
                
      		}//end of loop on pads
            
            
            //5:MINOS://////////// Fitting the filtered tracks in 3D (weight by 
charge, TMinuit) ////////////:MINOS://
            
			//For 1 track found or more (less than 5 in total)
        	//if(trackNbr_FINAL==1 || trackNbr_FINAL==2)
        	if(trackNbr_FINAL>=1)
        	{
                
				//////////Minimization in 3D to reconstruct track lines
        		allevt_2pfiltered++;
				padsleft2 = padsleft2 - npoint1 - npoint2;
                
    			Double_t pStart_1[4]={0,1,0,1};
    			Double_t pStart_2[4]={0,1,0,1};
    			Double_t chi1[2],chi2[2];
        		Int_t fitStatus[2];
            	
				// 1st track fitting _ 1p in MINOS _ in TMinuit
        		//cout<<"start Minuit "<<endl;
            	min = new TMinuit(4);
            	min->SetPrintLevel(-1);
            	Double_t arglist[10];
            	arglist[0] = 3;
            	Int_t iflag;
            	int nvpar,nparx;
            	double amin,edm, errdef;
            	double chi2res1, chi2res2;
            	
            	
            	
    			Tracking_functions->FindStart(pStart_1,chi1,fitStatus, 
grxz_1,gryz_1);
        		min->SetFCN(SumDistance1);
                // Set starting values and step sizes for parameters
        		min->mnparm(0,"x0",pStart_1[0],0.1,-500,500,iflag);
        		min->mnparm(1,"Ax",pStart_1[1],0.1,-10,10,iflag);
        		min->mnparm(2,"y0",pStart_1[2],0.1,-500,500,iflag);
        		min->mnparm(3,"Ay",pStart_1[3],0.1,-10,10,iflag);
        		arglist[0] = 100; // number of function calls
        		arglist[1] = 0.000001; // tolerance
        		min->mnexcm("MIGRAD",arglist,2,iflag); // minimization with 
MIGRAD
            	
        		min->mnstat(amin,edm,errdef,nvpar,nparx,iflag);  //returns 
current status of the minimization
        		// get fit parameters
        		for (int i = 0; i <4; i++) {
           			min->GetParameter(i,parFit_1[i],errFit_1[i]);
                    //cout << "before fit: " << i << "= " << pStart_1[i] << 
endl;
                    //cout << "after fit:  " << i << "= " << parFit_1[i] << " ; 
" << err_1[i] << endl;
        		}
                //AC temporaty fix before mapping phi angle of TPC
		//parFit_2_F8[0]=0;parFit_2_F8[1]=0;parFit_2_F8[2]=0;parFit_2_F8[3]=0;
		//cout<<"before if F8PPAC, x_F8PPAC = "<< x_F8PPAC<<"	y_F8PPAC "<< 
y_F8PPAC<<endl;

       //Now account for possible rotation angle between MINOS and 
F8PPAC/beamline
       double minosrot=TMath::DegToRad()*-30.0;
       //cout<<TMath::Cos(minosrot)<<"  "<<TMath::Sin(minosrot)<<endl;
       parFit_1r[0]=cos(minosrot)*parFit_1[0]-sin(minosrot)*parFit_1[2];
       parFit_1r[1]=cos(minosrot)*parFit_1[1]-sin(minosrot)*parFit_1[3];
       parFit_1r[2]=sin(minosrot)*parFit_1[0]+cos(minosrot)*parFit_1[2];
       parFit_1r[3]=sin(minosrot)*parFit_1[1]+cos(minosrot)*parFit_1[3];

      //cout<<parFit_1[0]<< "  "<<parFit_1r[0]<<endl;
      
                // Vertex with beam tracking : F8PPAC / DSSSD
		// FF (30/05/2015 2:30 AM)  commented if because it is never verified 
due to some bad PPAC calculation...              
		//if(abs(x_F8PPAC+10000.)>0.0001 && abs(y_F8PPAC+10000.)>0.0001) {
		    //cout<<"just before fit"<<parFit_2_F8[0]<< "  
"<<parFit_2_F8[2]<<endl;
			Tracking_functions->vertex(parFit_1r, parFit_2_F8, x_vertexPPAC, 
y_vertexPPAC, z_vertexPPAC);
			//cout<<"in if F8PPAC, xyz_vertexPPAC = ("<< x_vertexPPAC<<", "<< 
y_vertexPPAC<<", "<<z_vertexPPAC<<")"<<endl;
		//}         
		      
				//CASE: only 1p detected in MINOS - get PPAC/DSSSD tracking
				if(trackNbr_FINAL==1)
				{
					//PPAC - find 3Dline param.
					// FF (30/05/2015 2:30 AM)  commented if because it is never 
verified due to some bad PPAC calculation... 
					//if(abs(x_F8PPAC+10000.)>0.0001 && 
abs(y_F8PPAC+10000.)>0.0001) {
						for (int iii=0; iii<4; iii++) parFit_2[iii] = 
parFit_2_F8[iii];

						
					//}
                    
				}
                
				//CASE: 2p or more detected in MINOS - get the second proton 3D 
line track fitted
				else
				{
    				
Tracking_functions->FindStart(pStart_2,chi2,fitStatus,grxz_2,gryz_2);
        			min->SetFCN(SumDistance2);
        			// Set starting values and step sizes for parameters
        			min->mnparm(0,"x0",pStart_2[0],0.1,-500,500,iflag);
        			min->mnparm(1,"Ax",pStart_2[1],0.1,-10,10,iflag);
        			min->mnparm(2,"y0",pStart_2[2],0.1,-500,500,iflag);
        			min->mnparm(3,"Ay",pStart_2[3],0.1,-10,10,iflag);
        			arglist[0] = 100; // number of function calls
        			arglist[1] = 0.000001; // tolerance
        			min->mnexcm("MIGRAD",arglist,2,iflag); // minimization with 
MIGRAD
                    
        			min->mnstat(amin,edm,errdef,nvpar,nparx,iflag);  //returns 
current status of the minimization
        			// get fit parameters
        			for (int i = 0; i <4; i++) {
                        min->GetParameter(i,parFit_2[i],errFit_2[i]);
                        //cout << "before fit: " << i << "= " << pStart_2[i] << 
endl;
                        //cout << "after fit:  " << i << "= " << parFit_2[i] << 
" ; " << endl;	
        			}
            	}
                
            	int sumerr_2=0;
            	int sumerr_1=0;
                
            	for(int index=0; index<4; index++) {
                	sumerr_2+=errFit_2[index];
                	sumerr_1+=errFit_1[index];
            	}
            	
            	//cout<< "Entry " << evtOrig << " --- chi2 for 1st cluster=" << 
chi1[0]/npoint1 <<", "<< chi1[1]/npoint1 <<" -- chi2 for 2nd cluster="<<  
chi2[0]/npoint2 <<", "<< chi2[1]/npoint2 <<endl;
                
            	if(npoint1>0) chi2res1 = (chi1[0]+chi1[1])/npoint1;
            	else chi2res1 = -1000;
            	if(npoint2>0) chi2res2 = (chi2[0]+chi2[1])/npoint2;
            	else chi2res2 = -1000;
                
                //cout<<"parFit_2_F8[0] "<<parFit_2_F8[0]<<" <parFit_2[0]  
"<<parFit_2[0]<<endl;
                
                //Also rotate the second Track in the TPC (if it exists)
                if(trackNbr_FINAL>=2){
                
parFit_2r[0]=cos(minosrot)*parFit_2[0]-sin(minosrot)*parFit_2[2];
                
parFit_2r[1]=cos(minosrot)*parFit_2[1]-sin(minosrot)*parFit_2[3];
                
parFit_2r[2]=sin(minosrot)*parFit_2[0]+cos(minosrot)*parFit_2[2];
                
parFit_2r[3]=sin(minosrot)*parFit_2[1]+cos(minosrot)*parFit_2[3];
                
                //cout<<"after rotation  parFit_2_F8[0] "<<parFit_2_F8[0]<<" 
<parFit_2[0]  "<<parFit_2[0]<< "parFit_2r[0]  "<<parFit_2r[0]<<endl;
                }
				// Reconstruct vertex from 2 tracks
            	Tracking_functions->vertex(parFit_1r, parFit_2r, x_vertex, 
y_vertex, z_vertex);
                
            	r_vertex = sqrt(x_vertex*x_vertex + y_vertex*y_vertex);
                
				thetaz1 = acos(1/sqrt(1 + parFit_1r[1]*parFit_1r[1] + 
parFit_1r[3]*parFit_1r[3]))*180./PI;
				if(thetaz1>180) cerr << "thetaz1 cond &&&" << endl;
				thetaz2 = acos(1/sqrt(1 + parFit_2r[1]*parFit_2r[1] + 
parFit_2[3]*parFit_2r[3]))*180./PI;
				if(thetaz2>180) cerr << "thetaz2 cond @@@@" << endl;
                
				//phi_vertex: angle in 3D between 2 tracks
                phi_vertex = acos((parFit_1r[1]*parFit_2r[1] + 
parFit_1r[3]*parFit_2r[3]+1)/(sqrt(parFit_1r[1]*parFit_1r[1]+parFit_1r[3]*
parFit_1r[3]+1)*sqrt(parFit_2r[1]*parFit_2r[1]+parFit_2r[3]*parFit_2r[3]+1)))*
180./PI;
                
            	if(grxz_1->GetN()>1 && grxz_2->GetN()>1) {
    				if(z_vertex>=-500 && z_vertex<=500){
                    	allevt_2pvertex++;
                	}
            	}
                
      			delete min;
                
        	}// end if trackNbr_FINAL>=1
            
      		delete grxz_1;
      		delete gryz_1;
      		delete grxz_2;
      		delete gryz_2;
            
        	xin.clear();
        	yin.clear();
        	zin.clear();
        	qin.clear();
        	xout.clear();
        	yout.clear();
        	zout.clear();
        	qout.clear();
        	cluster.clear();
			//Hough_canvas.clear();
            
		}//loop for E(T) fits for less than 5 track evts
			  //Delta 
		
		//cout<<"Just before fill  "<<tpla->GetTLRaw()<<endl;

		tree->Fill();  // fill the tree in ALL cases to obtain the same number 
of evts for DALI2 analysis
        


		estore->ClearData();
    	neve ++;
  	}
    
    out<< rootfile << "  "<<InTrig<<endl;
    tree->Print();
    cout<<"Write..."<<endl;
    fout->Write();
    cout<<"Close..."<<endl;
    fout->Close();
    cout<<"Conversion to Root done!"<<endl;
    
    time(&stop);
    printf("Elapsed time: %.1f seconds\n",difftime(stop,start));
    
    
    return 0;
    
}

/// Functions to be minimized
void SumDistance1(int &, double *, double & sum, double * par,  int) {
    int nused=0;
    double qtot=0;
    sum = 0;
    //cout<<"sum "<<sum<<" over "<<npoints<<endl;
    //double factor;
    //cout<<"*************after fit "<<endl;
    for(int i=0; i<data_result.GetEntriesFast(); i++)
    {
    	minosdata_result = (TMinosResult*)data_result.At(i);
    	if(minosdata_result->n_Cluster==1)
    	{
        	float x=minosdata_result->x_mm;
        	float y=minosdata_result->y_mm;
        	float z=minosdata_result->z_mm;
        	float q=minosdata_result->Chargemax;
        	//if(nused<2)cout<<minosdata_result->n_Cluster<<" "<<x<<" "<<y<<" 
"<<z<<" "<<q<<endl;
        	double d = Tracking_functions->distance2(x, y, z, par);
      		sum += d*q;
      		nused++;
            qtot+=q;
    	}
    }
    //sum/=nused;
    sum/=qtot;
}

void SumDistance2(int &, double *, double & sum, double * par, int) {
    int nused=0;
    double qtot=0;
    sum = 0;
    //cout<<"sum "<<sum<<" over "<<npoints<<endl;
    //double factor;
    //cout<<"*************after fit "<<endl;
    for(int i=0; i<data_result.GetEntriesFast(); i++)
    {
    	minosdata_result = (TMinosResult*)data_result.At(i);
    	if(minosdata_result->n_Cluster==2) 
    	{
        	float x=minosdata_result->x_mm;
        	float y=minosdata_result->y_mm;
        	float z=minosdata_result->z_mm;
        	float q=minosdata_result->Chargemax;
        	//if(nused<2)cout<<minosdata_result->n_Cluster<<" "<<x<<" "<<y<<" 
"<<z<<" "<<q<<endl;
        	double d = Tracking_functions->distance2(x, y, z, par);
      		sum += d*q;
      		nused++;
            qtot+=q;
    	}
    }
    //sum/=nused;
    sum/=qtot;
 }   

