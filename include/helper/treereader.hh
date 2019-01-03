//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 25 09:20:01 2018 by ROOT version 6.12/06
// from TTree tree/tree
// found on file: build/output/out.root
//////////////////////////////////////////////////////////

#ifndef treereader_h
#define treereader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TNamed.h"
#include "TVectorT.h"
#include <iostream>

using std::vector;

class treereader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxEventInfo = 1;
   static constexpr Int_t kMaxBigRIPSPPAC = 37;
   static constexpr Int_t kMaxBigRIPSPlastic = 4;
   static constexpr Int_t kMaxBigRIPSIC = 2;
   static constexpr Int_t kMaxBigRIPSFocalPlane = 14;
   static constexpr Int_t kMaxBigRIPSRIPS = 4;
   static constexpr Int_t kMaxBigRIPSTOF = 2;
   static constexpr Int_t kMaxBigRIPSBeam = 5;

   // Declaration of leaf types
   Int_t           EventInfo_;
    UInt_t          EventInfo_fUniqueID[kMaxEventInfo];   //[EventInfo_]
    UInt_t          EventInfo_fBits[kMaxEventInfo];   //[EventInfo_]
    TString         EventInfo_fName[kMaxEventInfo];
    TString         EventInfo_fTitle[kMaxEventInfo];
    Int_t           EventInfo_runnumber[kMaxEventInfo];   //[EventInfo_]
    Int_t           EventInfo_eventnumber[kMaxEventInfo];   //[EventInfo_]
    TString         EventInfo_subsysname[kMaxEventInfo];
    ULong64_t       EventInfo_timestamp[kMaxEventInfo];   //[EventInfo_]
    Int_t           EventInfo_comp_val[kMaxEventInfo];   //[EventInfo_]
    UInt_t          EventInfo_fBit[kMaxEventInfo];   //[EventInfo_]
    Int_t           BigRIPSPPAC_;
    UInt_t          BigRIPSPPAC_fUniqueID[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    UInt_t          BigRIPSPPAC_fBits[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_id[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fpl[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    TString         BigRIPSPPAC_name[kMaxBigRIPSPPAC];
    Int_t           BigRIPSPPAC_fDataState[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_xzpos[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_yzpos[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fTX1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fTX2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fTY1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fTY2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fTARaw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fQX1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fQX2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fQY1Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fQY2Raw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPPAC_fQARaw[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTX1[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTX2[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTY1[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTY2[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTA[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTSumX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTSumY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTDiffX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fTDiffY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Double_t        BigRIPSPPAC_fY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Bool_t          BigRIPSPPAC_fFiredX[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Bool_t          BigRIPSPPAC_fFiredY[kMaxBigRIPSPPAC];   //[BigRIPSPPAC_]
    Int_t           BigRIPSPlastic_;
    UInt_t          BigRIPSPlastic_fUniqueID[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    UInt_t          BigRIPSPlastic_fBits[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_id[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fpl[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    TString         BigRIPSPlastic_name[kMaxBigRIPSPlastic];
    Int_t           BigRIPSPlastic_fDataState[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_zposition[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_zoffset[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_tcalL[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_tcalR[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fTLRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fTRRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fQLRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fQRRaw[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fQTCLRawWidth[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fQTCRRawWidth[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fQTCLRawStart[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSPlastic_fQTCRRawStart[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    vector<int>     BigRIPSPlastic_fTLRawArray[kMaxBigRIPSPlastic];
    vector<int>     BigRIPSPlastic_fTRRawArray[kMaxBigRIPSPlastic];
    Double_t        BigRIPSPlastic_fTime[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_fTimeL[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_fTimeR[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_fTimeLSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_fTimeRSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Double_t        BigRIPSPlastic_fTimeSlew[kMaxBigRIPSPlastic];   //[BigRIPSPlastic_]
    Int_t           BigRIPSIC_;
    UInt_t          BigRIPSIC_fUniqueID[kMaxBigRIPSIC];   //[BigRIPSIC_]
    UInt_t          BigRIPSIC_fBits[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Int_t           BigRIPSIC_id[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Int_t           BigRIPSIC_fpl[kMaxBigRIPSIC];   //[BigRIPSIC_]
    TString         BigRIPSIC_name[kMaxBigRIPSIC];
    Int_t           BigRIPSIC_fDataState[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_zcoef[kMaxBigRIPSIC][2];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_ionpair[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Int_t           BigRIPSIC_nhitchannel[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Int_t           BigRIPSIC_fADC[kMaxBigRIPSIC][32];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_fEnergy[kMaxBigRIPSIC][32];   //[BigRIPSIC_]
    Int_t           BigRIPSIC_fTDC[kMaxBigRIPSIC][32];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_fTime[kMaxBigRIPSIC][32];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_fRawADCSqSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_fRawADCAvSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_fCalMeVSqSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Double_t        BigRIPSIC_fCalMeVAvSum[kMaxBigRIPSIC];   //[BigRIPSIC_]
    Int_t           BigRIPSFocalPlane_;
    UInt_t          BigRIPSFocalPlane_fUniqueID[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    UInt_t          BigRIPSFocalPlane_fBits[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Int_t           BigRIPSFocalPlane_id[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Int_t           BigRIPSFocalPlane_fpl[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    TString         BigRIPSFocalPlane_name[kMaxBigRIPSFocalPlane];
    Int_t           BigRIPSFocalPlane_fDataState[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    TVectorT<double> BigRIPSFocalPlane_opt_vector[kMaxBigRIPSFocalPlane];
    Double_t        BigRIPSFocalPlane_X[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Double_t        BigRIPSFocalPlane_A[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Double_t        BigRIPSFocalPlane_Y[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Double_t        BigRIPSFocalPlane_B[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Int_t           BigRIPSFocalPlane_nfired_ppacx[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Int_t           BigRIPSFocalPlane_nfired_ppacy[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Double_t        BigRIPSFocalPlane_zpos[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Double_t        BigRIPSFocalPlane_zpos_offset[kMaxBigRIPSFocalPlane];   //[BigRIPSFocalPlane_]
    Int_t           BigRIPSRIPS_;
    UInt_t          BigRIPSRIPS_fUniqueID[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    UInt_t          BigRIPSRIPS_fBits[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Int_t           BigRIPSRIPS_id[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Int_t           BigRIPSRIPS_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    TString         BigRIPSRIPS_name[kMaxBigRIPSRIPS];
    Int_t           BigRIPSRIPS_fDataState[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Int_t           BigRIPSRIPS_upstream_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Int_t           BigRIPSRIPS_downstream_fpl[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Double_t        BigRIPSRIPS_center_brho[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Double_t        BigRIPSRIPS_brho[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Double_t        BigRIPSRIPS_length[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    TMatrixT<double> BigRIPSRIPS_matrix[kMaxBigRIPSRIPS];
    Double_t        BigRIPSRIPS_delta[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    Double_t        BigRIPSRIPS_angle[kMaxBigRIPSRIPS];   //[BigRIPSRIPS_]
    TString         BigRIPSRIPS_dipolename[kMaxBigRIPSRIPS];
    Int_t           BigRIPSTOF_;
    UInt_t          BigRIPSTOF_fUniqueID[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    UInt_t          BigRIPSTOF_fBits[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Int_t           BigRIPSTOF_id[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Int_t           BigRIPSTOF_fpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    TString         BigRIPSTOF_name[kMaxBigRIPSTOF];
    Int_t           BigRIPSTOF_fDataState[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Double_t        BigRIPSTOF_tof[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Double_t        BigRIPSTOF_clight[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Double_t        BigRIPSTOF_length[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Double_t        BigRIPSTOF_ulength[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Double_t        BigRIPSTOF_dlength[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    TString         BigRIPSTOF_upstream_plname[kMaxBigRIPSTOF];
    TString         BigRIPSTOF_downstream_plname[kMaxBigRIPSTOF];
    Int_t           BigRIPSTOF_upstream_plfpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Int_t           BigRIPSTOF_downstream_plfpl[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Double_t        BigRIPSTOF_time_offset[kMaxBigRIPSTOF];   //[BigRIPSTOF_]
    Int_t           BigRIPSBeam_;
    UInt_t          BigRIPSBeam_fUniqueID[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    UInt_t          BigRIPSBeam_fBits[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Int_t           BigRIPSBeam_id[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Int_t           BigRIPSBeam_fpl[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    TString         BigRIPSBeam_name[kMaxBigRIPSBeam];
    Int_t           BigRIPSBeam_fDataState[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Double_t        BigRIPSBeam_brho[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Double_t        BigRIPSBeam_aoq[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Double_t        BigRIPSBeam_zet[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Double_t        BigRIPSBeam_beta[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Double_t        BigRIPSBeam_clight[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Double_t        BigRIPSBeam_mnucleon[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    Int_t           BigRIPSBeam_nrips[kMaxBigRIPSBeam];   //[BigRIPSBeam_]
    TString         BigRIPSBeam_ripsname[2][kMaxBigRIPSBeam];
    TString         BigRIPSBeam_tofname[kMaxBigRIPSBeam];
    TString         BigRIPSBeam_icname[kMaxBigRIPSBeam];
    vector<double>  *MinosClustX;
    vector<double>  *MinosClustY;
    vector<double>  *MinosClustZ;
    vector<double>  *MinosClustT;
    vector<double>  *MinosClustQ;
    vector<double>  *MinosClustNCl;
    vector<double>  *MinosClustNP;
    vector<double>  *MinosClustChi2;
    vector<double>  *MinosResultX;
    vector<double>  *MinosResultY;
    vector<double>  *MinosResultZ;
    vector<double>  *MinosResultQ;
    vector<double>  *MinosResultNCl;
    vector<double>  *MinosResultNP;
    vector<double>  *MinosResultChi2;
    Int_t           trackNbr;
    Int_t           trackNbr_final;
    Int_t           evtOrig;
    Double_t        VDrift;
    Double_t        DelayTrig;
    Double_t        x_vertex;
    Double_t        y_vertex;
    Double_t        z_vertex;
    Double_t        r_vertex;
    Double_t        phi_vertex;
    Double_t        thetaz1;
    Double_t        thetaz2;
    vector<double>  *parFit_1;
    vector<double>  *parFit_2;
    Double_t        F3X;
    Double_t        F3A;
    Double_t        F5X;
    Double_t        F5A;
    Double_t        F7X;
    Double_t        F7A;
    Double_t        F8X;
    Double_t        F8A;
    Double_t        F9X;
    Double_t        F9A;
    Double_t        F11X;
    Double_t        F11A;
    Int_t           fgoodppacfocus[12];
    Int_t           fgoodppacfocusor[12];

   //List of event-enum
   UInt_t          currenteventnum =0;

   // List of branches
   TBranch        *b_EventInfo_;   //!
    TBranch        *b_EventInfo_fUniqueID;   //!
    TBranch        *b_EventInfo_fBits;   //!
    TBranch        *b_EventInfo_fName;   //!
    TBranch        *b_EventInfo_fTitle;   //!
    TBranch        *b_EventInfo_runnumber;   //!
    TBranch        *b_EventInfo_eventnumber;   //!
    TBranch        *b_EventInfo_subsysname;   //!
    TBranch        *b_EventInfo_timestamp;   //!
    TBranch        *b_EventInfo_comp_val;   //!
    TBranch        *b_EventInfo_fBit;   //!
    TBranch        *b_BigRIPSPPAC_;   //!
    TBranch        *b_BigRIPSPPAC_fUniqueID;   //!
    TBranch        *b_BigRIPSPPAC_fBits;   //!
    TBranch        *b_BigRIPSPPAC_id;   //!
    TBranch        *b_BigRIPSPPAC_fpl;   //!
    TBranch        *b_BigRIPSPPAC_name;   //!
    TBranch        *b_BigRIPSPPAC_fDataState;   //!
    TBranch        *b_BigRIPSPPAC_xzpos;   //!
    TBranch        *b_BigRIPSPPAC_yzpos;   //!
    TBranch        *b_BigRIPSPPAC_fTX1Raw;   //!
    TBranch        *b_BigRIPSPPAC_fTX2Raw;   //!
    TBranch        *b_BigRIPSPPAC_fTY1Raw;   //!
    TBranch        *b_BigRIPSPPAC_fTY2Raw;   //!
    TBranch        *b_BigRIPSPPAC_fTARaw;   //!
    TBranch        *b_BigRIPSPPAC_fQX1Raw;   //!
    TBranch        *b_BigRIPSPPAC_fQX2Raw;   //!
    TBranch        *b_BigRIPSPPAC_fQY1Raw;   //!
    TBranch        *b_BigRIPSPPAC_fQY2Raw;   //!
    TBranch        *b_BigRIPSPPAC_fQARaw;   //!
    TBranch        *b_BigRIPSPPAC_fTX1;   //!
    TBranch        *b_BigRIPSPPAC_fTX2;   //!
    TBranch        *b_BigRIPSPPAC_fTY1;   //!
    TBranch        *b_BigRIPSPPAC_fTY2;   //!
    TBranch        *b_BigRIPSPPAC_fTA;   //!
    TBranch        *b_BigRIPSPPAC_fTSumX;   //!
    TBranch        *b_BigRIPSPPAC_fTSumY;   //!
    TBranch        *b_BigRIPSPPAC_fTDiffX;   //!
    TBranch        *b_BigRIPSPPAC_fTDiffY;   //!
    TBranch        *b_BigRIPSPPAC_fX;   //!
    TBranch        *b_BigRIPSPPAC_fY;   //!
    TBranch        *b_BigRIPSPPAC_fFiredX;   //!
    TBranch        *b_BigRIPSPPAC_fFiredY;   //!
    TBranch        *b_BigRIPSPlastic_;   //!
    TBranch        *b_BigRIPSPlastic_fUniqueID;   //!
    TBranch        *b_BigRIPSPlastic_fBits;   //!
    TBranch        *b_BigRIPSPlastic_id;   //!
    TBranch        *b_BigRIPSPlastic_fpl;   //!
    TBranch        *b_BigRIPSPlastic_name;   //!
    TBranch        *b_BigRIPSPlastic_fDataState;   //!
    TBranch        *b_BigRIPSPlastic_zposition;   //!
    TBranch        *b_BigRIPSPlastic_zoffset;   //!
    TBranch        *b_BigRIPSPlastic_tcalL;   //!
    TBranch        *b_BigRIPSPlastic_tcalR;   //!
    TBranch        *b_BigRIPSPlastic_fTLRaw;   //!
    TBranch        *b_BigRIPSPlastic_fTRRaw;   //!
    TBranch        *b_BigRIPSPlastic_fQLRaw;   //!
    TBranch        *b_BigRIPSPlastic_fQRRaw;   //!
    TBranch        *b_BigRIPSPlastic_fQTCLRawWidth;   //!
    TBranch        *b_BigRIPSPlastic_fQTCRRawWidth;   //!
    TBranch        *b_BigRIPSPlastic_fQTCLRawStart;   //!
    TBranch        *b_BigRIPSPlastic_fQTCRRawStart;   //!
    TBranch        *b_BigRIPSPlastic_fTLRawArray;   //!
    TBranch        *b_BigRIPSPlastic_fTRRawArray;   //!
    TBranch        *b_BigRIPSPlastic_fTime;   //!
    TBranch        *b_BigRIPSPlastic_fTimeL;   //!
    TBranch        *b_BigRIPSPlastic_fTimeR;   //!
    TBranch        *b_BigRIPSPlastic_fTimeLSlew;   //!
    TBranch        *b_BigRIPSPlastic_fTimeRSlew;   //!
    TBranch        *b_BigRIPSPlastic_fTimeSlew;   //!
    TBranch        *b_BigRIPSIC_;   //!
    TBranch        *b_BigRIPSIC_fUniqueID;   //!
    TBranch        *b_BigRIPSIC_fBits;   //!
    TBranch        *b_BigRIPSIC_id;   //!
    TBranch        *b_BigRIPSIC_fpl;   //!
    TBranch        *b_BigRIPSIC_name;   //!
    TBranch        *b_BigRIPSIC_fDataState;   //!
    TBranch        *b_BigRIPSIC_zcoef;   //!
    TBranch        *b_BigRIPSIC_ionpair;   //!
    TBranch        *b_BigRIPSIC_nhitchannel;   //!
    TBranch        *b_BigRIPSIC_fADC;   //!
    TBranch        *b_BigRIPSIC_fEnergy;   //!
    TBranch        *b_BigRIPSIC_fTDC;   //!
    TBranch        *b_BigRIPSIC_fTime;   //!
    TBranch        *b_BigRIPSIC_fRawADCSqSum;   //!
    TBranch        *b_BigRIPSIC_fRawADCAvSum;   //!
    TBranch        *b_BigRIPSIC_fCalMeVSqSum;   //!
    TBranch        *b_BigRIPSIC_fCalMeVAvSum;   //!
    TBranch        *b_BigRIPSFocalPlane_;   //!
    TBranch        *b_BigRIPSFocalPlane_fUniqueID;   //!
    TBranch        *b_BigRIPSFocalPlane_fBits;   //!
    TBranch        *b_BigRIPSFocalPlane_id;   //!
    TBranch        *b_BigRIPSFocalPlane_fpl;   //!
    TBranch        *b_BigRIPSFocalPlane_name;   //!
    TBranch        *b_BigRIPSFocalPlane_fDataState;   //!
    TBranch        *b_BigRIPSFocalPlane_opt_vector;   //!
    TBranch        *b_BigRIPSFocalPlane_X;   //!
    TBranch        *b_BigRIPSFocalPlane_A;   //!
    TBranch        *b_BigRIPSFocalPlane_Y;   //!
    TBranch        *b_BigRIPSFocalPlane_B;   //!
    TBranch        *b_BigRIPSFocalPlane_nfired_ppacx;   //!
    TBranch        *b_BigRIPSFocalPlane_nfired_ppacy;   //!
    TBranch        *b_BigRIPSFocalPlane_zpos;   //!
    TBranch        *b_BigRIPSFocalPlane_zpos_offset;   //!
    TBranch        *b_BigRIPSRIPS_;   //!
    TBranch        *b_BigRIPSRIPS_fUniqueID;   //!
    TBranch        *b_BigRIPSRIPS_fBits;   //!
    TBranch        *b_BigRIPSRIPS_id;   //!
    TBranch        *b_BigRIPSRIPS_fpl;   //!
    TBranch        *b_BigRIPSRIPS_name;   //!
    TBranch        *b_BigRIPSRIPS_fDataState;   //!
    TBranch        *b_BigRIPSRIPS_upstream_fpl;   //!
    TBranch        *b_BigRIPSRIPS_downstream_fpl;   //!
    TBranch        *b_BigRIPSRIPS_center_brho;   //!
    TBranch        *b_BigRIPSRIPS_brho;   //!
    TBranch        *b_BigRIPSRIPS_length;   //!
    TBranch        *b_BigRIPSRIPS_matrix;   //!
    TBranch        *b_BigRIPSRIPS_delta;   //!
    TBranch        *b_BigRIPSRIPS_angle;   //!
    TBranch        *b_BigRIPSRIPS_dipolename;   //!
    TBranch        *b_BigRIPSTOF_;   //!
    TBranch        *b_BigRIPSTOF_fUniqueID;   //!
    TBranch        *b_BigRIPSTOF_fBits;   //!
    TBranch        *b_BigRIPSTOF_id;   //!
    TBranch        *b_BigRIPSTOF_fpl;   //!
    TBranch        *b_BigRIPSTOF_name;   //!
    TBranch        *b_BigRIPSTOF_fDataState;   //!
    TBranch        *b_BigRIPSTOF_tof;   //!
    TBranch        *b_BigRIPSTOF_clight;   //!
    TBranch        *b_BigRIPSTOF_length;   //!
    TBranch        *b_BigRIPSTOF_ulength;   //!
    TBranch        *b_BigRIPSTOF_dlength;   //!
    TBranch        *b_BigRIPSTOF_upstream_plname;   //!
    TBranch        *b_BigRIPSTOF_downstream_plname;   //!
    TBranch        *b_BigRIPSTOF_upstream_plfpl;   //!
    TBranch        *b_BigRIPSTOF_downstream_plfpl;   //!
    TBranch        *b_BigRIPSTOF_time_offset;   //!
    TBranch        *b_BigRIPSBeam_;   //!
    TBranch        *b_BigRIPSBeam_fUniqueID;   //!
    TBranch        *b_BigRIPSBeam_fBits;   //!
    TBranch        *b_BigRIPSBeam_id;   //!
    TBranch        *b_BigRIPSBeam_fpl;   //!
    TBranch        *b_BigRIPSBeam_name;   //!
    TBranch        *b_BigRIPSBeam_fDataState;   //!
    TBranch        *b_BigRIPSBeam_brho;   //!
    TBranch        *b_BigRIPSBeam_aoq;   //!
    TBranch        *b_BigRIPSBeam_zet;   //!
    TBranch        *b_BigRIPSBeam_beta;   //!
    TBranch        *b_BigRIPSBeam_clight;   //!
    TBranch        *b_BigRIPSBeam_mnucleon;   //!
    TBranch        *b_BigRIPSBeam_nrips;   //!
    TBranch        *b_BigRIPSBeam_ripsname;   //!
    TBranch        *b_BigRIPSBeam_tofname;   //!
    TBranch        *b_BigRIPSBeam_icname;   //!
    TBranch        *b_MinosClustX;   //!
    TBranch        *b_MinosClustY;   //!
    TBranch        *b_MinosClustZ;   //!
    TBranch        *b_MinosClustT;   //!
    TBranch        *b_MinosClustQ;   //!
    TBranch        *b_MinosClustNCl;   //!
    TBranch        *b_MinosClustNP;   //!
    TBranch        *b_MinosClustChi2;   //!
    TBranch        *b_MinosResultX;   //!
    TBranch        *b_MinosResultY;   //!
    TBranch        *b_MinosResultZ;   //!
    TBranch        *b_MinosResultQ;   //!
    TBranch        *b_MinosResultNCl;   //!
    TBranch        *b_MinosResultNP;   //!
    TBranch        *b_MinosResultChi2;   //!
    TBranch        *b_trackNbr;   //!
    TBranch        *b_trackNbr_final;   //!
    TBranch        *b_evtOrig;   //!
    TBranch        *b_MINOS_VDrift;   //!
    TBranch        *b_DelayTrigger;   //!
    TBranch        *b_x_vertex;   //!
    TBranch        *b_y_vertex;   //!
    TBranch        *b_z_vertex;   //!
    TBranch        *b_r_vertex;   //!
    TBranch        *b_phi_vertex;   //!
    TBranch        *b_thetaz1;   //!
    TBranch        *b_thetaz2;   //!
    TBranch        *b_parFit_1;   //!
    TBranch        *b_parFit_2;   //!
    TBranch        *b_F3X;   //!
    TBranch        *b_F3A;   //!
    TBranch        *b_F5X;   //!
    TBranch        *b_F5A;   //!
    TBranch        *b_F7X;   //!
    TBranch        *b_F7A;   //!
    TBranch        *b_F8X;   //!
    TBranch        *b_F8A;   //!
    TBranch        *b_F9X;   //!
    TBranch        *b_F9A;   //!
    TBranch        *b_F11X;   //!
    TBranch        *b_F11A;   //!
    TBranch        *b_fGoodPPACFocus;   //!
    TBranch        *b_fGoodPPACFocusOr;   //!

   treereader(TTree *tree=0);
   virtual ~treereader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual bool     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

    bool singleloop();
    Long64_t NumEntries();
    void setloopkeys(std::vector<std::string> Vals);
    bool getevent(int eventno);
};

#endif

#ifdef treereader_cxx
treereader::treereader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("build/output/out.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("build/output/out.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

treereader::~treereader()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile(); - stack objects don't need this
}

Int_t treereader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t treereader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void treereader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

    // Set object pointer
    MinosClustX = 0;
    MinosClustY = 0;
    MinosClustZ = 0;
    MinosClustT = 0;
    MinosClustQ = 0;
    MinosClustNCl = 0;
    MinosClustNP = 0;
    MinosClustChi2 = 0;
    MinosResultX = 0;
    MinosResultY = 0;
    MinosResultZ = 0;
    MinosResultQ = 0;
    MinosResultNCl = 0;
    MinosResultNP = 0;
    MinosResultChi2 = 0;
    parFit_1 = 0;
    parFit_2 = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

    fChain->SetBranchAddress("EventInfo", &EventInfo_, &b_EventInfo_);
    fChain->SetBranchAddress("EventInfo.fUniqueID", EventInfo_fUniqueID, &b_EventInfo_fUniqueID);
    fChain->SetBranchAddress("EventInfo.fBits", EventInfo_fBits, &b_EventInfo_fBits);
    fChain->SetBranchAddress("EventInfo.fName", EventInfo_fName, &b_EventInfo_fName);
    fChain->SetBranchAddress("EventInfo.fTitle", EventInfo_fTitle, &b_EventInfo_fTitle);
    fChain->SetBranchAddress("EventInfo.runnumber", EventInfo_runnumber, &b_EventInfo_runnumber);
    fChain->SetBranchAddress("EventInfo.eventnumber", EventInfo_eventnumber, &b_EventInfo_eventnumber);
    fChain->SetBranchAddress("EventInfo.subsysname", EventInfo_subsysname, &b_EventInfo_subsysname);
    fChain->SetBranchAddress("EventInfo.timestamp", EventInfo_timestamp, &b_EventInfo_timestamp);
    fChain->SetBranchAddress("EventInfo.comp_val", EventInfo_comp_val, &b_EventInfo_comp_val);
    fChain->SetBranchAddress("EventInfo.fBit", EventInfo_fBit, &b_EventInfo_fBit);
    fChain->SetBranchAddress("BigRIPSPPAC", &BigRIPSPPAC_, &b_BigRIPSPPAC_);
    fChain->SetBranchAddress("BigRIPSPPAC.fUniqueID", BigRIPSPPAC_fUniqueID, &b_BigRIPSPPAC_fUniqueID);
    fChain->SetBranchAddress("BigRIPSPPAC.fBits", BigRIPSPPAC_fBits, &b_BigRIPSPPAC_fBits);
    fChain->SetBranchAddress("BigRIPSPPAC.id", BigRIPSPPAC_id, &b_BigRIPSPPAC_id);
    fChain->SetBranchAddress("BigRIPSPPAC.fpl", BigRIPSPPAC_fpl, &b_BigRIPSPPAC_fpl);
    fChain->SetBranchAddress("BigRIPSPPAC.name", BigRIPSPPAC_name, &b_BigRIPSPPAC_name);
    fChain->SetBranchAddress("BigRIPSPPAC.fDataState", BigRIPSPPAC_fDataState, &b_BigRIPSPPAC_fDataState);
    fChain->SetBranchAddress("BigRIPSPPAC.xzpos", BigRIPSPPAC_xzpos, &b_BigRIPSPPAC_xzpos);
    fChain->SetBranchAddress("BigRIPSPPAC.yzpos", BigRIPSPPAC_yzpos, &b_BigRIPSPPAC_yzpos);
    fChain->SetBranchAddress("BigRIPSPPAC.fTX1Raw", BigRIPSPPAC_fTX1Raw, &b_BigRIPSPPAC_fTX1Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fTX2Raw", BigRIPSPPAC_fTX2Raw, &b_BigRIPSPPAC_fTX2Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fTY1Raw", BigRIPSPPAC_fTY1Raw, &b_BigRIPSPPAC_fTY1Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fTY2Raw", BigRIPSPPAC_fTY2Raw, &b_BigRIPSPPAC_fTY2Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fTARaw", BigRIPSPPAC_fTARaw, &b_BigRIPSPPAC_fTARaw);
    fChain->SetBranchAddress("BigRIPSPPAC.fQX1Raw", BigRIPSPPAC_fQX1Raw, &b_BigRIPSPPAC_fQX1Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fQX2Raw", BigRIPSPPAC_fQX2Raw, &b_BigRIPSPPAC_fQX2Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fQY1Raw", BigRIPSPPAC_fQY1Raw, &b_BigRIPSPPAC_fQY1Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fQY2Raw", BigRIPSPPAC_fQY2Raw, &b_BigRIPSPPAC_fQY2Raw);
    fChain->SetBranchAddress("BigRIPSPPAC.fQARaw", BigRIPSPPAC_fQARaw, &b_BigRIPSPPAC_fQARaw);
    fChain->SetBranchAddress("BigRIPSPPAC.fTX1", BigRIPSPPAC_fTX1, &b_BigRIPSPPAC_fTX1);
    fChain->SetBranchAddress("BigRIPSPPAC.fTX2", BigRIPSPPAC_fTX2, &b_BigRIPSPPAC_fTX2);
    fChain->SetBranchAddress("BigRIPSPPAC.fTY1", BigRIPSPPAC_fTY1, &b_BigRIPSPPAC_fTY1);
    fChain->SetBranchAddress("BigRIPSPPAC.fTY2", BigRIPSPPAC_fTY2, &b_BigRIPSPPAC_fTY2);
    fChain->SetBranchAddress("BigRIPSPPAC.fTA", BigRIPSPPAC_fTA, &b_BigRIPSPPAC_fTA);
    fChain->SetBranchAddress("BigRIPSPPAC.fTSumX", BigRIPSPPAC_fTSumX, &b_BigRIPSPPAC_fTSumX);
    fChain->SetBranchAddress("BigRIPSPPAC.fTSumY", BigRIPSPPAC_fTSumY, &b_BigRIPSPPAC_fTSumY);
    fChain->SetBranchAddress("BigRIPSPPAC.fTDiffX", BigRIPSPPAC_fTDiffX, &b_BigRIPSPPAC_fTDiffX);
    fChain->SetBranchAddress("BigRIPSPPAC.fTDiffY", BigRIPSPPAC_fTDiffY, &b_BigRIPSPPAC_fTDiffY);
    fChain->SetBranchAddress("BigRIPSPPAC.fX", BigRIPSPPAC_fX, &b_BigRIPSPPAC_fX);
    fChain->SetBranchAddress("BigRIPSPPAC.fY", BigRIPSPPAC_fY, &b_BigRIPSPPAC_fY);
    fChain->SetBranchAddress("BigRIPSPPAC.fFiredX", BigRIPSPPAC_fFiredX, &b_BigRIPSPPAC_fFiredX);
    fChain->SetBranchAddress("BigRIPSPPAC.fFiredY", BigRIPSPPAC_fFiredY, &b_BigRIPSPPAC_fFiredY);
    fChain->SetBranchAddress("BigRIPSPlastic", &BigRIPSPlastic_, &b_BigRIPSPlastic_);
    fChain->SetBranchAddress("BigRIPSPlastic.fUniqueID", BigRIPSPlastic_fUniqueID, &b_BigRIPSPlastic_fUniqueID);
    fChain->SetBranchAddress("BigRIPSPlastic.fBits", BigRIPSPlastic_fBits, &b_BigRIPSPlastic_fBits);
    fChain->SetBranchAddress("BigRIPSPlastic.id", BigRIPSPlastic_id, &b_BigRIPSPlastic_id);
    fChain->SetBranchAddress("BigRIPSPlastic.fpl", BigRIPSPlastic_fpl, &b_BigRIPSPlastic_fpl);
    fChain->SetBranchAddress("BigRIPSPlastic.name", BigRIPSPlastic_name, &b_BigRIPSPlastic_name);
    fChain->SetBranchAddress("BigRIPSPlastic.fDataState", BigRIPSPlastic_fDataState, &b_BigRIPSPlastic_fDataState);
    fChain->SetBranchAddress("BigRIPSPlastic.zposition", BigRIPSPlastic_zposition, &b_BigRIPSPlastic_zposition);
    fChain->SetBranchAddress("BigRIPSPlastic.zoffset", BigRIPSPlastic_zoffset, &b_BigRIPSPlastic_zoffset);
    fChain->SetBranchAddress("BigRIPSPlastic.tcalL", BigRIPSPlastic_tcalL, &b_BigRIPSPlastic_tcalL);
    fChain->SetBranchAddress("BigRIPSPlastic.tcalR", BigRIPSPlastic_tcalR, &b_BigRIPSPlastic_tcalR);
    fChain->SetBranchAddress("BigRIPSPlastic.fTLRaw", BigRIPSPlastic_fTLRaw, &b_BigRIPSPlastic_fTLRaw);
    fChain->SetBranchAddress("BigRIPSPlastic.fTRRaw", BigRIPSPlastic_fTRRaw, &b_BigRIPSPlastic_fTRRaw);
    fChain->SetBranchAddress("BigRIPSPlastic.fQLRaw", BigRIPSPlastic_fQLRaw, &b_BigRIPSPlastic_fQLRaw);
    fChain->SetBranchAddress("BigRIPSPlastic.fQRRaw", BigRIPSPlastic_fQRRaw, &b_BigRIPSPlastic_fQRRaw);
    fChain->SetBranchAddress("BigRIPSPlastic.fQTCLRawWidth", BigRIPSPlastic_fQTCLRawWidth, &b_BigRIPSPlastic_fQTCLRawWidth);
    fChain->SetBranchAddress("BigRIPSPlastic.fQTCRRawWidth", BigRIPSPlastic_fQTCRRawWidth, &b_BigRIPSPlastic_fQTCRRawWidth);
    fChain->SetBranchAddress("BigRIPSPlastic.fQTCLRawStart", BigRIPSPlastic_fQTCLRawStart, &b_BigRIPSPlastic_fQTCLRawStart);
    fChain->SetBranchAddress("BigRIPSPlastic.fQTCRRawStart", BigRIPSPlastic_fQTCRRawStart, &b_BigRIPSPlastic_fQTCRRawStart);
    fChain->SetBranchAddress("BigRIPSPlastic.fTLRawArray", BigRIPSPlastic_fTLRawArray, &b_BigRIPSPlastic_fTLRawArray);
    fChain->SetBranchAddress("BigRIPSPlastic.fTRRawArray", BigRIPSPlastic_fTRRawArray, &b_BigRIPSPlastic_fTRRawArray);
    fChain->SetBranchAddress("BigRIPSPlastic.fTime", BigRIPSPlastic_fTime, &b_BigRIPSPlastic_fTime);
    fChain->SetBranchAddress("BigRIPSPlastic.fTimeL", BigRIPSPlastic_fTimeL, &b_BigRIPSPlastic_fTimeL);
    fChain->SetBranchAddress("BigRIPSPlastic.fTimeR", BigRIPSPlastic_fTimeR, &b_BigRIPSPlastic_fTimeR);
    fChain->SetBranchAddress("BigRIPSPlastic.fTimeLSlew", BigRIPSPlastic_fTimeLSlew, &b_BigRIPSPlastic_fTimeLSlew);
    fChain->SetBranchAddress("BigRIPSPlastic.fTimeRSlew", BigRIPSPlastic_fTimeRSlew, &b_BigRIPSPlastic_fTimeRSlew);
    fChain->SetBranchAddress("BigRIPSPlastic.fTimeSlew", BigRIPSPlastic_fTimeSlew, &b_BigRIPSPlastic_fTimeSlew);
    fChain->SetBranchAddress("BigRIPSIC", &BigRIPSIC_, &b_BigRIPSIC_);
    fChain->SetBranchAddress("BigRIPSIC.fUniqueID", BigRIPSIC_fUniqueID, &b_BigRIPSIC_fUniqueID);
    fChain->SetBranchAddress("BigRIPSIC.fBits", BigRIPSIC_fBits, &b_BigRIPSIC_fBits);
    fChain->SetBranchAddress("BigRIPSIC.id", BigRIPSIC_id, &b_BigRIPSIC_id);
    fChain->SetBranchAddress("BigRIPSIC.fpl", BigRIPSIC_fpl, &b_BigRIPSIC_fpl);
    fChain->SetBranchAddress("BigRIPSIC.name", BigRIPSIC_name, &b_BigRIPSIC_name);
    fChain->SetBranchAddress("BigRIPSIC.fDataState", BigRIPSIC_fDataState, &b_BigRIPSIC_fDataState);
    fChain->SetBranchAddress("BigRIPSIC.zcoef[2]", BigRIPSIC_zcoef, &b_BigRIPSIC_zcoef);
    fChain->SetBranchAddress("BigRIPSIC.ionpair", BigRIPSIC_ionpair, &b_BigRIPSIC_ionpair);
    fChain->SetBranchAddress("BigRIPSIC.nhitchannel", BigRIPSIC_nhitchannel, &b_BigRIPSIC_nhitchannel);
    fChain->SetBranchAddress("BigRIPSIC.fADC[32]", BigRIPSIC_fADC, &b_BigRIPSIC_fADC);
    fChain->SetBranchAddress("BigRIPSIC.fEnergy[32]", BigRIPSIC_fEnergy, &b_BigRIPSIC_fEnergy);
    fChain->SetBranchAddress("BigRIPSIC.fTDC[32]", BigRIPSIC_fTDC, &b_BigRIPSIC_fTDC);
    fChain->SetBranchAddress("BigRIPSIC.fTime[32]", BigRIPSIC_fTime, &b_BigRIPSIC_fTime);
    fChain->SetBranchAddress("BigRIPSIC.fRawADCSqSum", BigRIPSIC_fRawADCSqSum, &b_BigRIPSIC_fRawADCSqSum);
    fChain->SetBranchAddress("BigRIPSIC.fRawADCAvSum", BigRIPSIC_fRawADCAvSum, &b_BigRIPSIC_fRawADCAvSum);
    fChain->SetBranchAddress("BigRIPSIC.fCalMeVSqSum", BigRIPSIC_fCalMeVSqSum, &b_BigRIPSIC_fCalMeVSqSum);
    fChain->SetBranchAddress("BigRIPSIC.fCalMeVAvSum", BigRIPSIC_fCalMeVAvSum, &b_BigRIPSIC_fCalMeVAvSum);
    fChain->SetBranchAddress("BigRIPSFocalPlane", &BigRIPSFocalPlane_, &b_BigRIPSFocalPlane_);
    fChain->SetBranchAddress("BigRIPSFocalPlane.fUniqueID", BigRIPSFocalPlane_fUniqueID, &b_BigRIPSFocalPlane_fUniqueID);
    fChain->SetBranchAddress("BigRIPSFocalPlane.fBits", BigRIPSFocalPlane_fBits, &b_BigRIPSFocalPlane_fBits);
    fChain->SetBranchAddress("BigRIPSFocalPlane.id", BigRIPSFocalPlane_id, &b_BigRIPSFocalPlane_id);
    fChain->SetBranchAddress("BigRIPSFocalPlane.fpl", BigRIPSFocalPlane_fpl, &b_BigRIPSFocalPlane_fpl);
    fChain->SetBranchAddress("BigRIPSFocalPlane.name", BigRIPSFocalPlane_name, &b_BigRIPSFocalPlane_name);
    fChain->SetBranchAddress("BigRIPSFocalPlane.fDataState", BigRIPSFocalPlane_fDataState, &b_BigRIPSFocalPlane_fDataState);
    fChain->SetBranchAddress("BigRIPSFocalPlane.opt_vector", BigRIPSFocalPlane_opt_vector, &b_BigRIPSFocalPlane_opt_vector);
    fChain->SetBranchAddress("BigRIPSFocalPlane.X", BigRIPSFocalPlane_X, &b_BigRIPSFocalPlane_X);
    fChain->SetBranchAddress("BigRIPSFocalPlane.A", BigRIPSFocalPlane_A, &b_BigRIPSFocalPlane_A);
    fChain->SetBranchAddress("BigRIPSFocalPlane.Y", BigRIPSFocalPlane_Y, &b_BigRIPSFocalPlane_Y);
    fChain->SetBranchAddress("BigRIPSFocalPlane.B", BigRIPSFocalPlane_B, &b_BigRIPSFocalPlane_B);
    fChain->SetBranchAddress("BigRIPSFocalPlane.nfired_ppacx", BigRIPSFocalPlane_nfired_ppacx, &b_BigRIPSFocalPlane_nfired_ppacx);
    fChain->SetBranchAddress("BigRIPSFocalPlane.nfired_ppacy", BigRIPSFocalPlane_nfired_ppacy, &b_BigRIPSFocalPlane_nfired_ppacy);
    fChain->SetBranchAddress("BigRIPSFocalPlane.zpos", BigRIPSFocalPlane_zpos, &b_BigRIPSFocalPlane_zpos);
    fChain->SetBranchAddress("BigRIPSFocalPlane.zpos_offset", BigRIPSFocalPlane_zpos_offset, &b_BigRIPSFocalPlane_zpos_offset);
    fChain->SetBranchAddress("BigRIPSRIPS", &BigRIPSRIPS_, &b_BigRIPSRIPS_);
    fChain->SetBranchAddress("BigRIPSRIPS.fUniqueID", BigRIPSRIPS_fUniqueID, &b_BigRIPSRIPS_fUniqueID);
    fChain->SetBranchAddress("BigRIPSRIPS.fBits", BigRIPSRIPS_fBits, &b_BigRIPSRIPS_fBits);
    fChain->SetBranchAddress("BigRIPSRIPS.id", BigRIPSRIPS_id, &b_BigRIPSRIPS_id);
    fChain->SetBranchAddress("BigRIPSRIPS.fpl", BigRIPSRIPS_fpl, &b_BigRIPSRIPS_fpl);
    fChain->SetBranchAddress("BigRIPSRIPS.name", BigRIPSRIPS_name, &b_BigRIPSRIPS_name);
    fChain->SetBranchAddress("BigRIPSRIPS.fDataState", BigRIPSRIPS_fDataState, &b_BigRIPSRIPS_fDataState);
    fChain->SetBranchAddress("BigRIPSRIPS.upstream_fpl", BigRIPSRIPS_upstream_fpl, &b_BigRIPSRIPS_upstream_fpl);
    fChain->SetBranchAddress("BigRIPSRIPS.downstream_fpl", BigRIPSRIPS_downstream_fpl, &b_BigRIPSRIPS_downstream_fpl);
    fChain->SetBranchAddress("BigRIPSRIPS.center_brho", BigRIPSRIPS_center_brho, &b_BigRIPSRIPS_center_brho);
    fChain->SetBranchAddress("BigRIPSRIPS.brho", BigRIPSRIPS_brho, &b_BigRIPSRIPS_brho);
    fChain->SetBranchAddress("BigRIPSRIPS.length", BigRIPSRIPS_length, &b_BigRIPSRIPS_length);
    fChain->SetBranchAddress("BigRIPSRIPS.matrix", BigRIPSRIPS_matrix, &b_BigRIPSRIPS_matrix);
    fChain->SetBranchAddress("BigRIPSRIPS.delta", BigRIPSRIPS_delta, &b_BigRIPSRIPS_delta);
    fChain->SetBranchAddress("BigRIPSRIPS.angle", BigRIPSRIPS_angle, &b_BigRIPSRIPS_angle);
    fChain->SetBranchAddress("BigRIPSRIPS.dipolename", BigRIPSRIPS_dipolename, &b_BigRIPSRIPS_dipolename);
    fChain->SetBranchAddress("BigRIPSTOF", &BigRIPSTOF_, &b_BigRIPSTOF_);
    fChain->SetBranchAddress("BigRIPSTOF.fUniqueID", BigRIPSTOF_fUniqueID, &b_BigRIPSTOF_fUniqueID);
    fChain->SetBranchAddress("BigRIPSTOF.fBits", BigRIPSTOF_fBits, &b_BigRIPSTOF_fBits);
    fChain->SetBranchAddress("BigRIPSTOF.id", BigRIPSTOF_id, &b_BigRIPSTOF_id);
    fChain->SetBranchAddress("BigRIPSTOF.fpl", BigRIPSTOF_fpl, &b_BigRIPSTOF_fpl);
    fChain->SetBranchAddress("BigRIPSTOF.name", BigRIPSTOF_name, &b_BigRIPSTOF_name);
    fChain->SetBranchAddress("BigRIPSTOF.fDataState", BigRIPSTOF_fDataState, &b_BigRIPSTOF_fDataState);
    fChain->SetBranchAddress("BigRIPSTOF.tof", BigRIPSTOF_tof, &b_BigRIPSTOF_tof);
    fChain->SetBranchAddress("BigRIPSTOF.clight", BigRIPSTOF_clight, &b_BigRIPSTOF_clight);
    fChain->SetBranchAddress("BigRIPSTOF.length", BigRIPSTOF_length, &b_BigRIPSTOF_length);
    fChain->SetBranchAddress("BigRIPSTOF.ulength", BigRIPSTOF_ulength, &b_BigRIPSTOF_ulength);
    fChain->SetBranchAddress("BigRIPSTOF.dlength", BigRIPSTOF_dlength, &b_BigRIPSTOF_dlength);
    fChain->SetBranchAddress("BigRIPSTOF.upstream_plname", BigRIPSTOF_upstream_plname, &b_BigRIPSTOF_upstream_plname);
    fChain->SetBranchAddress("BigRIPSTOF.downstream_plname", BigRIPSTOF_downstream_plname, &b_BigRIPSTOF_downstream_plname);
    fChain->SetBranchAddress("BigRIPSTOF.upstream_plfpl", BigRIPSTOF_upstream_plfpl, &b_BigRIPSTOF_upstream_plfpl);
    fChain->SetBranchAddress("BigRIPSTOF.downstream_plfpl", BigRIPSTOF_downstream_plfpl, &b_BigRIPSTOF_downstream_plfpl);
    fChain->SetBranchAddress("BigRIPSTOF.time_offset", BigRIPSTOF_time_offset, &b_BigRIPSTOF_time_offset);
    fChain->SetBranchAddress("BigRIPSBeam", &BigRIPSBeam_, &b_BigRIPSBeam_);
    fChain->SetBranchAddress("BigRIPSBeam.fUniqueID", BigRIPSBeam_fUniqueID, &b_BigRIPSBeam_fUniqueID);
    fChain->SetBranchAddress("BigRIPSBeam.fBits", BigRIPSBeam_fBits, &b_BigRIPSBeam_fBits);
    fChain->SetBranchAddress("BigRIPSBeam.id", BigRIPSBeam_id, &b_BigRIPSBeam_id);
    fChain->SetBranchAddress("BigRIPSBeam.fpl", BigRIPSBeam_fpl, &b_BigRIPSBeam_fpl);
    fChain->SetBranchAddress("BigRIPSBeam.name", BigRIPSBeam_name, &b_BigRIPSBeam_name);
    fChain->SetBranchAddress("BigRIPSBeam.fDataState", BigRIPSBeam_fDataState, &b_BigRIPSBeam_fDataState);
    fChain->SetBranchAddress("BigRIPSBeam.brho", BigRIPSBeam_brho, &b_BigRIPSBeam_brho);
    fChain->SetBranchAddress("BigRIPSBeam.aoq", BigRIPSBeam_aoq, &b_BigRIPSBeam_aoq);
    fChain->SetBranchAddress("BigRIPSBeam.zet", BigRIPSBeam_zet, &b_BigRIPSBeam_zet);
    fChain->SetBranchAddress("BigRIPSBeam.beta", BigRIPSBeam_beta, &b_BigRIPSBeam_beta);
    fChain->SetBranchAddress("BigRIPSBeam.clight", BigRIPSBeam_clight, &b_BigRIPSBeam_clight);
    fChain->SetBranchAddress("BigRIPSBeam.mnucleon", BigRIPSBeam_mnucleon, &b_BigRIPSBeam_mnucleon);
    fChain->SetBranchAddress("BigRIPSBeam.nrips", BigRIPSBeam_nrips, &b_BigRIPSBeam_nrips);
    fChain->SetBranchAddress("BigRIPSBeam.ripsname[2]", BigRIPSBeam_ripsname, &b_BigRIPSBeam_ripsname);
    fChain->SetBranchAddress("BigRIPSBeam.tofname", BigRIPSBeam_tofname, &b_BigRIPSBeam_tofname);
    fChain->SetBranchAddress("BigRIPSBeam.icname", BigRIPSBeam_icname, &b_BigRIPSBeam_icname);
    fChain->SetBranchAddress("MinosClustX", &MinosClustX, &b_MinosClustX);
    fChain->SetBranchAddress("MinosClustY", &MinosClustY, &b_MinosClustY);
    fChain->SetBranchAddress("MinosClustZ", &MinosClustZ, &b_MinosClustZ);
    fChain->SetBranchAddress("MinosClustT", &MinosClustT, &b_MinosClustT);
    fChain->SetBranchAddress("MinosClustQ", &MinosClustQ, &b_MinosClustQ);
    fChain->SetBranchAddress("MinosClustNCl", &MinosClustNCl, &b_MinosClustNCl);
    fChain->SetBranchAddress("MinosClustNP", &MinosClustNP, &b_MinosClustNP);
    fChain->SetBranchAddress("MinosClustChi2", &MinosClustChi2, &b_MinosClustChi2);
    fChain->SetBranchAddress("MinosResultX", &MinosResultX, &b_MinosResultX);
    fChain->SetBranchAddress("MinosResultY", &MinosResultY, &b_MinosResultY);
    fChain->SetBranchAddress("MinosResultZ", &MinosResultZ, &b_MinosResultZ);
    fChain->SetBranchAddress("MinosResultQ", &MinosResultQ, &b_MinosResultQ);
    fChain->SetBranchAddress("MinosResultNCl", &MinosResultNCl, &b_MinosResultNCl);
    fChain->SetBranchAddress("MinosResultNP", &MinosResultNP, &b_MinosResultNP);
    fChain->SetBranchAddress("MinosResultChi2", &MinosResultChi2, &b_MinosResultChi2);
    fChain->SetBranchAddress("trackNbr", &trackNbr, &b_trackNbr);
    fChain->SetBranchAddress("trackNbr_final", &trackNbr_final, &b_trackNbr_final);
    fChain->SetBranchAddress("evtOrig", &evtOrig, &b_evtOrig);
    fChain->SetBranchAddress("VDrift", &VDrift, &b_MINOS_VDrift);
    fChain->SetBranchAddress("DelayTrig", &DelayTrig, &b_DelayTrigger);
    fChain->SetBranchAddress("x_vertex", &x_vertex, &b_x_vertex);
    fChain->SetBranchAddress("y_vertex", &y_vertex, &b_y_vertex);
    fChain->SetBranchAddress("z_vertex", &z_vertex, &b_z_vertex);
    fChain->SetBranchAddress("r_vertex", &r_vertex, &b_r_vertex);
    fChain->SetBranchAddress("phi_vertex", &phi_vertex, &b_phi_vertex);
    fChain->SetBranchAddress("thetaz1", &thetaz1, &b_thetaz1);
    fChain->SetBranchAddress("thetaz2", &thetaz2, &b_thetaz2);
    fChain->SetBranchAddress("parFit_1", &parFit_1, &b_parFit_1);
    fChain->SetBranchAddress("parFit_2", &parFit_2, &b_parFit_2);
    fChain->SetBranchAddress("F3X", &F3X, &b_F3X);
    fChain->SetBranchAddress("F3A", &F3A, &b_F3A);
    fChain->SetBranchAddress("F5X", &F5X, &b_F5X);
    fChain->SetBranchAddress("F5A", &F5A, &b_F5A);
    fChain->SetBranchAddress("F7X", &F7X, &b_F7X);
    fChain->SetBranchAddress("F7A", &F7A, &b_F7A);
    fChain->SetBranchAddress("F8X", &F8X, &b_F8X);
    fChain->SetBranchAddress("F8A", &F8A, &b_F8A);
    fChain->SetBranchAddress("F9X", &F9X, &b_F9X);
    fChain->SetBranchAddress("F9A", &F9A, &b_F9A);
    fChain->SetBranchAddress("F11X", &F11X, &b_F11X);
    fChain->SetBranchAddress("F11A", &F11A, &b_F11A);
    fChain->SetBranchAddress("fgoodppacfocus", fgoodppacfocus, &b_fGoodPPACFocus);
    fChain->SetBranchAddress("fgoodppacfocusor", fgoodppacfocusor, &b_fGoodPPACFocusOr);
   Notify();
}

Bool_t treereader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void treereader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t treereader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef treereader_cxx
