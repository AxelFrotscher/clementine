#define treereader_cxx
#include "treereader.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

bool treereader::Loop()
{
//   In a ROOT session, you can do:
//      root> .L treereader.C
//      root> treereader t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
    fChain->SetBranchStatus("*",0);  // disable all branches
    fChain->SetBranchStatus("BigRIPSIC.fADC[32]",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return false;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) return false;
      nb = fChain->GetEntry(jentry);
      std::cout << BigRIPSIC_fADC[0][4] << std::endl;
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
   return true;
}

bool treereader::singleloop(){
    // This method aims at providing a neat interface for gettings data values
    Long64_t size = LoadTree(currenteventnum);
    if(size < 0 ) {
        return false;
    }
    fChain->GetEntry(currenteventnum);
    currenteventnum++;
    return true;
}

bool treereader::getevent(int eventno){
    Long64_t size = LoadTree(eventno);
    if(size<0) return false;
    fChain->GetEntry(eventno);

    return true;
}

void treereader::setloopkeys(std::vector <std::string> Vals){
    // Set Branches we want to read
    currenteventnum = 0; // reset event number on change of readout
    if(!fChain) __throw_exception_again;

    fChain->SetBranchStatus("*",false);
    for(auto name : Vals)
        fChain->SetBranchStatus(name.c_str(),true);

    Long64_t nentries = fChain->GetEntriesFast();
}

Long64_t treereader::NumEntries(){
    return fChain->GetEntries();
}