#include "histograms.hh"
#include "libconstant.h"
#include "chargestatecut.h"
#include "ICcut.h"
#include "plasticcut.h"
#include "ppaccut.h"
#include "targetcut.h"
#include "triggercut.h"
#include "higherorder.h"
#include "PID/pid.h"
#include "txtwriter.h"
#include <numeric>
#include "zdssetting.h"
#include "TGraphErrors.h"
#include <filesystem>

using std::cout, std::endl, std::string, std::vector, std::atomic, std::to_string,
      std::__throw_invalid_argument;

 void makehistograms(vector<string> input) {
    /// Generating an outputfile that matches names with the input file
    /// Determine run type:
    TChain chain("tree");
    for(auto &i:input) chain.Add(i.c_str());

    //setting set;
    setting::setcountno((int)chain.GetEntries());
    const string settingname = setting::getsetname();
    const string modename = setting::getmodename();
    //const bool minosbool = setting::getminos();

    /// lambda to generate expression
    auto gentxt = [input, settingname, modename](auto suffix, const bool plastic,
            const bool cs){
        if(!plastic){
            //Part for generating file path for plastic analysis root-file (w/o the tree)
            if(setting::getminos() && !cs){
                return "/home/afrotscher/Clementine/build/output/"+ settingname +
                       "/" + modename + "_MINOS_" + input.at(0).substr(34,9) +
                       suffix;
            }
            else if(!cs) return "/home/afrotscher/Clementine/build/output/"+ settingname +
                        "/" + modename + "_CS_" + input.at(0).substr(34,9) + suffix;
            //Part for the cross sections
            else{
                if(setting::getminos()){
                    return "/home/afrotscher/Clementine/build/output/"+ settingname +
                           "/" + modename + "_MINOS" + suffix;
                }
                else return "/home/afrotscher/Clementine/build/output/"+ settingname +
                            "/" + modename + "_CS" + suffix;
            }
        }
        else{
            if(setting::getminos()){
                return "/d/d02-1/ag_ob/afrotscher/SEASTAR/MINOS/" +
                       settingname + suffix;
            }
            else return "/d/d02-1/ag_ob/afrotscher/SEASTAR/CS/" +
                        settingname + suffix;
        }
    };

    /// Initialize ROOT and txt outputfile
    txtwriter writetotxt(gentxt(".txt", false, true)); // Writer class
    
    bool skip_cuts = false;
    vector<string> preanalysed_files;
    if(std::filesystem::exists(gentxt(".root", true, false))){
        preanalysed_files.push_back(gentxt(".root", true, false));
        //Find out if other files exist
        int i =1;
        while(std::filesystem::exists(
                gentxt("_"+ to_string(i) +".root", true, false))){
            preanalysed_files.push_back(gentxt("_"+ to_string(i) +".root", true, false));
            i++;
        }
        cout << preanalysed_files.size() << " Preanalysed file(s): " << endl;
        for(auto &j:preanalysed_files) cout << j << endl;
        cout << "already exists. Skip cut process? yes [1], no [0]" << endl;
        if(!(std::cin >> skip_cuts)) __throw_invalid_argument("Not 0 or 1!");
    }
    
    if(!skip_cuts) {
        cout << "Beginning reconstruction of " << chain.GetEntries()
             << " Elements." << endl;
        
        preanalysed_files.clear(); //delete old entries, as we want new ones
        auto outputfile = new TFile(gentxt(".root", false, false).c_str(),
                                    "RECREATE");
        assert(outputfile->IsOpen());
    
        vector<vector<atomic<bool>>> goodevents;
        goodevents.reserve(chain.GetEntries());
        for(uint i=0; i<(uint)chain.GetEntries(); i++){
            goodevents.emplace_back(vector<atomic<bool>>(2));
        }
        for(auto &i:goodevents) for(auto &j:i) j.exchange(true);
        
        cout << "Making cuts..." << endl;
        triggercut(input, goodevents);
        ccsc(input, goodevents, outputfile);
        targetcut(input, goodevents, outputfile);
        ppaccut(input, goodevents, outputfile);
        plasticcut(input, goodevents, outputfile);
        iccut(input, goodevents, outputfile);
        higherorder(input, goodevents, outputfile);
    
        // Write out file
        writeroot(input, gentxt(".root", true, false), goodevents, setting::getminos());
        outputfile->Close();
        //Push back newly analysed file to the vector.
        preanalysed_files.push_back(gentxt(".root", true, false));
    
        //Check if there have been unwanted file splits, and treat them as well
        int i =1;
        while(std::filesystem::exists(
                gentxt("_"+ to_string(i) +".root", true, false))){
            preanalysed_files.push_back(gentxt("_"+ to_string(i) +".root", true, false));
            i++;
        }
    }

    //Set new input on file we just generated
    input = {preanalysed_files};
    treereader::PID_only = true; // Load a limited selection of branches
    treereader::additional_minos = setting::getminos();
    
    auto csfile = new TFile(gentxt(".root", false, true).c_str(), "RECREATE");
    assert(csfile->IsOpen());
    
    /// Get Z vs. A/Q
    const vector<string> reactionmodes = setting::getreactions();

    /// All reaction modes:
    TGraphErrors crosssection;
    crosssection.SetName(Form("cs%s", settingname.c_str()));
    crosssection.SetTitle("Cross Sections Frotscher 2019");

    for(auto &i: reactionmodes){
        PID(input,csfile,i, crosssection);
    }

    TGraphErrors nancytcs = nancycs(setting::getsetnumber());
    csfile->cd();
    nancytcs.Write();
    crosssection.Write();

    printf("Made PID histograms in %s\n", gentxt(".root",false, true).c_str());

    writetotxt.writetofile();
    csfile->Close();
}

TGraphErrors nancycs(const int &setnumber){
     // Create a TGraph with Nancies values
     TGraphErrors temp;
     temp.SetTitle("Cross Sections Hupin 2019");
     temp.SetName("csnancy");

     vector<vector<double>> massnumbers{
         {110.,110.,111.,111.,112.,113.},  //nb,mo,nb,mo,mo,tc
         {},
         {95,95,96,97},  //br,kr,kr,rb
         {100},
         {67,68,68,69,70}, // fe,fe,co,co,ni
         {72,73,74,74,75}, // ni,ni,ni,cu,cu
         {78,79,81}        // zn,zn,ga
     };

     vector<vector<double>> crosssections{
         {3.0,8,4.3,5.9,7.4,6.5},
         {},
         {2.5,6,6.7,4.7},
         {9},
         {7,9.2,5.4,8.0,12},
         {10,7.0,7.3,4.8,6.6},
         {8,5.3,4.7}
     };

     vector<vector<double>> crosssectione{
         {0.4,1,0.9,0.7,0.8,0.7},
         {},
         {0.7,3,0.7,0.6},
         {1},
         {0.3,0.6,0.5,0.5,1},
         {1,0.4,0.4,0.3,0.5},
         {1,0.4,0.3}
     };

     for(unsigned long i=0; i<massnumbers.at(setnumber).size(); i++){
         temp.SetPoint(temp.GetN(),massnumbers.at(setnumber).at(i),
                                   crosssections.at(setnumber).at(i));
         temp.SetPointError(temp.GetN()-1, 0, crosssectione.at(setnumber).at(i));
     }
     temp.SetPoint(temp.GetN(), 50,0);
     temp.SetPoint(temp.GetN(), 120,0);
     return temp;
 }
 
void writeroot(const vector<string> &input, const string &out,
               const vector<vector<atomic<bool>>> &goodevents, const bool minosbool){
    /// Write all elements to a separate root file:
    
    cout << "Writing cut events to file: " << endl
         << out << endl;
    
    auto inter_file = new TFile(out.c_str(), "RECREATE");
    assert(inter_file->IsOpen());
    inter_file->SetCompressionLevel(6);
    inter_file->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kZLIB);
    
    treereader tree(input);
    
    if(minosbool){ //Tree with events considered for MINOS
        tree.setloopkeys_minos();  //Only write branches necessary for further analysis.
        TTree *inter_tree = tree.fChain->CloneTree(0);
        progressbar progress(tree.fChain->GetEntries(), 0);
        for(long long i=0; i<tree.fChain->GetEntries(); i++){
            if(goodevents.at(i).at(0) && goodevents.at(i).at(1)){
                tree.GetEntry(i);
                tree.EventInfo_fUniqueID[0] = 1; //F1-11 event-flag
                inter_tree->Fill();
            }
            progress.increaseevent();
            if(!(i%5000)) progress.draw();
        }
        progressbar::reset();
    }
    else{ // Tree with events considered for cross sections
        tree.setloopkeys_cs();  //Only write branches necessary for further analysis.
        TTree *inter_tree = tree.fChain->CloneTree(0);
        progressbar progress(tree.fChain->GetEntries(), 0);
        for(long long i=0; i<tree.fChain->GetEntries(); i++){
            if(goodevents.at(i).at(0)){
                tree.GetEntry(i);
                if(goodevents.at(i).at(1)) tree.EventInfo_fUniqueID[0] = 1;
                else tree.EventInfo_fUniqueID[0] = 0; // not F1-F11
                inter_tree->Fill();
            }
            progress.increaseevent();
            if(!(i%50000)) progress.draw();
        }
        progressbar::reset();
        cout << "Done with F1-11 events." << endl;
    }
    //inter_file->Write();
    inter_file->Close();
    cout << "... Done!" << endl;
 }