#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility> 
#include <typeinfo>
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "GenAnalysis.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
GenAnalysis::GenAnalysis()
  : AnaBase(),
    _dumpEvent(false)
{
 
}
// ----------
// Destructor
// ----------
GenAnalysis::~GenAnalysis() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool GenAnalysis::beginJob() 
{ 
  AnaBase::beginJob();

  // Open the output ROOT file
  histf()->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void GenAnalysis::bookHistograms() 
{
 new TH1F("genHmass", "Gen Higgs Mass", 200, 0., 200.);
 new TH1F("genHmassfromZ", "Gen Higgs Mass from Z", 200, 0., 200.);
 new TH1F("genHmassfromZdau", "Gen Higgs Mass from Z daughters", 400, 0., 200.);
 new TH1F("genZmother", "Z mother Id", 30, 0, 30);
 new TH1F("genZmass", "Gen Z Mass", 150, 0., 150.);
 new TH1F("genZ1mass", "Gen Z1 Mass", 150, 0., 150.);
 new TH1F("genZ2mass", "Gen Z2 Mass", 150, 0., 150.);
 new TH1F("nZdauflav", "Flavour of daughters Z ", 10, -0.5, 9.5);
 new TH1F("nZ1dauflav", "Flavour of daughters Z1 ", 10, -0.5, 9.5);
 new TH1F("nZ2dauflav", "Flavour of daughters Z2 ", 10, -0.5, 9.5);
 new TH1F("nEvntflav", "Flavour of ZZ decay ", 10, -0.5, 9.5);
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------
void GenAnalysis::clearLists() {
 vtxList.clear();
 genZ.clear();
}
// -------------------
// The main event loop
// -------------------
void GenAnalysis::eventLoop()
{
    // Initialize analysis
    if (!beginJob()) return;
    int nPrint = max(10000, nEvents()/1000);

    Options op;
    op.verbose = false;
    op.usesbit = true;  // Crucial
    op.printselected = false;
  
  
    // --------------------
    // Start the event loop
    // --------------------
    string lastFile;
    std::cout<<"Bunch Crossing>>>>"<<bunchCrossing()<<std::endl;
    int n4mu=0,n2e2mu=0,n4e=0,nZ1clash=0,ntauev=0,nHtoZ=0;
    for (int ev = 0; ev < nEvents(); ++ev) {
        clearEvent();
        clearLists();
        int lflag = chain()->LoadTree(ev);
        int nbytes = getEntry(lflag);    // returns total bytes read
        string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));
        const Event& evt = eventColl()->at(0);

        histf()->cd();

        //For data or for MC without pileup
        puevWt_ = 1;
        /*
         if (isMC()) {
         int npu = 0;
         puevWt_ = wtPileUp(npu);
         }
        */
        // Show status of the run
        int run   = evt.run;
        int event = evt.event;
        int lumis = evt.lumis;
        
        // Show status of the run
        if (currentFile != lastFile)
            cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
            << " ==> " << currentFile
            << " <<< Run# " << run
            << " Lumis# " << lumis
            << " Event# " << setw(8) << event << " >>> "
            << " Events proc. " << setw(8) << ev
            << endl;
        lastFile = currentFile;

        // Show the status
        if (ev%nPrint == 0)
            cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
            << " ==> " << chain()->GetCurrentFile()->GetName()
            << " <<< Run# " << run
            << " Lumis# " << lumis
            << " Event# " << setw(8) << event << " >>> "
            << " Events proc. " << setw(8) << ev
            << endl;
        //good vertex finding
        op.verbose = (logOption() >> 1 & 0x1);
        findVtxInfo(vtxList, op, fLog());
        double ngoodVtx=vtxList.size();
        //Gen Analysis study
        for (auto jt = genParticleColl()->begin(); jt != genParticleColl()->end(); ++jt) {
          const GenParticle& gp = (*jt);
          int pdgid = fabs(gp.pdgId);
          int status = gp.status;
          if( pdgid == 25 ) {
            int mmid = -1;
            int index = getMotherId(gp, mmid);
            //cstd::cout<<"Higgs status>>"<< status << "  Mother>>" << index <<std::endl;
            TLorentzVector l;
            l.SetPtEtaPhiE(gp.pt,gp.eta,gp.phi,gp.energy);
            AnaUtil::fillHist1D("genHmass",l.M(),1);
          }
          if (pdgid == 23 ) {
            int mmid = -1;
            int index = getMotherId(gp, mmid);
            AnaUtil::fillHist1D("genZmother",mmid,1);
            if (index < 0 || fabs(mmid) != 25) continue; // H decay to Higgs
            nHtoZ++;
            genZ.push_back(gp);
            TLorentzVector l1;
            l1.SetPtEtaPhiE(gp.pt,gp.eta,gp.phi,gp.energy);
            AnaUtil::fillHist1D("genZmass",l1.M(),1);
           //check e mu decay of Z
           vector<int> d = gp.daughtIndices;
           int ndau = 0;
           for (auto it = d.begin(); it != d.end(); ++it) {
             int di = (*it);
             if (di >= ngenparticle()) continue;
             const GenParticle& dgp = genParticleColl()->at(di);
             int pid = fabs(dgp.pdgId);
             if (pid == 13) AnaUtil::fillHist1D("nZdauflav",1,1);
             else if (pid == 11) AnaUtil::fillHist1D("nZdauflav",3,1);
             else if (pid == 15) AnaUtil::fillHist1D("nZdauflav",5,1);
             else AnaUtil::fillHist1D("nZdauflav",7,1);
           }                
          }
        }//gen particle loop
        std::sort(genZ.begin(),genZ.end(),MassComparator<vhtm::GenParticle>());
        AnaUtil::fillHist1D("genZ1mass",getP4(genZ[0]).M(),1); 
        AnaUtil::fillHist1D("genZ2mass",getP4(genZ[1]).M(),1);
        AnaUtil::fillHist1D("genHmassfromZ",(getP4(genZ[0])+getP4(genZ[1])).M(),1);
        //Z1 daughters
        if (getGenDauPgd(genZ[0]) == 13 ) AnaUtil::fillHist1D("nZ1dauflav",1,1);
        else if (getGenDauPgd(genZ[0]) == 11) AnaUtil::fillHist1D("nZ1dauflav",3,1);
        else if (getGenDauPgd(genZ[0]) == 15) AnaUtil::fillHist1D("nZ1dauflav",5,1);
        else AnaUtil::fillHist1D("nZ1dauflav",7,1);
        //Z2 daughter
        if (getGenDauPgd(genZ[1]) == 13 ) AnaUtil::fillHist1D("nZ2dauflav",1,1);
        else if (getGenDauPgd(genZ[1]) == 11) AnaUtil::fillHist1D("nZ2dauflav",3,1);
        else if (getGenDauPgd(genZ[1]) == 15) AnaUtil::fillHist1D("nZ2dauflav",5,1);
        else AnaUtil::fillHist1D("nZ2dauflav",7,1);
        //Event Flavour
        if( getGenDauPgd(genZ[0]) == 13 && getGenDauPgd(genZ[1]) == 13 ) {
          AnaUtil::fillHist1D("nEvntflav",1,1);
          AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
        }
        else if( getGenDauPgd(genZ[0]) == 11 && getGenDauPgd(genZ[1]) == 11 ) {
          AnaUtil::fillHist1D("nEvntflav",3,1);
         AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
        }
        else if( (getGenDauPgd(genZ[0]) == 13 && getGenDauPgd(genZ[1]) == 11) || 
                 (getGenDauPgd(genZ[0]) == 11 && getGenDauPgd(genZ[1]) == 13 ) ) {
          AnaUtil::fillHist1D("nEvntflav",5,1);
          AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
        } else {
          AnaUtil::fillHist1D("nEvntflav",7,1);
          //AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
        }
    }//event loop
 std::cout<<"Higgs to Z"<<nHtoZ<<std::endl;
 endJob();
}
///////////////////////////
template <typename T>
TLorentzVector GenAnalysis::getP4(const T& obj) {
   TLorentzVector l1;
   l1.SetPtEtaPhiE(obj.pt,obj.eta,obj.phi,obj.energy);
   return l1;
}

int GenAnalysis::getGenDauPgd(const GenParticle gp) {
 vector<int> d = gp.daughtIndices;
 vector<int> dpdgid;
 int ndau = 0;
 for (auto it = d.begin(); it != d.end(); ++it) {
   int di = (*it);
   if (di >= ngenparticle()) continue;
   const GenParticle& dgp = genParticleColl()->at(di);
   int pid = fabs(dgp.pdgId);
   if(pid==23) continue;
   dpdgid.push_back(pid);
 }
 if(!dpdgid.empty()) 
   return dpdgid[0];
 else return -1;
}

double GenAnalysis::getHmassfromZdau(const GenParticle Z1,const GenParticle Z2) {
 vector<int> d1 = Z1.daughtIndices;
 vector<int> d2 = Z2.daughtIndices;
 TLorentzVector lZ1,lZ2;
 lZ1.SetPtEtaPhiE(0.,0.,0.,0.);
 lZ2.SetPtEtaPhiE(0.,0.,0.,0.);
 //std::cout<<"Z1 dau size"<<d1.size()<<std::endl;
 for (auto it = d1.begin(); it != d1.end(); ++it) {
   int di = (*it);
   if (di >= ngenparticle()) continue;
   const GenParticle& dgp = genParticleColl()->at(di);
   int pid = fabs(dgp.pdgId);
   //std::cout<<"Z1 dau pdgId"<<pid<<std::endl; 
   if(pid==23) continue;
   lZ1+=getP4(dgp);
 }

 //std::cout<<"Z2 dau size"<<d2.size()<<std::endl;
 for (auto it = d2.begin(); it != d2.end(); ++it) {
   int di = (*it);
   if (di >= ngenparticle()) continue;
   const GenParticle& dgp = genParticleColl()->at(di);
   int pid = fabs(dgp.pdgId);
   //std::cout<<"Z2 dau pdgId"<<pid<<std::endl; 
   if(pid==23) continue;
   lZ2+=getP4(dgp);
 }
 //std::cout<<"H mass="<<(lZ1+lZ2).M()<<std::endl;
 return (lZ1+lZ2).M();
}
///////////////////////////
void GenAnalysis::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}



