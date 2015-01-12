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

#include "MuonIsolation.h"
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
MuonIsolation::MuonIsolation()
  : AnaBase(),
    _dumpEvent(false)
{
  
}
// ----------
// Destructor
// ----------
MuonIsolation::~MuonIsolation() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MuonIsolation::beginJob() 
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
void MuonIsolation::bookHistograms() 
{
 new TH1D("nRecomuon","Number of Reconstructed Muons(without cut)",20,0,20);
 new TH1D("nGoodmuon","Number of Good Muons(with selection cuts)",20,0,20);
 new TH1F("muon1Pt", "Pt of muon", 150, 0., 150.);
// new TH1F("muon1Pt", "Pt of sub-leading muon", 150, 0., 150.);
 new TH1F("muon1Eta", "Eta of muon", 30, -3., 3.);
//Charged Isolation
 new TH1F("muonChIso_015","Muon Charged Isolation for cone size 0.15",100,0.,10.);
 new TH1F("muonChIso_020","Muon Charged Isolation for cone size 0.20",100,0.,10.);
 new TH1F("muonChIso_025","Muon Charged Isolation for cone size 0.25",100,0.,10.);
 new TH1F("muonChIso_030","Muon Charged Isolation for cone size 0.30",100,0.,10.);
 new TH1F("muonChIso_035","Muon Charged Isolation for cone size 0.35",100,0.,10.);
 new TH1F("muonChIso_040","Muon Charged Isolation for cone size 0.40",100,0.,10.);
//Neutral Isolation
 new TH1F("muonNIso_015","Muon Neutral Isolation for cone size 0.15",100,0.,10.);
 new TH1F("muonNIso_020","Muon Neutral Isolation for cone size 0.20",100,0.,10.);
 new TH1F("muonNIso_025","Muon Neutral Isolation for cone size 0.25",100,0.,10.);
 new TH1F("muonNIso_030","Muon Neutral Isolation for cone size 0.30",100,0.,10.);
 new TH1F("muonNIso_035","Muon Neutral Isolation for cone size 0.35",100,0.,10.);
 new TH1F("muonNIso_040","Muon Neutral Isolation for cone size 0.40",100,0.,10.);
//Pileup 
 new TH1F("muonPuIso_015","Muon Isolation from Pileup for cone size 0.15",100,0.,10.);
 new TH1F("muonPuIso_020","Muon Isolation from Pileup for cone size 0.20",100,0.,10.);
 new TH1F("muonPuIso_025","Muon Isolation from Pileup for cone size 0.25",100,0.,10.);
 new TH1F("muonPuIso_030","Muon Isolation from Pileup for cone size 0.30",100,0.,10.);
 new TH1F("muonPuIso_035","Muon Isolation from Pileup for cone size 0.35",100,0.,10.);
 new TH1F("muonPuIso_040","Muon Isolation from Pileup for cone size 0.40",100,0.,10.);
//Relative Isolation
 new TH1F("muonrelIso_015","Muon Relative Isolation for cone size 0.15",100,0.,10.);
 new TH1F("muonrelIso_020","Muon Relative Isolation for cone size 0.20",100,0.,10.);
 new TH1F("muonrelIso_025","Muon Relative Isolation for cone size 0.25",100,0.,10.);
 new TH1F("muonrelIso_030","Muon Relative Isolation for cone size 0.30",100,0.,10.);
 new TH1F("muonrelIso_035","Muon Relative Isolation for cone size 0.35",100,0.,10.);
 new TH1F("muonrelIso_040","Muon Relative Isolation for cone size 0.40",100,0.,10.);
//Rho Values
 new TH1F("fGridRhoAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAllCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralChargedPileUp","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralNeutral","Grid Rho for event",100,0.,100.);
//Neutral Isolation vs Rho
 new TProfile("NIso015vsrho", "Neutral Isolation 015 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso020vsrho", "Neutral Isolation 020 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso025vsrho", "Neutral Isolation 025 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso030vsrho", "Neutral Isolation 030 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso035vsrho", "Neutral Isolation 035 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso040vsrho", "Neutral Isolation 040 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.); 
}
// -------------------
// The main event loop
// -------------------
void MuonIsolation::clearLists() {
 eventMuons.clear();
}
void MuonIsolation::eventLoop() 
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
  //long int nRecoMuons=0,nGoodRecoMuons=0;
  for (int ev = 0; ev < nEvents(); ++ev) {
    clearEvent();
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

   // AnaUtil::fillHist1D("evcounter", 0, _puevWt);
   AnaUtil::fillHist1D("fGridRhoAll",evt.fGridRhoAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAll",evt.fGridRhoFastjetAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAllCalo",evt.fGridRhoFastjetAllCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralCalo",evt.fGridRhoFastjetCentralCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralChargedPileUp",evt.fGridRhoFastjetCentralChargedPileUp,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralNeutral",evt.fGridRhoFastjetCentralNeutral,1);
   if(!muonColl()->empty())
     AnaUtil::fillHist1D("nRecomuon",muonColl()->size(),1);
   int ngoodmu=0;
   for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
     const Muon& muon = (*it);
     if (abs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta"))                      continue;
     if (muon.pt < AnaUtil::cutValue(muonCutMap(), "pt"))                            continue;
     if (!muon.isTrackerMuon)                                                         continue;
     if (!muon.isGlobalMuonPromptTight)                                               continue;
     if (muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers"))                continue;
     if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches"))                  continue;
     if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations"))  continue; 
     if (muon.pixHits < AnaUtil::cutValue(muonCutMap(), "pixHits"))                    continue;
     if (muon.trkHits < AnaUtil::cutValue(muonCutMap(), "trkHits"))                    continue;
     //if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap_,"globalChi2"))              continue;
     //if (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap_,"trkD0"))                   continue;
     ngoodmu++;
     AnaUtil::fillHist1D("muon1Pt",muon.pt,1);
     AnaUtil::fillHist1D("muon1Eta",muon.eta,1);
     if(muon.isolationMap.find(0.15) != muon.isolationMap.end() ) {
       AnaUtil::fillHist1D("muonChIso_015",muon.isolationMap.at(0.15).at(0)/muon.pt,1);
       AnaUtil::fillHist1D("muonNIso_015",muon.isolationMap.at(0.15).at(1)/muon.pt,1);
       AnaUtil::fillHist1D("muonPuIso_015",muon.isolationMap.at(0.15).at(2)/muon.pt,1);
       // do deltaBeta
       double iso = muon.isolationMap.at(0.15).at(0) + std::max(0.0, muon.isolationMap.at(0.15).at(1)-0.5*muon.isolationMap.at(0.15).at(2));
       AnaUtil::fillHist1D("muonrelIso_015",iso/muon.pt,1);
       AnaUtil::fillProfile("NIso015vsrho",evt.fGridRhoFastjetAll,muon.isolationMap.at(0.15).at(1));
     }
     if(muon.isolationMap.find(0.20) != muon.isolationMap.end() ) {
       AnaUtil::fillHist1D("muonChIso_020",muon.isolationMap.at(0.20).at(0)/muon.pt,1);
       AnaUtil::fillHist1D("muonNIso_020",muon.isolationMap.at(0.20).at(1)/muon.pt,1);
       AnaUtil::fillHist1D("muonPuIso_020",muon.isolationMap.at(0.20).at(2)/muon.pt,1);
       double iso = muon.isolationMap.at(0.20).at(0) + std::max(0.0, muon.isolationMap.at(0.20).at(1)-0.5*muon.isolationMap.at(0.20).at(2)); 
       AnaUtil::fillHist1D("muonrelIso_020",iso/muon.pt,1);
       AnaUtil::fillProfile("NIso020vsrho",evt.fGridRhoFastjetAll,muon.isolationMap.at(0.20).at(1));
    }
     if(muon.isolationMap.find(0.25) != muon.isolationMap.end() ) {
       AnaUtil::fillHist1D("muonChIso_025",muon.isolationMap.at(0.25).at(0)/muon.pt,1);
       AnaUtil::fillHist1D("muonNIso_025",muon.isolationMap.at(0.25).at(1)/muon.pt,1);
       AnaUtil::fillHist1D("muonPuIso_025",muon.isolationMap.at(0.25).at(2)/muon.pt,1);
       double iso = muon.isolationMap.at(0.25).at(0) + std::max(0.0, muon.isolationMap.at(0.25).at(1)-0.5*muon.isolationMap.at(0.25).at(2));
       AnaUtil::fillHist1D("muonrelIso_025",iso/muon.pt,1);
       AnaUtil::fillProfile("NIso025vsrho",evt.fGridRhoFastjetAll,muon.isolationMap.at(0.25).at(1));
     }
     if(muon.isolationMap.find(0.30) != muon.isolationMap.end() ) {
       AnaUtil::fillHist1D("muonChIso_030",muon.isolationMap.at(0.30).at(0)/muon.pt,1);
       AnaUtil::fillHist1D("muonNIso_030",muon.isolationMap.at(0.30).at(1)/muon.pt,1);
       AnaUtil::fillHist1D("muonPuIso_030",muon.isolationMap.at(0.30).at(2)/muon.pt,1);
       double iso = muon.isolationMap.at(0.30).at(0) + std::max(0.0, muon.isolationMap.at(0.30).at(1)-0.5*muon.isolationMap.at(0.30).at(2));
       AnaUtil::fillHist1D("muonrelIso_030",iso/muon.pt,1);
       AnaUtil::fillProfile("NIso030vsrho",evt.fGridRhoFastjetAll,muon.isolationMap.at(0.30).at(1));
     }
     if(muon.isolationMap.find(0.35) != muon.isolationMap.end() ) {
       AnaUtil::fillHist1D("muonChIso_035",muon.isolationMap.at(0.35).at(0)/muon.pt,1);
       AnaUtil::fillHist1D("muonNIso_035",muon.isolationMap.at(0.35).at(1)/muon.pt,1);
       AnaUtil::fillHist1D("muonPuIso_035",muon.isolationMap.at(0.35).at(2)/muon.pt,1);     
       double iso = muon.isolationMap.at(0.35).at(0) + std::max(0.0, muon.isolationMap.at(0.35).at(1)-0.5*muon.isolationMap.at(0.35).at(2));
       AnaUtil::fillHist1D("muonrelIso_035",iso/muon.pt,1);
       AnaUtil::fillProfile("NIso035vsrho",evt.fGridRhoFastjetAll,muon.isolationMap.at(0.35).at(1));
     }
     if(muon.isolationMap.find(0.40) != muon.isolationMap.end() ) {
       AnaUtil::fillHist1D("muonChIso_040",muon.isolationMap.at(0.4).at(0)/muon.pt,1);
       AnaUtil::fillHist1D("muonNIso_040",muon.isolationMap.at(0.4).at(1)/muon.pt,1);
       AnaUtil::fillHist1D("muonPuIso_040",muon.isolationMap.at(0.40).at(2)/muon.pt,1);
       double iso = muon.isolationMap.at(0.40).at(0) + std::max(0.0, muon.isolationMap.at(0.40).at(1)-0.5*muon.isolationMap.at(0.40).at(2));
       AnaUtil::fillHist1D("muonrelIso_040",iso/muon.pt,1);
       AnaUtil::fillProfile("NIso040vsrho",evt.fGridRhoFastjetAll,muon.isolationMap.at(0.40).at(1));
     }
   }
   if(ngoodmu > 0)
     AnaUtil::fillHist1D("nGoodmuon",ngoodmu,1);
  }  
  // Analysis is over
  endJob();
}

void MuonIsolation::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}
