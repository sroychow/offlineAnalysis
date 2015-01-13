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

#include "ElectronIsolation.h"
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
ElectronIsolation::ElectronIsolation()
  : AnaBase(),
    _dumpEvent(false)
{
  
}
// ----------
// Destructor
// ----------
ElectronIsolation::~ElectronIsolation() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ElectronIsolation::beginJob() 
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
void ElectronIsolation::bookHistograms() 
{
 new TH1D("nRecoelectron","Number of Reconstructed electrons(without cut)",20,0,20);
 new TH1D("nGoodelectron","Number of Good electrons(with selection cuts)",20,0,20);
 new TH1F("electron1Pt", "Pt of electron", 150, 0., 150.);
// new TH1F("electron1Pt", "Pt of sub-leading electron", 150, 0., 150.);
 new TH1F("electron1Eta", "Eta of electron", 30, -3., 3.);
//Charged Isolation
 new TH1F("electronChIso_015","electron Charged Isolation for cone size 0.15",100,0.,10.);
 new TH1F("electronChIso_020","electron Charged Isolation for cone size 0.20",100,0.,10.);
 new TH1F("electronChIso_025","electron Charged Isolation for cone size 0.25",100,0.,10.);
 new TH1F("electronChIso_030","electron Charged Isolation for cone size 0.30",100,0.,10.);
 new TH1F("electronChIso_035","electron Charged Isolation for cone size 0.35",100,0.,10.);
 new TH1F("electronChIso_040","electron Charged Isolation for cone size 0.40",100,0.,10.);
//Neutral Isolation
 new TH1F("electronNIso_015","electron Neutral Isolation for cone size 0.15",100,0.,10.);
 new TH1F("electronNIso_020","electron Neutral Isolation for cone size 0.20",100,0.,10.);
 new TH1F("electronNIso_025","electron Neutral Isolation for cone size 0.25",100,0.,10.);
 new TH1F("electronNIso_030","electron Neutral Isolation for cone size 0.30",100,0.,10.);
 new TH1F("electronNIso_035","electron Neutral Isolation for cone size 0.35",100,0.,10.);
 new TH1F("electronNIso_040","electron Neutral Isolation for cone size 0.40",100,0.,10.);
//Pileup 
 new TH1F("electronPuIso_015","electron Isolation from Pileup for cone size 0.15",100,0.,10.);
 new TH1F("electronPuIso_020","electron Isolation from Pileup for cone size 0.20",100,0.,10.);
 new TH1F("electronPuIso_025","electron Isolation from Pileup for cone size 0.25",100,0.,10.);
 new TH1F("electronPuIso_030","electron Isolation from Pileup for cone size 0.30",100,0.,10.);
 new TH1F("electronPuIso_035","electron Isolation from Pileup for cone size 0.35",100,0.,10.);
 new TH1F("electronPuIso_040","electron Isolation from Pileup for cone size 0.40",100,0.,10.);
//Relative Isolation
 new TH1F("electronrelIso_015","electron Relative Isolation for cone size 0.15",100,0.,10.);
 new TH1F("electronrelIso_020","electron Relative Isolation for cone size 0.20",100,0.,10.);
 new TH1F("electronrelIso_025","electron Relative Isolation for cone size 0.25",100,0.,10.);
 new TH1F("electronrelIso_030","electron Relative Isolation for cone size 0.30",100,0.,10.);
 new TH1F("electronrelIso_035","electron Relative Isolation for cone size 0.35",100,0.,10.);
 new TH1F("electronrelIso_040","electron Relative Isolation for cone size 0.40",100,0.,10.);
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
void ElectronIsolation::clearLists() {
 eventElectrons.clear();
}
void ElectronIsolation::eventLoop() 
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
  //long int nRecoelectrons=0,nGoodRecoelectrons=0;
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
   if(!electronColl()->empty())
     AnaUtil::fillHist1D("nRecoelectron",electronColl()->size(),1);
   int ngoodmu=0;
   for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it) {
     const Electron& electron = (*it);

     bool quality_EB_loose = (abs(electron.eta) <= 1.4442 
                         &&  electron.sigmaEtaEta < 0.01 
                         &&  electron.deltaEtaTrkSC < 0.007 
                         &&  electron.deltaPhiTrkSC < 0.8
                         &&  electron.hoe < 0.15);
     bool quality_EE_loose = (abs(electron.eta) >= 1.566 
                         &&  electron.sigmaEtaEta < 0.03
                         &&  electron.deltaEtaTrkSC < 0.01 
                         &&  electron.deltaPhiTrkSC < 0.7
                         &&  electron.hoe < 0.07);
     bool quality_loose = quality_EB_loose || quality_EE_loose;  

     //if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt"))                  continue;
     if (abs(electron.eta) >= AnaUtil::cutValue(electronCutMap(), "eta"))           continue;
     if (!quality_loose)                                                       continue;
     if (electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits")) continue;
    
     ngoodmu++;
     AnaUtil::fillHist1D("electron1Pt",electron.pt,1);
     AnaUtil::fillHist1D("electron1Eta",electron.eta,1);
     if(electron.isolationMap.find(0.15) != electron.isolationMap.end() ) {
       AnaUtil::fillHist1D("electronChIso_015",electron.isolationMap.at(0.15).at(0)/electron.pt,1);
       AnaUtil::fillHist1D("electronNIso_015",electron.isolationMap.at(0.15).at(1)/electron.pt,1);
       AnaUtil::fillHist1D("electronPuIso_015",electron.isolationMap.at(0.15).at(2)/electron.pt,1);
       // do deltaBeta
       double iso = electron.isolationMap.at(0.15).at(0) + std::max(0.0, electron.isolationMap.at(0.15).at(1)-0.5*electron.isolationMap.at(0.15).at(2));
       AnaUtil::fillHist1D("electronrelIso_015",iso/electron.pt,1);
       AnaUtil::fillProfile("NIso015vsrho",evt.fGridRhoFastjetAll,electron.isolationMap.at(0.15).at(1));
     }
     if(electron.isolationMap.find(0.20) != electron.isolationMap.end() ) {
       AnaUtil::fillHist1D("electronChIso_020",electron.isolationMap.at(0.20).at(0)/electron.pt,1);
       AnaUtil::fillHist1D("electronNIso_020",electron.isolationMap.at(0.20).at(1)/electron.pt,1);
       AnaUtil::fillHist1D("electronPuIso_020",electron.isolationMap.at(0.20).at(2)/electron.pt,1);
       double iso = electron.isolationMap.at(0.20).at(0) + std::max(0.0, electron.isolationMap.at(0.20).at(1)-0.5*electron.isolationMap.at(0.20).at(2)); 
       AnaUtil::fillHist1D("electronrelIso_020",iso/electron.pt,1);
       AnaUtil::fillProfile("NIso020vsrho",evt.fGridRhoFastjetAll,electron.isolationMap.at(0.20).at(1));
    }
     if(electron.isolationMap.find(0.25) != electron.isolationMap.end() ) {
       AnaUtil::fillHist1D("electronChIso_025",electron.isolationMap.at(0.25).at(0)/electron.pt,1);
       AnaUtil::fillHist1D("electronNIso_025",electron.isolationMap.at(0.25).at(1)/electron.pt,1);
       AnaUtil::fillHist1D("electronPuIso_025",electron.isolationMap.at(0.25).at(2)/electron.pt,1);
       double iso = electron.isolationMap.at(0.25).at(0) + std::max(0.0, electron.isolationMap.at(0.25).at(1)-0.5*electron.isolationMap.at(0.25).at(2));
       AnaUtil::fillHist1D("electronrelIso_025",iso/electron.pt,1);
       AnaUtil::fillProfile("NIso025vsrho",evt.fGridRhoFastjetAll,electron.isolationMap.at(0.25).at(1));
     }
     if(electron.isolationMap.find(0.30) != electron.isolationMap.end() ) {
       AnaUtil::fillHist1D("electronChIso_030",electron.isolationMap.at(0.30).at(0)/electron.pt,1);
       AnaUtil::fillHist1D("electronNIso_030",electron.isolationMap.at(0.30).at(1)/electron.pt,1);
       AnaUtil::fillHist1D("electronPuIso_030",electron.isolationMap.at(0.30).at(2)/electron.pt,1);
       double iso = electron.isolationMap.at(0.30).at(0) + std::max(0.0, electron.isolationMap.at(0.30).at(1)-0.5*electron.isolationMap.at(0.30).at(2));
       AnaUtil::fillHist1D("electronrelIso_030",iso/electron.pt,1);
       AnaUtil::fillProfile("NIso030vsrho",evt.fGridRhoFastjetAll,electron.isolationMap.at(0.30).at(1));
     }
     if(electron.isolationMap.find(0.35) != electron.isolationMap.end() ) {
       AnaUtil::fillHist1D("electronChIso_035",electron.isolationMap.at(0.35).at(0)/electron.pt,1);
       AnaUtil::fillHist1D("electronNIso_035",electron.isolationMap.at(0.35).at(1)/electron.pt,1);
       AnaUtil::fillHist1D("electronPuIso_035",electron.isolationMap.at(0.35).at(2)/electron.pt,1);     
       double iso = electron.isolationMap.at(0.35).at(0) + std::max(0.0, electron.isolationMap.at(0.35).at(1)-0.5*electron.isolationMap.at(0.35).at(2));
       AnaUtil::fillHist1D("electronrelIso_035",iso/electron.pt,1);
       AnaUtil::fillProfile("NIso035vsrho",evt.fGridRhoFastjetAll,electron.isolationMap.at(0.35).at(1));
     }
     if(electron.isolationMap.find(0.40) != electron.isolationMap.end() ) {
       AnaUtil::fillHist1D("electronChIso_040",electron.isolationMap.at(0.4).at(0)/electron.pt,1);
       AnaUtil::fillHist1D("electronNIso_040",electron.isolationMap.at(0.4).at(1)/electron.pt,1);
       AnaUtil::fillHist1D("electronPuIso_040",electron.isolationMap.at(0.40).at(2)/electron.pt,1);
       double iso = electron.isolationMap.at(0.40).at(0) + std::max(0.0, electron.isolationMap.at(0.40).at(1)-0.5*electron.isolationMap.at(0.40).at(2));
       AnaUtil::fillHist1D("electronrelIso_040",iso/electron.pt,1);
       AnaUtil::fillProfile("NIso040vsrho",evt.fGridRhoFastjetAll,electron.isolationMap.at(0.40).at(1));
     }
   }
   if(ngoodmu > 0)
     AnaUtil::fillHist1D("nGoodelectron",ngoodmu,1);
  }  
  // Analysis is over
  endJob();
}

void ElectronIsolation::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}
