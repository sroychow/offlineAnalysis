//
//  FSRStudy.cc
//  
//
//  Created by Suvankar Roy Chowdhury on 13/03/15.
//
//
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

#include "AnaUtil.h"
#include "PhysicsObjects.h"
#include "FSRStudy.h"

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
FSRStudy::FSRStudy()
  : PhysicsObjSelector(),
    _dumpEvent(false)
{
}
// ----------
// Destructor
// ----------
FSRStudy::~FSRStudy() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool FSRStudy::beginJob() 
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
void FSRStudy::bookHistograms() 
{
//-------muon-object selection plots
  new TH1F("mupt", "muon Pt", 150, 0., 150.);
//-------Electron Selection Plots-----------------------
  new TH1F("elept", "electron Pt", 150, 0., 150.);
//-------Photon Selection plots---------------------------------
  new TH1F("phopt", "Preselected pfphoton Pt", 150, 0., 150.);
//------------------------------------------------------
  new TH1D("nGoodmuon","Number of Good muons(with selection cuts) per event",20,-0.5,19.5);
  new TH1D("nGoodelectron","Number of Good electrons(with selection cuts) per event",20,-0.5,19.5);
    new TH1D("nGoodprePFPhoFSR","Number of pre-selected FSR photon per event",20,-0.5,19.5);
    new TH1D("nPhoFSRele","Number of FSR photon matched to an electron per event",20,-0.5,19.5);
    new TH1D("nPhoFSRmu","Number of FSR photon matched to a muon per event",20,-0.5,19.5);
    
//Event Level Cuts
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------
void FSRStudy::clearLists() {

 loosemuVec.clear();
 tightmuVec.clear();
 looseeleVec.clear();
 tighteleVec.clear();
 pfPhoFSR.clear();
 eleFSRPhopair.clear();
 muFSRPhopair.clear();
 PhysicsObjSelector::clear();
}
// -------------------
// The main event loop
// -------------------
void FSRStudy::eventLoop()
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
    //int n4mu=0,n2e2mu=0,n4e=0,nZ1clash=0,ntauev=0;
    for (int ev = 0; ev < nEvents(); ++ev) {
        
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
        //if(!muonColl()->empty())
            //  AnaUtil::fillHist1D("nRecomuon",muonColl()->size(),1);
          // std::cout<<"Start muon loop"<<std::endl;
        if( ngoodVtx < 1 ) continue;
        //PhysicsObjSelector pObj;
        std::cout<< " Event# " << setw(8) << ev << std::endl;
        std::cout << "FSR PAT mu/e" << muonColl()->size() << "/" << electronColl()->size()<<std::endl;
        PhysicsObjSelector::electronSelector();//std::cout<<" point 1>>";
        PhysicsObjSelector::muonSelector();//std::cout<<" point 2>>";
        PhysicsObjSelector::pfphotonSelector();//std::cout<<" point 3>>";
        pfPhoFSR = PhysicsObjSelector::getPFPhotonFSR();//std::cout<<" point 4>>";
        loosemuVec = PhysicsObjSelector::getLooseMu();//std::cout<<" point 5>>";
        looseeleVec = PhysicsObjSelector::getLooseEle();//std::cout<<" point 6>>";
        //std::cout<<"Nele="<<looseeleVec.size()<<std::endl;
        //std::cout<<"NPobjPATele="<<electronColl()->size()<<std::endl;
       // std::cout << "FSRLoose/Tight mu" << loosemuVec.size() << "/" << tightmuVec.size() << std::endl;
       // std::cout << "FSRLoose/Tight e" << looseeleVec.size() << "/" << tighteleVec.size() << std::endl;
        AnaUtil::fillHist1D("nGoodmuon",loosemuVec.size(),1);
        AnaUtil::fillHist1D("nGoodelectron",looseeleVec.size(),1);
        AnaUtil::fillHist1D("nGoodprePFPhoFSR",pfPhoFSR.size(),1);
        for( auto& ele : looseeleVec ) {
          TLorentzVector eleP4;
          eleP4.SetPtEtaPhiE(ele.pt,ele.eta,ele.phi,ele.energy);
          AnaUtil::fillHist1D("elept",ele.pt,1);
        }std::cout<<" point 7>>";
        for( auto& mu : loosemuVec ) {
          TLorentzVector muP4;
          muP4.SetPtEtaPhiE(mu.pt,mu.eta,mu.phi,mu.energy);
          AnaUtil::fillHist1D("mupt",mu.pt,1);
        }std::cout<<" point 8>>";
        for( auto& pho : pfPhoFSR ) {
           AnaUtil::fillHist1D("phopt",pho.pt,1);
        }std::cout<<" point 9>>";
        eleFSRPhopair = PhysicsObjSelector::getlooseElePhoFSR();
        muFSRPhopair = PhysicsObjSelector::getlooseMuPhoFSR();
        for( auto& elepho: eleFSRPhopair ) {
           AnaUtil::fillHist1D("nPhoFSRele",elepho.second.size(),1);
        }std::cout<<" point 10>>";
        for( auto& mupho: muFSRPhopair ) {
            AnaUtil::fillHist1D("nPhoFSRmu",mupho.second.size(),1);
        }std::cout<<" point 11>>";
        std::cout<<std::endl;
    }//event loop
  // Analysis is over
  endJob();
}
//////////New///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void FSRStudy::endJob() {

  closeFiles();
  
  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}



