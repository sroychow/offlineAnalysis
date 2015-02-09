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

#include "EventSelection.h"
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
EventSelection::EventSelection()
  : AnaBase(),
    _dumpEvent(false)
{
 //Mzdiff=9999.;
 //Z1present=false;
 //Z2present=false;
 //Z1tomu=false;
 //Z1toele=false;
 n4mu=0;
 n2e2mu=0;
 n4e=0;
 nZ1clash=0;
 ntauev=0;
}
// ----------
// Destructor
// ----------
EventSelection::~EventSelection() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool EventSelection::beginJob() 
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
void EventSelection::bookHistograms() 
{
//-------muon-object selection plots
  new TH1F("mupt", "muon Pt", 150, 0., 150.);
  new TH1F("muptid", "muon Pt(id passed)", 150, 0., 150.);
  new TH1F("muptidIso03", "muon Pt(id + pfRelIso03<0.4 passed)", 150, 0., 150.);
  new TH1F("muptidIso04", "muon Pt(id + pfRelIso04<0.4 passed)", 150, 0., 150.);
  new TH1F("muPFRelIso03","Muon PF Relative Isolation 03",100,0.,1.);
  new TH1F("muPFRelIso04","Muon PF Relative Isolation 04",100,0.,1.);
  new TH1F("muAbsSIP3D","Muon |SIP3D| ",30,0.,15.);
  new TH1D("muCutFlow","Muon cut flow",50,-0.5,49.5);
//------------------------------------------------------
//-------Electron Selection Plots-----------------------
  new TH1F("elept", "electron Pt", 150, 0., 150.);
  new TH1F("eleptid", "electron Pt(id passed)", 150, 0., 150.);
  new TH1F("eleptidIso", "muon Pt(id + pfRelIso<0.4 passed)", 150, 0., 150.);
  new TH1F("elePFRelIso","electron PF Relative Isolation",100,0.,1.);
  new TH1F("eleAbsSIP3D","Electron |SIP3D| ",30,0.,15.);
  new TH1D("eleCutFlow","Electron cut flow",50,-0.5,49.5);
//------------------------------------------------------
  new TH1D("nGoodmuon","Number of Good muons(with selection cuts) per event",20,0,20);
  new TH1D("nGoodelectron","Number of Good electrons(with selection cuts) per event",20,0,20);
  new TH1F("nl1Pt", "Pt of highest Pt lepton", 150, 0., 150.);
  new TH1F("nl2Pt", "Pt of second-highest Pt lepton", 150, 0., 150.);
  new TH1F("nl3Pt", "Pt of third-highest Pt lepton", 150, 0., 150.);
  new TH1F("nl4Pt", "Pt of lowest Pt lepton", 150, 0., 150.);
  new TH1F("nZ1cand", "Num of all possible Z1 candidates", 10, -0.5, 9.5);
  new TH1F("newmassZ1", "Mass of all possible Z1 candidates", 150, 0., 150.);
  new TH1F("newmassZ2", "Mass of all possible Z2 candidates", 150, 0., 150.);
  new TH1F("nZ2cand", "Num of all possible Z2 candidates", 10, -0.5, 9.5);
  new TH1F("newmassZ1select", "Mass of selected Z1 candidate(Z2 selection pending)", 150, 0., 150.);
  new TH1F("newmassZ1Z2select", "Mass of selected Z1 candidate when Z2 is selected", 150, 0., 150.);
  new TH1F("newmassZ2select", "Mass of selected Z2 candidates", 150, 0., 150.);
  new TH1F("newmassH", "Invariant Mass of Higgs", 500, 0., 500.);
  new TH1F("nZ1flav", "Flavour of selected Z1 candidates", 10, -0.5, 9.5);
  new TH1F("nZ2flav", "Flavour of selected Z2 candidates", 10, -0.5, 9.5);
  new TH1F("nEvflav", "Flavour of selected event", 10, -0.5, 9.5);
//GEN Level histograms
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
  //new TH1D("genPass","genPass",10,-0.5,9.5);
//Event Level Cuts
  new TH1D("evntCutFlow","Event cut flow",50,-0.5,49.5);
  new TH1D("dRgenZ1recoZ1","Event cut flow",200,0.,1.);
  new TH1D("dRgenZ2recoZ2","Event cut flow",200,0.,1.);
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------
void EventSelection::clearLists() {
 vtxList.clear();
 eleVec.clear();
 muVec.clear();
 tauVec.clear();
 myZ1.clear();
 myZ2.clear();
 genZ.clear();
}
// -------------------
// The main event loop
// -------------------
void EventSelection::eventLoop()
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
        //if(!muonColl()->empty())
            //  AnaUtil::fillHist1D("nRecomuon",muonColl()->size(),1);
          // std::cout<<"Start muon loop"<<std::endl;
        AnaUtil::fillHist1D("evntCutFlow",0,1);
        bool genPass = genOk();
        if(!genPass) continue;
        AnaUtil::fillHist1D("evntCutFlow",1,1);
        muonSelector();
        electronSelector();
        tauSelector();
        if( eleVec.size() + muVec.size() > 4 )   AnaUtil::fillHist1D("evntCutFlow",2,1); 
        if( !eleVec.empty() || !muVec.empty() )   AnaUtil::fillHist1D("evntCutFlow",3,1);
        if( !muVec.empty() )   AnaUtil::fillHist1D("evntCutFlow",4,1); 
        if( !eleVec.empty() )   AnaUtil::fillHist1D("evntCutFlow",5,1); 
        if( !tauVec.empty() )   AnaUtil::fillHist1D("evntCutFlow",6,1); 
        if( muVec.size() >= 2 )   AnaUtil::fillHist1D("evntCutFlow",7,1); 
        if( eleVec.size() >= 2 )   AnaUtil::fillHist1D("evntCutFlow",8,1);
        if( muVec.size() >=2 && eleVec.size() >=2 )   AnaUtil::fillHist1D("evntCutFlow",9,1); 
        if( muVec.size() >=4 )   AnaUtil::fillHist1D("evntCutFlow",10,1);
        if( eleVec.size() >=4 )   AnaUtil::fillHist1D("evntCutFlow",11,1);
        if( !muVec.empty() ) {
          AnaUtil::fillHist1D("nGoodmuon",muVec.size(),1);
          std::sort(muVec.begin(),muVec.end(),PtComparator<vhtm::Muon>());
          ZSelector<vhtm::Muon>(muVec);
        }
        if( !eleVec.empty() ) {
          AnaUtil::fillHist1D("nGoodelectron",eleVec.size(),1);
          std::sort(eleVec.begin(),eleVec.end(),PtComparator<vhtm::Electron>());
          ZSelector<vhtm::Electron>(eleVec);
        }
        plotHmass();
    }
  // Analysis is over
  endJob();
}
//----
//Function to select muons
//id + iso condtions applied
//----
void EventSelection::muonSelector() {
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);
    AnaUtil::fillHist1D("muCutFlow",0,1);
    AnaUtil::fillHist1D("mupt",muon.pt,1);
    if (abs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta"))                      continue;
    AnaUtil::fillHist1D("muCutFlow",1,1);
    if (muon.pt < AnaUtil::cutValue(muonCutMap(), "pt"))                              continue;
    AnaUtil::fillHist1D("muCutFlow",2,1);
    if (!muon.isTrackerMuon)                                                          continue;
    AnaUtil::fillHist1D("muCutFlow",3,1);
    if (!muon.isGlobalMuonPromptTight)                                                continue;
    AnaUtil::fillHist1D("muCutFlow",4,1);
    if (muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers"))                continue;
    AnaUtil::fillHist1D("muCutFlow",5,1);
    if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches"))                  continue;
    AnaUtil::fillHist1D("muCutFlow",6,1);
    if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations"))  continue;
    AnaUtil::fillHist1D("muCutFlow",7,1);
    if (muon.pixHits < AnaUtil::cutValue(muonCutMap(), "pixHits"))                    continue;
    AnaUtil::fillHist1D("muCutFlow",8,1);
    if (muon.trkHits < AnaUtil::cutValue(muonCutMap(), "trkHits"))                    continue;
    AnaUtil::fillHist1D("muCutFlow",9,1);
    if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2"))              continue;
    AnaUtil::fillHist1D("muCutFlow",10,1);
    if (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))                   continue;
    AnaUtil::fillHist1D("muCutFlow",11,1);
    AnaUtil::fillHist1D("muptid",muon.pt,1);
    AnaUtil::fillHist1D("muPFRelIso03",muon.pfRelIso03,1);
    AnaUtil::fillHist1D("muPFRelIso04",muon.pfRelIso04,1);
    if( muon.pfRelIso03 < 0.4)
      AnaUtil::fillHist1D("muptidIso03",muon.pt,1);
    if( muon.pfRelIso04 < 0.4)
      AnaUtil::fillHist1D("muptidIso04",muon.pt,1);
    if( muon.pfRelIso04 > 0.4 )                                                       continue;
    AnaUtil::fillHist1D("muCutFlow",12,1);
    AnaUtil::fillHist1D("muAbsSIP3D",std::fabs(muon.dB3D/muon.edB3D),1);
    if(std::fabs(muon.dB3D/muon.edB3D) > 4. ) continue; 
    AnaUtil::fillHist1D("muCutFlow",13,1);
     muVec.push_back(muon);
  }

}
//----
//Function to select electrons
//id + iso condtions applied
//----
void EventSelection::electronSelector() {
  for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it) {
    const Electron& electron = (*it);
    AnaUtil::fillHist1D("eleCutFlow",0,1);
    AnaUtil::fillHist1D("elept",electron.pt,1);
    if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt")) continue;
    AnaUtil::fillHist1D("eleCutFlow",1,1);
    double eleta  = std::fabs(electron.eta);
    bool   etaCut = (eleta >= AnaUtil::cutValue(electronCutMap(), "etaLow") &&
	             eleta <= AnaUtil::cutValue(electronCutMap(), "etaUp")) ||
                     eleta  >= AnaUtil::cutValue(electronCutMap(), "eta");
    if (etaCut) continue;
    AnaUtil::fillHist1D("eleCutFlow",2,1);
    if(!eleId(electron,bunchCrossing())) continue; 
    AnaUtil::fillHist1D("eleCutFlow",3,1);
    if (std::fabs(electron.trkD0) >= AnaUtil::cutValue(electronCutMap(),"trkD0")) continue;
    AnaUtil::fillHist1D("eleCutFlow",4,1);
    if (std::fabs(electron.vtxDistZ) >= AnaUtil::cutValue(electronCutMap(), "dz")) continue;
    AnaUtil::fillHist1D("eleCutFlow",5,1);
    AnaUtil::fillHist1D("eleptid",electron.pt,1);
    AnaUtil::fillHist1D("elePFRelIso",electron.pfRelIso,1);
    if (electron.pfRelIso > 0.4) continue;
    AnaUtil::fillHist1D("eleCutFlow",6,1);
    AnaUtil::fillHist1D("eleptidIso",electron.pt,1);
    AnaUtil::fillHist1D("eleAbsSIP3D",std::fabs(electron.dB3D/electron.edB3D),1);
    if( std::fabs(electron.dB3D/electron.edB3D)  > 4. ) continue; 
    AnaUtil::fillHist1D("eleCutFlow",7,1);
    eleVec.push_back(electron);
  }
}
//Tau Selector
void EventSelection::tauSelector() {
  int taucounters[] = {0,0,0,0,0,0};
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);
    
    //    std::cout << "ptTau1 = " << AnaUtil::cutValue(_evselCutMap,"ptTau1")) << std::endl;
    if (tau.pt < AnaUtil::cutValue(evselCutMap(),"ptTau1")) continue;
    ++taucounters[0];
    if( std::fabs(tau.eta) >= AnaUtil::cutValue(evselCutMap(),"etaTau1")) continue;
    ++taucounters[1];
    
    if (tau.decayModeFinding != 1.0) continue;
    ++taucounters[2];
    
    if(tau.chargedIsoPtSum >= 2) continue; /// loose combined Isolation
    ++taucounters[3];
    
    if (tau.againstMuonLoose3 <= 0.5) continue;
    ++taucounters[4];
    
    if (tau.againstElectronLooseMVA5 <= 0.5) continue;
    ++taucounters[5];
    
    
    tauVec.push_back(tau);
  }
}
//////////New///////////////////////////////////////////////////////////////////////////////////
template <typename T>
void EventSelection::ZSelector(const std::vector<T>& lepVec) {
//Z1 selection loop:: Z1 selection is based on how close mupair mass is close
//to nominal Z mass
 bool Z1found=false,Z2found=false,dup;
 for(auto& mu1: lepVec ) {
   dup=false;
   TLorentzVector mu1P4;
   mu1P4.SetPtEtaPhiE(mu1.pt,mu1.eta,mu1.phi,mu1.energy);
   for(auto& mu2: lepVec ) {
     if( mu1.charge + mu2.charge != 0) continue;
      //std::cout<<"Zselector point 1>>>"<<std::endl;
     TLorentzVector mu2P4;
     mu2P4.SetPtEtaPhiE(mu2.pt,mu2.eta,mu2.phi,mu2.energy);
       if(!myZ1.empty()){
           Z1found=true;//already at least one Z1 is there. maybe of same flavour maybe different
           for(int i=0;i<myZ1.size();i++) {
               if( (mu1P4 == myZ1[i].l1P && mu2P4 == myZ1[i].l2P) || (mu1P4 == myZ1[i].l2P && mu2P4 == myZ1[i].l1P) ) {
                   dup=true;
                   break;
               }
           }
       }
     if(dup) continue;
     double MZ1temp=(mu1P4+mu2P4).M();
     if( MZ1temp >= AnaUtil::cutValue(evselCutMap(),"mZ1lowbound") && MZ1temp <= AnaUtil::cutValue(evselCutMap(),"mZ1upperbound") ) {
         //std::cout<<"Zselector point 2>>>"<<std::endl;
       Z1found=true;
       Zcand ztmp;
       ztmp.l1P=mu1P4;
       ztmp.l2P=mu2P4;
       ztmp.chargel1=mu1.charge;
       ztmp.chargel2=mu2.charge;
       if( typeid(mu1)==typeid(vhtm::Muon) ) 
         ztmp.flavour=0;
       else if( typeid(mu1)==typeid(vhtm::Electron) ) 
         ztmp.flavour=1;
       else
         ztmp.flavour=2; 
         ztmp.mass=MZ1temp;
         ztmp.massDiff=std::fabs(MZ1temp-MZnominal) ;
         myZ1.push_back(ztmp);
     }
   }
 }
 //std::cout<<"Zselector point 3>>>"<<std::endl;
 if(Z1found) {//Z2selector
   AnaUtil::fillHist1D("nZ1cand",myZ1.size(),1);  
   AnaUtil::fillHist1D("evntCutFlow",12,1); 
   for(auto& z:myZ1) {
     
     AnaUtil::fillHist1D("newmassZ1",(z.l1P+z.l2P).M(),1);
   }
   std::sort(myZ1.begin(),myZ1.end(),ZSorter);//sort the z1 found according to mass diff with nominal Z mass
     AnaUtil::fillHist1D("newmassZ1select",(myZ1[0].l1P+myZ1[0].l2P).M(),1);
     for(int i =0;i<myZ1.size();i++)
     {
         std::cout<<"index="<<i<<">>>Mass="<<(myZ1[i].l1P+myZ1[i].l2P).M()<<">>P41="<<myZ1[i].l1P.Pt()<<">>P41="<<myZ1[i].l2P.Pt()<<std::endl;
     }
     //std::cout<<"Zselector point 4>>>"<<std::endl;
   for(auto& mu1: lepVec ) {
     TLorentzVector mu1P4;
     mu1P4.SetPtEtaPhiE(mu1.pt,mu1.eta,mu1.phi,mu1.energy);
     if( mu1P4==myZ1[0].l1P || mu1P4==myZ1[0].l2P ) continue;
       std::cout<<"Zselector point 5>>>"<<std::endl;
     for(auto& mu2: lepVec ) {
       //std::cout<<"Opp charged lepton found>>>"<<std::endl;
       TLorentzVector mu2P4;
       mu2P4.SetPtEtaPhiE(mu2.pt,mu2.eta,mu2.phi,mu2.energy);
       if( (mu2P4==myZ1[0].l1P || mu2P4==myZ1[0].l2P) || mu1.charge + mu2.charge != 0) continue;
         std::cout<<"Zselector point 6>>>"<<std::endl;
       double MZ2temp=(mu1P4+mu2P4).M();
       if( MZ2temp >= AnaUtil::cutValue(evselCutMap(),"mZ2lowbound") && MZ2temp <= AnaUtil::cutValue(evselCutMap(),"mZ2upperbound") ) {
           std::cout<<"Zselector point 7>>>"<<std::endl;
         Z2found=true;
         Zcand ztmp;
         ztmp.l1P=mu1P4;
         ztmp.l2P=mu2P4;
         ztmp.chargel1=mu1.charge;
         ztmp.chargel2=mu2.charge;
         if( typeid(mu1)==typeid(vhtm::Muon) ) 
           ztmp.flavour=0;
         else if( typeid(mu1)==typeid(vhtm::Electron) ) 
           ztmp.flavour=1;
         else
           ztmp.flavour=2; 
           ztmp.mass=MZ2temp;
           ztmp.massDiff=std::fabs(MZ2temp-MZnominal) ;
           myZ2.push_back(ztmp);
       }       
     }  
   }
 }//Z2selector
 else AnaUtil::fillHist1D("nZ1cand",0,1);  
    if(Z2found)  {
        AnaUtil::fillHist1D("nZ2cand",myZ2.size(),1);
        AnaUtil::fillHist1D("newmassZ1Z2select",(myZ1[0].l1P+myZ1[0].l2P).M(),1);
        AnaUtil::fillHist1D("evntCutFlow",13,1); 
        for(auto& z:myZ2) {
            
            AnaUtil::fillHist1D("newmassZ2",(z.l1P+z.l2P).M(),1);
        }
    }
 else
   AnaUtil::fillHist1D("nZ2cand",0,1);
    
    if( myZ2.size() > 1 ) {
        std::sort(myZ2.begin(),myZ2.end(),ZSorter);
        Z2actual=myZ2[0];
        std::cout<<"Z2actual>>>Mass="<<(Z2actual.l1P+Z2actual.l2P).M()<<">>P41="<<Z2actual.l1P.Pt()<<">>P42="<<Z2actual.l2P.Pt()<<std::endl;
        std::cout<<"Zselector point 8>>>"<<std::endl;
        for( int i = 1; i<myZ2.size(); i++ ) {
          TLorentzVector mu1P4,mu2P4;
            mu1P4=myZ2[i].l1P;
            mu2P4=myZ2[i].l2P;
         if(   ( myZ2[i].chargel1 == Z2actual.chargel1 && myZ2[i].l1P.Pt() > Z2actual.l1P.Pt() ) ||
               ( myZ2[i].chargel1 == Z2actual.chargel2 && myZ2[i].l1P.Pt() > Z2actual.l2P.Pt() ) ||
               ( myZ2[i].chargel2 == Z2actual.chargel1 && myZ2[i].l2P.Pt() > Z2actual.l1P.Pt() ) ||
               ( myZ2[i].chargel2 == Z2actual.chargel2 && myZ2[i].l2P.Pt() > Z2actual.l2P.Pt() )
            ) {
             Z2actual=myZ2[i];
             std::cout<<"Zselector point 9>>>"<<std::endl;
         }
        }
        AnaUtil::fillHist1D("newmassZ2select",(Z2actual.l1P+Z2actual.l2P).M(),1);
    }
    std::cout<<"Zselector point 10>>>"<<std::endl;
}

void EventSelection::plotHmass() {
    if(!myZ1.empty() && !myZ2.empty()){
        std::cout<<"plotHmass() point 1>>>"<<std::endl;
        std::vector<TLorentzVector> eventLep;
        eventLep.push_back(myZ1[0].l1P);
        eventLep.push_back(myZ1[0].l2P);
        eventLep.push_back(Z2actual.l1P);
        eventLep.push_back(Z2actual.l2P);
        std::sort(eventLep.begin(),eventLep.end(),PtComparatorTL<TLorentzVector>());
        std::cout<<"plotHmass() point 2>>>"<<std::endl;
        if ( eventLep[0].Pt() > 20. && eventLep[1].Pt() > 10. ) {
            std::cout<<"plotHmass() point 3>>>"<<std::endl;
            AnaUtil::fillHist1D("evntCutFlow",14,1); 
            std::cout<<"lep pts>>>"<<eventLep[0].Pt()<<":"<<eventLep[1].Pt()<<":"<<eventLep[2].Pt()<<":"<<eventLep[3].Pt()<<std::endl;
            AnaUtil::fillHist1D("nl1Pt",eventLep[0].Pt(),1);
            AnaUtil::fillHist1D("nl2Pt",eventLep[1].Pt(),1);
            AnaUtil::fillHist1D("nl3Pt",eventLep[2].Pt(),1);
            AnaUtil::fillHist1D("nl4Pt",eventLep[3].Pt(),1);
            AnaUtil::fillHist1D("newmassH",(myZ1[0].l1P+myZ1[0].l2P+Z2actual.l1P+Z2actual.l2P).M(),1);
            AnaUtil::fillHist1D("nZ1flav",myZ1[0].flavour,1);
            AnaUtil::fillHist1D("nZ2flav",Z2actual.flavour,1);

            AnaUtil::fillHist1D("dRgenZ1recoZ1",(myZ1[0].l1P+myZ1[0].l2P).DeltaR(getP4(genZ[0])),1);
            AnaUtil::fillHist1D("dRgenZ2recoZ2",(Z2actual.l1P+Z2actual.l2P).DeltaR(getP4(genZ[1])),1);
            if( myZ1[0].flavour == 0 && Z2actual.flavour == 0)
              AnaUtil::fillHist1D("nEvflav",1,1);//4mu
            else if( myZ1[0].flavour == 1 && Z2actual.flavour == 1)
              AnaUtil::fillHist1D("nEvflav",3,1);//4e
            else if( (myZ1[0].flavour == 0 && Z2actual.flavour == 1) || (myZ1[0].flavour == 1 && Z2actual.flavour == 0) )
              AnaUtil::fillHist1D("nEvflav",5,1);//2e2mu
            else
              AnaUtil::fillHist1D("nEvflav",7,1);//wrong 
        }
    }
    std::cout<<"plotHmass() point 4   >>>"<<std::endl;
}
///////////////////////////////////
//Gen Level Functions
//////////////////////////////////
template <typename T>
TLorentzVector EventSelection::getP4(const T& obj) {
   TLorentzVector l1;
   l1.SetPtEtaPhiE(obj.pt,obj.eta,obj.phi,obj.energy);
   return l1;
}

int EventSelection::getGenDauPgd(const GenParticle gp) {
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

double EventSelection::getHmassfromZdau(const GenParticle Z1,const GenParticle Z2) {
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
bool EventSelection::genOk() {
  bool genPass=false;
  int nHtoZ=0;
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
//          if (pdgid == 23 && status==2 ) {
          if (pdgid == 23 && status == 3 ) {
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
          genPass=true;
        }
        else if( getGenDauPgd(genZ[0]) == 11 && getGenDauPgd(genZ[1]) == 11 ) {
          AnaUtil::fillHist1D("nEvntflav",3,1);
          AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
          genPass=true;
        }
        else if( (getGenDauPgd(genZ[0]) == 13 && getGenDauPgd(genZ[1]) == 11) || 
                 (getGenDauPgd(genZ[0]) == 11 && getGenDauPgd(genZ[1]) == 13 ) ) {
          AnaUtil::fillHist1D("nEvntflav",5,1);
          AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
          genPass=true;
        } else {
          AnaUtil::fillHist1D("nEvntflav",7,1);
          //AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
        }
  return genPass;
}
//////////////////////////////////////////////////////////////////////////////////
void EventSelection::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}


