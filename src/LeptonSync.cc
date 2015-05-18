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

#include "LeptonSync.h"
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
 LeptonSync::LeptonSync()
: AnaBase(),
 _dumpEvent(false)
{
 evtype = -1;
 lepsyncDumper_.open("lepsync.txt");
}
// ----------
// Destructor
// ----------
LeptonSync::~LeptonSync() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool LeptonSync::beginJob() 
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
void LeptonSync::bookHistograms()
{
    new TH1D("nGoodvtx","Number of Vertices",50,0,50);
    new TH1D("muCutFlow","",30,-0.5,29.5);
    new TH1D("eleCutFlow","",30,-0.5,29.5);
    new TH1D("mud0","",500,0,1.);
    new TH1D("mud0PV","",500,0,1.);
    new TH1D("mudz","",200,-20.,20.);
    new TH1D("mudzPV","",200,-20.,20.);
    new TH1D("eled0","",500,0,1.);
    new TH1D("eled0PV","",500,0,1.);
    new TH1D("eledz","",200,-20.,20.);
    new TH1D("eledzPV","",200,-20.,20.);
    new TH1D("nloosemu_ev4mu","No of loose mu selected where gen level is 4mu state",8,-0.5,7.5);
    new TH1D("nlooseele_ev4e","No of loose ele selected where gen level is 4ele state",8,-0.5,7.5);
    new TH1D("nloosemu_ev2e2mu","No of loose mu selected where gen level is 2e2mu state",8,-0.5,7.5);
    new TH1D("nlooseele_ev2e2mu","No of loose ele selected where gen level is 2e2mu state",8,-0.5,7.5);
    new TH1D("formZZ_ev4mu","N events where nloose mu >= 4 when gen is 4mu",2,0,2);
    new TH1D("formZZ_ev4e","N events where nloose ele >= 4 when gen is 4ele",2,0,2);
    new TH1D("formZZ_ev2e2mu","N events where nloose (mu >= 2 && ele>=2) when gen is 2e2mu",2,0,2);
}


// -------------------
// The main event loop
// -------------------
void LeptonSync::clearLists() {
 evtype = -1;
 vtxList.clear();
 genZ.clear();
    looseMuPhoFSRpair_.clear();
    looseElePhoFSRpair_.clear();
    loosemuVec_.clear();
    tightmuVec_.clear();
    looseeleVec_.clear();
    tighteleVec_.clear();
    pfphotonVec_.clear();
}
void LeptonSync::eventLoop() 
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
 std::cout<<"Bunch Crossing>>>>"<<bunchCrossing()<<std::endl;
    
 for (int ev = 0; ev < nEvents(); ++ev){
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
  run_ = run; 
  lumi_ = lumis; 
  event_ = event;
  rho_  = evt.fGridRhoFastjetAll;
  //lepsyncDumper_ << run << ":" << lumis << ":" << event << ":" << evt.fGridRhoFastjetAll << endl;

  op.verbose =(logOption() >> 1 & 0x1);
  findVtxInfo(vtxList, op, fLog());
  if( !genOk() ) continue;
  double ngoodVtx=vtxList.size();
  AnaUtil::fillHist1D("nGoodvtx",ngoodVtx,1);
  if( ngoodVtx < 1 ) continue;
  if (useTrigger() && !isTriggered(true,false) ) continue;
  
     std::cout<<"Event="<<ev<<endl;

     muonSelector();
     electronSelector();
     pfphotonSelector(); 
//    if( !tightmuVec_.empty() ) {
    if( !looseMuPhoFSRpair_.empty() ) {
       for(auto& mu : looseMuPhoFSRpair_ ) {
        /*
         lepsyncDumper_ << "MuonInfoMiniAOD>>>" 
              << mu.first.pfChargedIsoR04 << ":" 
                               << mu.first.pfChargedIsoR04emu    << ":"
                               << mu.first.pfNeutralHadIsoR04    << ":"
                               << mu.first.pfPhotonIso04    << ":"
                               << mu.first.sumPUPt04    << endl;
         lepsyncDumper_ << "MuonInfoTreeVec>>>" 
              << mu.first.isolationMap["c40"].at(0) << ":" 
              << mu.first.isolationMap["c40"].at(1) << ":" 
              << mu.first.isolationMap["c40"].at(2) << ":" 
              << mu.first.isolationMap["c40"].at(3) << ":" 
              << mu.first.isolationMap["c40"].at(4) << endl;
       */
         if(mu.first.isPFMuon)
           dumpMuon(mu);
       }
     }  
    if( !looseElePhoFSRpair_.empty() ) {
     for(auto& e : looseElePhoFSRpair_ ) {
        double elePt = e.first.pt;
        double fSCeta = std::fabs(e.first.scEta);
        double BDT=e.first.BDT;
        bool isBDT = (elePt <= 10 && (( fSCeta < 0.8 && BDT > 0.47)   ||
                                      (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > 0.004) ||
                                      (fSCeta >= 1.479               && BDT > 0.295))) ||
        (elePt >  10 && ((fSCeta < 0.8 && BDT > -0.34)      ||
                         (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > -0.65) ||
                         (fSCeta >= 1.479               && BDT > 0.6)));
        if( isBDT ) 
          dumpElectron(e);
     }
   }
     if( evtype == EventType::mmmm ) {
         int size = looseMuPhoFSRpair_.size();
         AnaUtil::fillHist1D("nloosemu_ev4mu",size,1);
       if(size >= 4) AnaUtil::fillHist1D("formZZ_ev4mu",1,1);
     } else if( evtype == EventType::eeee ) {
         int size = looseElePhoFSRpair_.size();
         AnaUtil::fillHist1D("nlooseele_ev4e",size,1);
         if(size >= 4) AnaUtil::fillHist1D("formZZ_ev4e",1,1);
     } else if( evtype == EventType::eemm ) {
         int size1 = looseMuPhoFSRpair_.size();
         int size2 = looseElePhoFSRpair_.size();
         AnaUtil::fillHist1D("nloosemu_ev2e2mu",size1,1);
         AnaUtil::fillHist1D("nlooseele_ev2e2mu",size2,1);
         if(size1 + size2 >= 4) AnaUtil::fillHist1D("formZZ_ev2e2mu",1,1);
     }
     
 }  
 // Analysis is over
 endJob();
}

//
//loose muon////
void LeptonSync::muonSelector() {
    //std::cout<< AnaUtil::cutValue(muonCutMap(), "pt") << std::endl;
    std::cout<<"Muon Collection size="<<muonColl()->size()<<std::endl;
    for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
        const vhtm::Muon& muon = (*it);
        AnaUtil::fillHist1D("muCutFlow",0,1);
        
        //bool muType = !muon.isGlobalMuon && ( !muon.isTrackerMuon || !muon.isAllArbitrated );
        bool muType = muon.isGlobalMuon || ( muon.isTrackerMuon && muon.isAllArbitrated );
        if( muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt") )   continue;
        AnaUtil::fillHist1D("muCutFlow",1,1);
        if( std::fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta") )  continue;
        AnaUtil::fillHist1D("muCutFlow",2,1);
        AnaUtil::fillHist1D("mud0",muon.trkD0,1);
        AnaUtil::fillHist1D("mud0PV",muon.dxyPV,1);
        AnaUtil::fillHist1D("mudz",muon.trkDz,1);
        AnaUtil::fillHist1D("mudzPV",muon.dzPV,1);
        if( muon.dxyPV >= AnaUtil::cutValue(muonCutMap(), "dxyPV") )        continue;
        AnaUtil::fillHist1D("muCutFlow",3,1);
        if( std::fabs(muon.dzPV) >= AnaUtil::cutValue(muonCutMap(), "dzPV") ) continue;
        AnaUtil::fillHist1D("muCutFlow",4,1);
        if( !muType )                           continue;
        AnaUtil::fillHist1D("muCutFlow",5,1);
        if( std::fabs(muon.dB3D/muon.edB3D) > AnaUtil::cutValue(muonCutMap(), "SIP3D") )  continue;
        AnaUtil::fillHist1D("muCutFlow",6,1);
        loosemuVec_.push_back(muon);
        if( muon.isPFMuon ) {
            AnaUtil::fillHist1D("muCutFlow",7,1);
            tightmuVec_.push_back(muon);
        }
    }
    if( !loosemuVec_.empty() )
        std::sort(loosemuVec_.begin(),loosemuVec_.end(),PtComparator<vhtm::Muon>());
    if( !tightmuVec_.empty() )
        std::sort(tightmuVec_.begin(),tightmuVec_.end(),PtComparator<vhtm::Muon>());
    ///searchedMu_ = true;
    // std::cout << "POBJLoose/Tight mu=" << loosemuVec_.size() << "/" << tightmuVec_.size() << std::endl;
}
//////////////////electron selecton /////////////////
void LeptonSync::electronSelector() {
    for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it) {
        AnaUtil::fillHist1D("eleCutFlow",0,1);
        const vhtm::Electron& electron = (*it);
        double elePt = electron.pt;
        
        
        double eleta  = std::fabs(electron.eta);
        bool   etaCut = (eleta >= AnaUtil::cutValue(electronCutMap(), "etaLow") &&
                         eleta <= AnaUtil::cutValue(electronCutMap(), "etaUp")) ||
                         eleta >= AnaUtil::cutValue(electronCutMap(), "eta");
        
        if ( elePt <= AnaUtil::cutValue(electronCutMap(), "pt")  )          continue;
        AnaUtil::fillHist1D("eleCutFlow",1,1);
        if( etaCut )     continue;
        AnaUtil::fillHist1D("eleCutFlow",2,1);
        AnaUtil::fillHist1D("eled0",electron.trkD0,1);
        AnaUtil::fillHist1D("eled0PV",electron.dxyPV,1);
        AnaUtil::fillHist1D("eledz",electron.trkDz,1);
        AnaUtil::fillHist1D("eledzPV",electron.dzPV,1);
        if( electron.dxyPV >= AnaUtil::cutValue(electronCutMap(), "dxyPV") )      continue;
        AnaUtil::fillHist1D("eleCutFlow",3,1);
        if( std::fabs(electron.dzPV) >= AnaUtil::cutValue(electronCutMap(), "dzPV") )    continue;
        AnaUtil::fillHist1D("eleCutFlow",4,1);
        if( electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits") ) continue;
        AnaUtil::fillHist1D("eleCutFlow",5,1);
        if( std::fabs(electron.dB3D/electron.edB3D)  > AnaUtil::cutValue(electronCutMap(), "SIP3D") ) continue;
        AnaUtil::fillHist1D("eleCutFlow",6,1);
        
        looseeleVec_.push_back(electron);
        
        double fSCeta = std::fabs(electron.scEta);
        double BDT=electron.BDT;
        bool isBDT = (elePt <= 10 && (( fSCeta < 0.8 && BDT > 0.47)   ||
                                      (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > 0.004) ||
                                      (fSCeta >= 1.479               && BDT > 0.295))) ||
        (elePt >  10 && ((fSCeta < 0.8 && BDT > -0.34)      ||
                         (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > -0.65) ||
                         (fSCeta >= 1.479               && BDT > 0.6)));
        if( isBDT ) {
            AnaUtil::fillHist1D("eleCutFlow",7,1);
            tighteleVec_.push_back(electron);
        }
        
    }
    if( !looseeleVec_.empty() )
        std::sort(looseeleVec_.begin(),looseeleVec_.end(),PtComparator<vhtm::Electron>());
    if( !tighteleVec_.empty() )
        std::sort(tighteleVec_.begin(),tighteleVec_.end(),PtComparator<vhtm::Electron>());
    //searchedEle_ = true;
    //std::cout << "POBJLoose/Tight e=" << looseeleVec_.size() << "/" << tighteleVec_.size() << std::endl;
}
//packed pf photon///////
void LeptonSync::pfphotonSelector() {
    // if( !searchedEle_ ) {
    ///    electronSelector();//std::cout<<" POBJ point 1>>";
    //  }
    //  if( !searchedMu_ ){
    // std::cout<<" POBJ point 2>>";
    //     muonSelector();
    // }
    //std::cout<<"mu size>>"<<loosemuVec_.size()<<">>>>>elesize>>>"<<looseeleVec_.size()<<std::endl;
    //Pre-selection of PF Photons
    std::vector<TLorentzVector> eleP4vec;
    if( !looseeleVec_.empty() ) {
        for( auto& e : looseeleVec_ ) {
            TLorentzVector eP4;
            eP4.SetPtEtaPhiE(e.pt,e.eta,e.phi,e.energy);
            eleP4vec.push_back(eP4);
        }
    }
    //std::cout<<" POBJ point 3>>";
    std::vector<vhtm::PackedPFCandidate> pfphotonVectemp;//preselected PFPhotons
    //std::cout<<"ppfsize>>>>"<<packedPFCandidateColl()->size()<<endl;
    for (auto it = packedPFCandidateColl()->begin(); it != packedPFCandidateColl()->end(); ++it) {
        const vhtm::PackedPFCandidate& pfphoton = (*it);
        // std::cout<<" POBJ point 3a>>";
        if( pfphoton.pdgId != 22 || pfphoton.pt < 2. || std::abs(pfphoton.eta) > 2.4 )   continue;
        TLorentzVector pfPhoP4;
        pfPhoP4.SetPtEtaPhiE(pfphoton.pt,pfphoton.eta,pfphoton.phi,pfphoton.energy);
        bool eleVeto = true;
        //std::cout<<" POBJ point 3b>>";
        for( auto& eP4 : eleP4vec ) {
            if( eP4.DeltaR(pfPhoP4) < 0.15 ) {
                eleVeto = false;
                break;
            }
        }
        //std::cout<<" POBJ point 3c>>";
        if( !eleVeto )   continue;
        pfphotonVectemp.push_back( pfphoton );//photon preselection
    }
    //std::cout<<" POBJ point 4>>";
    //Now select FSR photons for each lepton
    //Note->Multiple photons can be attached to a lepton
    //Screening for those will be done later down the analysis chain
    //std::cout<<"POBJ loosemuVec_ size="<<loosemuVec_.size() << "\n";
    //std::cout<<"POBJ pfphotonVectemp size="<<pfphotonVectemp.size() << "\n";
    for( auto& mu : loosemuVec_ ) {
        TLorentzVector muP4;
        muP4.SetPtEtaPhiE(mu.pt,mu.eta,mu.phi,mu.energy);
        std::vector<vhtm::PackedPFCandidate> photemp;
        for( auto& pfPho : pfphotonVectemp ) {
            TLorentzVector phoP4;
            phoP4.SetPtEtaPhiE(pfPho.pt,pfPho.eta,pfPho.phi,pfPho.energy);
            if( muP4.DeltaR(phoP4) < 0.07 && pfPho.pt > 2. ) {
                photemp.push_back(pfPho);
                pfphotonVec_.push_back(pfPho);
            } else if( muP4.DeltaR(phoP4) >= 0.07 && muP4.DeltaR(phoP4) < 0.5 && pfPho.pt > 4. ) {
                //double PFdbiso = ch + std::max(0.0, nu-0.5*pu);
                double PFdbiso = pfPho.isolationMap.at("c30").at(0) +
                std::max(0.0, pfPho.isolationMap.at("c30").at(1) -
                         0.5*pfPho.isolationMap.at("c30").at(2));
                if( PFdbiso/pfPho.pt < 1.) {
                    photemp.push_back(pfPho);
                    pfphotonVec_.push_back(pfPho);
                }
            }
        }


        std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > mufsrtemp = std::make_pair(mu,photemp);
        looseMuPhoFSRpair_.push_back( mufsrtemp );
    }
    std::cout<<"POBJ loosemuVec//looseMuPhoFSRpair size="<<loosemuVec_.size()<<"//"<<looseMuPhoFSRpair_.size() << "\n";
    //For electron
    for( auto& ele : looseeleVec_ ) {
        TLorentzVector eleP4;
        eleP4.SetPtEtaPhiE(ele.pt,ele.eta,ele.phi,ele.energy);
        std::vector<vhtm::PackedPFCandidate> photemp;
        for( auto& pfPho : pfphotonVectemp ) {
            TLorentzVector phoP4;
            phoP4.SetPtEtaPhiE(pfPho.pt,pfPho.eta,pfPho.phi,pfPho.energy);
            if( eleP4.DeltaR(phoP4) < 0.07 && pfPho.pt > 2. ) {
                photemp.push_back(pfPho);
                pfphotonVec_.push_back(pfPho);
            } else if( eleP4.DeltaR(phoP4) >= 0.07 && eleP4.DeltaR(phoP4) < 0.5 && pfPho.pt > 4. ) {
                //double PFdbiso = ch + std::max(0.0, nu-0.5*pu);
                double PFdbiso = pfPho.isolationMap.at("c30").at(0) +
                std::max(0.0, pfPho.isolationMap.at("c30").at(1) -
                         0.5*pfPho.isolationMap.at("c30").at(2));
                if( PFdbiso/pfPho.pt < 1.)  {
                    photemp.push_back(pfPho);
                    pfphotonVec_.push_back(pfPho);
                }
            }
        }
        looseElePhoFSRpair_.push_back( std::make_pair(ele,photemp) );
    }
//    std::cout<<"POBJ looseeleVec//looseElePhoFSRpair size="<<looseeleVec_.size()<<"//"<<looseElePhoFSRpair_.size() << "\n";
    
    
    // }
}
//////////////////////////////////
//Gen Level Functions
//////////////////////////////////
template <typename T>
TLorentzVector LeptonSync::getP4(const T& obj) {
   TLorentzVector l1;
   l1.SetPtEtaPhiE(obj.pt,obj.eta,obj.phi,obj.energy);
   return l1;
}

int LeptonSync::getGenDauPgd(const GenParticle gp) {
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

double LeptonSync::getHmassfromZdau(const GenParticle Z1,const GenParticle Z2) {
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
bool LeptonSync::genOk() {
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
           // AnaUtil::fillHist1D("genHmass",l.M(),1);
          }
//          if (pdgid == 23 && status==2 ) {
          else if (pdgid == 23 && status == 3 ) {
            int mmid = -1;
            int index = getMotherId(gp, mmid);
           // AnaUtil::fillHist1D("genZmother",mmid,1);
            if (index < 0 || fabs(mmid) != 25) continue; // H decay to Higgs
            nHtoZ++;
            genZ.push_back(gp);
            TLorentzVector l1;
            l1.SetPtEtaPhiE(gp.pt,gp.eta,gp.phi,gp.energy);
           // AnaUtil::fillHist1D("genZmass",l1.M(),1);
           //check e mu decay of Z
           vector<int> d = gp.daughtIndices;
           int ndau = 0;
           for (auto it = d.begin(); it != d.end(); ++it) {
             int di = (*it);
             if (di >= ngenparticle()) continue;
             const GenParticle& dgp = genParticleColl()->at(di);
             int pid = fabs(dgp.pdgId);
            // if (pid == 13) AnaUtil::fillHist1D("nZdauflav",1,1); 
            // else if (pid == 11) AnaUtil::fillHist1D("nZdauflav",3,1);
           //  else if (pid == 15) AnaUtil::fillHist1D("nZdauflav",5,1);
            // else AnaUtil::fillHist1D("nZdauflav",7,1);
           }                
          }
        }//gen particle loop
          std::sort(genZ.begin(),genZ.end(),MassComparator<vhtm::GenParticle>());
         // AnaUtil::fillHist1D("genZ1mass",getP4(genZ[0]).M(),1); 
         // AnaUtil::fillHist1D("genZ2mass",getP4(genZ[1]).M(),1);
         // AnaUtil::fillHist1D("genHmassfromZ",(getP4(genZ[0])+getP4(genZ[1])).M(),1);
          //Event Flavour
          if( getGenDauPgd(genZ[0]) == 13 && getGenDauPgd(genZ[1]) == 13 ) {
         ////   AnaUtil::fillHist1D("nEvntflav",1,1);
          //  AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
            genPass=true;
            evtype=EventType::mmmm;
          }
          else if( getGenDauPgd(genZ[0]) == 11 && getGenDauPgd(genZ[1]) == 11 ) {
           // AnaUtil::fillHist1D("nEvntflav",3,1);
           // AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
            genPass=true;
            evtype=EventType::eeee;
          }
          else if( (getGenDauPgd(genZ[0]) == 13 && getGenDauPgd(genZ[1]) == 11) || 
                   (getGenDauPgd(genZ[0]) == 11 && getGenDauPgd(genZ[1]) == 13 ) ) {
           // AnaUtil::fillHist1D("nEvntflav",5,1);
           // AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
            genPass=true;
            evtype=EventType::eemm;
          } else {
           // AnaUtil::fillHist1D("nEvntflav",7,1);
          //AnaUtil::fillHist1D("genHmassfromZdau",getHmassfromZdau(genZ[0],genZ[1]),1);
          }
       if(genPass) {   
         //std::cout<<"Event Type="<<evtype<<std::endl;
         //Z1 daughters
          /*
          if (getGenDauPgd(genZ[0]) == 13 ) AnaUtil::fillHist1D("nZ1dauflav",1,1);
          else if (getGenDauPgd(genZ[0]) == 11) AnaUtil::fillHist1D("nZ1dauflav",3,1);
          else if (getGenDauPgd(genZ[0]) == 15) AnaUtil::fillHist1D("nZ1dauflav",5,1);
          else AnaUtil::fillHist1D("nZ1dauflav",7,1);
          //Z2 daughter
          if (getGenDauPgd(genZ[1]) == 13 ) AnaUtil::fillHist1D("nZ2dauflav",1,1);
          else if (getGenDauPgd(genZ[1]) == 11) AnaUtil::fillHist1D("nZ2dauflav",3,1);
          else if (getGenDauPgd(genZ[1]) == 15) AnaUtil::fillHist1D("nZ2dauflav",5,1);
          else AnaUtil::fillHist1D("nZ2dauflav",7,1);
          //AnaUtil::fillHist1D("nEvntflavvv",evtype,1);
          */
       }
   return genPass;
}
//////////////////////Gen Level/////////////////////
//dumpers for sync exercise
void LeptonSync::dumpMuon(std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > m) {
  int pdg =  ( m.first.charge == -1 ) ? 13 : -13;
  lepsyncDumper_ << run_ << ":"
                 << lumi_ << ":" 
                 << event_ << ":"
                 << pdg << ":"
                 << std::fixed << std::setprecision(2) 
                 << m.first.pt << ":" 
                 << m.first.eta << ":" 
                 << m.first.phi << ":"
                 << std::fabs(m.first.dB3D/m.first.edB3D) << ":"
                 << m.first.pfChargedIsoR04 << ":" 
                 << m.first.pfNeutralHadIsoR04 << ":"
                 << m.first.pfPhotonIso04 << ":"
                 << m.first.sumPUPt04 << ":" 
                 << std::fixed << std::setprecision(3) 
                 << m.first.pfRelIso04 << ":" 
                 << 0. << ":";
              if( !m.second.empty() ){
                 std::sort(m.second.begin(),m.second.end(),PtComparator<vhtm::PackedPFCandidate>());
                 TLorentzVector muP4,phoP4;
                 muP4.SetPtEtaPhiE(m.first.pt,m.first.eta,m.first.phi,m.first.energy);
                 phoP4.SetPtEtaPhiE( m.second.at(0).pt,
                                       m.second.at(0).eta,
                                       m.second.at(0).phi,
                                       m.second.at(0).energy );
                 lepsyncDumper_ << std::fixed << std::setprecision(2) 
                           << phoP4.Pt() << ":"
                           << phoP4.DeltaR(muP4) 
         		   //<< ":" 
                           //<< 0. 
                           << endl;
               } else {
                lepsyncDumper_  << 0. << ":"
                 << 0.
                // << ":" 
                // << 0. 
                 << endl;  
               }
      
}

void LeptonSync::dumpElectron( std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > e ) {
  int pdg =  ( e.first.charge == -1 ) ? 11 : -11;
  double relIso = calcRhoLepIso(e.first,rho_)/e.first.pt;
  lepsyncDumper_ << run_ << ":"
                 << lumi_ << ":" 
                 << event_ << ":"
                 << pdg << ":" 
                 << std::fixed << std::setprecision(2) 
                 << e.first.pt << ":" 
                 << e.first.eta << ":" 
                 << e.first.phi << ":"
                 << std::fabs(e.first.dB3D/e.first.edB3D) << ":"
                 << e.first.chargedHadronIso << ":" 
                 << e.first.neutralHadronIso << ":"
                 << e.first.photonIso << ":"
                 << rho_ << ":"
                 << std::fixed << std::setprecision(3) 
                 << relIso << ":"
                 << e.first.BDT << ":";
              if( !e.second.empty() ){
                 std::sort(e.second.begin(),e.second.end(),PtComparator<vhtm::PackedPFCandidate>());
                 TLorentzVector muP4,phoP4;
                 muP4.SetPtEtaPhiE(e.first.pt,e.first.eta,e.first.phi,e.first.energy);
                 phoP4.SetPtEtaPhiE(   e.second.at(0).pt,
                                       e.second.at(0).eta,
                                       e.second.at(0).phi,
                                       e.second.at(0).energy );
                 lepsyncDumper_ << std::fixed << std::setprecision(2) 
                           << phoP4.Pt() << ":"
                           << phoP4.DeltaR(muP4) 
                           //<< ":" 
                           //<< 0. 
                           << endl;
               } else {
                lepsyncDumper_  << 0. << ":"
                 << 0. 
                //<< ":" 
               //  << 0. 
                 << endl;  
               }

}

double LeptonSync::calcRhoLepIso(vhtm::Electron lep,double rho) {
   // double iso = lep.isolationMap.at("c40").at(0) + std::max(0.0, lep.isolationMap.at("c40").at(1) - phoEnergy
   //- effArea*rho);
    double effArea = getEleRhoEffectiveArea( std::fabs(lep.scEta));  
    double iso = lep.chargedHadronIso + std::max(0.,lep.neutralHadronIso + lep.photonIso - effArea*rho );
   // double iso = lep.chargedHadronIso + std::max(0.,lep.neutralHadronIso - effArea*rho );
    return iso;
}

double LeptonSync::getEleRhoEffectiveArea( double fSCeta) {
    if (fSCeta >= 0.0 && fSCeta < 0.8 ) return 0.1830;
    else if (fSCeta >= 0.8 && fSCeta < 1.3 ) return 0.1734;
    else if (fSCeta >= 1.3 && fSCeta < 2.0 ) return 0.1077;
    else if (fSCeta >= 2.0 && fSCeta < 2.2 ) return 0.1565;
    else if (fSCeta >= 2.2) return 0.2680;
}


void LeptonSync::endJob() {
 lepsyncDumper_.close();
 
 closeFiles();

 histf()->cd();
 histf()->Write();
 histf()->Close();
 delete histf();
}
