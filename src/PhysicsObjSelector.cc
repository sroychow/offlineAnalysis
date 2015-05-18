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
#include "PhysicsObjSelector.h"


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
PhysicsObjSelector::PhysicsObjSelector()
: AnaBase()
{
 searchedEle_ = false; 
 searchedMu_ = false; 
 searchedPhoton_ = false;
    
}
void PhysicsObjSelector::bookHistograms() {
  cout<<"ddd";
}
//loose muon////
void PhysicsObjSelector::muonSelector() {
    //std::cout<< AnaUtil::cutValue(muonCutMap(), "pt") << std::endl;
 for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const vhtm::Muon& muon = (*it);
    //bool muType = !muon.isGlobalMuon && ( !muon.isTrackerMuon || !muon.isAllArbitrated );
    bool muType = muon.isGlobalMuon || ( muon.isTrackerMuon && muon.isAllArbitrated );
    if ( muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")                            ||
         std::fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta")               ||
         muon.dxyPV >= AnaUtil::cutValue(muonCutMap(), "dxyPV")                      ||
         //muon.trkD0 >= AnaUtil::cutValue(muonCutMap(), "dxyPV")                      ||
         std::fabs(muon.dzPV) >= AnaUtil::cutValue(muonCutMap(), "dzPV")             ||
         //std::fabs(muon.trkDz) >= AnaUtil::cutValue(muonCutMap(), "dzPV")             ||
         !muType                                                                    ||
         std::fabs(muon.dB3D/muon.edB3D) >= AnaUtil::cutValue(muonCutMap(), "SIP3D")
       )  continue;
    loosemuVec_.push_back(muon);
    if( muon.isPFMuon ) tightmuVec_.push_back(muon);
 }
 if( !loosemuVec_.empty() )           
   std::sort(loosemuVec_.begin(),loosemuVec_.end(),PtComparator<vhtm::Muon>());
 if( !tightmuVec_.empty() )           
   std::sort(tightmuVec_.begin(),tightmuVec_.end(),PtComparator<vhtm::Muon>());
 searchedMu_ = true;
   // std::cout << "POBJLoose/Tight mu=" << loosemuVec_.size() << "/" << tightmuVec_.size() << std::endl;
}
//////////////////electron selecton /////////////////
void PhysicsObjSelector::electronSelector() {
  for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it) {
    const vhtm::Electron& electron = (*it); 
    double elePt = electron.pt;
    if (elePt <= AnaUtil::cutValue(electronCutMap(), "pt")) continue;
    double eleta  = std::fabs(electron.eta);
    //bool   etaCut = (eleta >= AnaUtil::cutValue(electronCutMap(), "etaLow") &&
   //	             eleta <= AnaUtil::cutValue(electronCutMap(), "etaUp")) ||
     bool etaCut =  eleta  >= AnaUtil::cutValue(electronCutMap(), "eta");
    if ( elePt <= AnaUtil::cutValue(electronCutMap(), "pt")                         ||
         etaCut                                                                     ||
         electron.dxyPV >= AnaUtil::cutValue(electronCutMap(), "dxyPV")              ||
         //electron.trkD0 >= AnaUtil::cutValue(electronCutMap(), "dxyPV")              ||
         std::fabs(electron.dzPV) >= AnaUtil::cutValue(electronCutMap(), "dzPV")                ||
         //std::fabs(electron.trkDz) >= AnaUtil::cutValue(electronCutMap(), "dzPV")                ||
         electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits")  ||
         std::fabs(electron.dB3D/electron.edB3D)  >= AnaUtil::cutValue(electronCutMap(), "SIP3D")      
       ) continue;
      
      looseeleVec_.push_back(electron);

      double fSCeta = std::fabs(electron.scEta);
      double BDT=electron.BDT;
      bool isBDT = (elePt <= 10 && (( fSCeta < 0.8 && BDT > 0.47)   ||
                   (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > 0.004) ||
                   (fSCeta >= 1.479               && BDT > 0.295))) || 
                   (elePt >  10 && ((fSCeta < 0.8 && BDT > -0.34)      ||
                   (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > -0.65) || 
                   (fSCeta >= 1.479               && BDT > 0.6)));
      if( isBDT ) tighteleVec_.push_back(electron); 
 
  }
 if( !looseeleVec_.empty() ) 
   std::sort(looseeleVec_.begin(),looseeleVec_.end(),PtComparator<vhtm::Electron>());
 if( !tighteleVec_.empty() ) 
   std::sort(tighteleVec_.begin(),tighteleVec_.end(),PtComparator<vhtm::Electron>());
 searchedEle_ = true;
    //std::cout << "POBJLoose/Tight e=" << looseeleVec_.size() << "/" << tighteleVec_.size() << std::endl;
}
//packed pf photon///////
void PhysicsObjSelector::pfphotonSelector() {
    if( !searchedEle_ ) {
        electronSelector();//std::cout<<" POBJ point 1>>";
    }
    if( !searchedMu_ ){
       // std::cout<<" POBJ point 2>>";
        muonSelector();
    }
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
       // if ( photemp.empty() ) {
       //     vhtm::PackedPFCandidate pho;
       //     photemp.push_back(pho);
       // }
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
std::cout<<"POBJ looseeleVec//looseElePhoFSRpair size="<<looseeleVec_.size()<<"//"<<looseElePhoFSRpair_.size() << "\n";

 
 // }
}
////////////////////clear lists/////////////////
void PhysicsObjSelector::clear() {
  clearEvent();  
  loosemuVec_.clear();
  tightmuVec_.clear();
  looseeleVec_.clear();
  tighteleVec_.clear();
  pfphotonVec_.clear();
  looseMuPhoFSRpair_.clear();
  looseElePhoFSRpair_.clear();
  searchedEle_ = false; 
  searchedMu_ = false; 
  searchedPhoton_ = false;
 }















