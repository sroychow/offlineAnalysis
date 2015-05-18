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
using std::setprecision;
using std::setw;

const double PhysicsObjSelector::MZnominal = 91.1876;

PhysicsObjSelector::PhysicsObjSelector()
  : AnaBase(),
    searchedEle_(false),
    searchedMu_(false),
    searchedPhoton_(false),
    nJets_(0),
    nbJets_(0)
{
}
bool PhysicsObjSelector::beginJob() {
  AnaBase::beginJob();
  histf()->cd();
  histf()->mkdir("ObjectSelection");
  return true;
}
void PhysicsObjSelector::bookHistograms() {
  histf()->cd();
  histf()->cd("ObjectSelection");
  new TH1D("muCutFlow", "Muon Cut Flow", 10, -0.5, 10.5);
  new TH1D("eleCutFlow", "Electron Cut Flow", 9, -0.5, 8.5);
  new TH1D("jetCutFlow", "Jet Cut Flow", 9, -0.5, 8.5); 
  new TH1D("nJets", "Number of jets cleaned w.r.t tight leptons passing iso per event", 20, 0, 20);
  new TH1D("btagfactor", "bDiscriminator", 50, 0., 2.);
  new TH1D("nbJets", "Number of b-jets per event", 20, 0, 20);
}
// muon selection
void PhysicsObjSelector::muonSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");

  for (const auto& muon: *muonColl()) {
    AnaUtil::fillHist1D("muCutFlow", 0, wt);

    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt"))                             continue;
    AnaUtil::fillHist1D("muCutFlow", 1, wt);

    if (std::fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta"))                continue;
    AnaUtil::fillHist1D("muCutFlow", 2, wt);

    if (std::fabs(muon.dxyPV) >= AnaUtil::cutValue(muonCutMap(), "dxyPV") )           continue;
    AnaUtil::fillHist1D("muCutFlow", 3, wt);

    if (std::fabs(muon.dzPV) >= AnaUtil::cutValue(muonCutMap(), "dzPV") )             continue;
    AnaUtil::fillHist1D("muCutFlow", 4, wt);

    bool muType = muon.isGlobalMuon || (muon.isTrackerMuon && muon.matches > 0);
    if (!muType)                                                                      continue;
    AnaUtil::fillHist1D("muCutFlow", 5, wt);

    if (muon.muonBestTrackType == 2)                                                  continue;
    AnaUtil::fillHist1D("muCutFlow", 6, wt);

    if (std::fabs(muon.dB3D/muon.edB3D) >= AnaUtil::cutValue(muonCutMap(), "SIP3D"))  continue;
    AnaUtil::fillHist1D("muCutFlow", 7, wt);

    if (!muon.isghostCleaned)                                                         continue;
    AnaUtil::fillHist1D("muCutFlow", 8, wt);

    // loose muon 
    loosemuVec_.push_back(muon);

    // attach an empty PackedPFCandidate vector to each loose Muon
    std::vector<vhtm::PackedPFCandidate> phov;
    looseMuPhotonPairVec_.push_back({muon, phov});

    // tight muon 
    if (muon.isPFMuon) {
      AnaUtil::fillHist1D("muCutFlow", 9, wt);
      tightmuVec_.push_back(muon);
    }
  }
  if (loosemuVec_.size() > 1)
    std::sort(loosemuVec_.begin(), loosemuVec_.end(), PtComparator<vhtm::Muon>());
  if (tightmuVec_.size() > 1)           
    std::sort(tightmuVec_.begin(), tightmuVec_.end(), PtComparator<vhtm::Muon>());

  searchedMu_ = true;
}
// electron selecton
void PhysicsObjSelector::electronSelector(double wt) {
  histf()->cd();
  histf()->cd("ObjectSelection");

  if (!searchedMu_) muonSelector();

  for (const auto& electron: *electronColl()) {
    AnaUtil::fillHist1D("eleCutFlow", 0, wt);
    
    if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt"))                                continue;
    AnaUtil::fillHist1D("eleCutFlow", 1, wt);

    if (std::fabs(electron.eta) >= AnaUtil::cutValue(electronCutMap(), "eta"))                   continue;
    AnaUtil::fillHist1D("eleCutFlow", 2, wt);

    if (std::fabs(electron.dxyPV) >= AnaUtil::cutValue(electronCutMap(), "dxyPV"))               continue;
    AnaUtil::fillHist1D("eleCutFlow", 3, wt);

    if (std::fabs(electron.dzPV) >= AnaUtil::cutValue(electronCutMap(), "dzPV"))                 continue;
    AnaUtil::fillHist1D("eleCutFlow", 4, wt);

    if (electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits"))               continue;
    AnaUtil::fillHist1D("eleCutFlow", 5, wt);

    if (std::fabs(electron.dB3D/electron.edB3D) >= AnaUtil::cutValue(electronCutMap(), "SIP3D")) continue;
    AnaUtil::fillHist1D("eleCutFlow", 6, wt);
    
    // Cross cleaning
    if (!crossCleaned(electron))  continue;
    AnaUtil::fillHist1D("eleCutFlow", 7, wt);

    // loose electrons
    looseeleVec_.push_back(electron);
    
    // attach an empty PackedPFCandidate vector to each loose electron
    std::vector<vhtm::PackedPFCandidate> phov;
    looseElePhotonPairVec_.push_back({electron, phov});

    // tight electrons
    if (electronBDT(electron)) {
      AnaUtil::fillHist1D("eleCutFlow", 8, wt);
      tighteleVec_.push_back(electron);
    }
  }
  if (looseeleVec_.size() > 1)
    std::sort(looseeleVec_.begin(), looseeleVec_.end(), PtComparator<vhtm::Electron>());
  if (tighteleVec_.size() > 1) 
    std::sort(tighteleVec_.begin(), tighteleVec_.end(), PtComparator<vhtm::Electron>());

  searchedEle_ = true;
}
// photon from PFCandidates
void PhysicsObjSelector::photonSelector() {
  if (!searchedMu_)  muonSelector();
  if (!searchedEle_) electronSelector();

  // -- from twiki --
  // FSR Photon selection: for each photon, consider the closest lepton:
  // Note: Multiple photons can be attached to a lepton
  // Muons

  // Pre-selection of PF Photons
  if (looseeleVec_.size() || loosemuVec_.size()) {
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pdgId != 22 || pfcand.pt <= 2. || std::fabs(pfcand.eta) >= 2.4) continue;
      
      // Supercluster veto
      if (!passedSuperClusterVeto(pfcand)) continue;
      
      int muindx = -1, elindx = -1;
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      bool leptonMatched = false;
      if (dRmin < 0.07) {
	if (pfcand.pt > 2.) leptonMatched = true;
      }
      else if (dRmin < 0.5) {
	if (pfcand.pt > 4.) {
	  double iso = pfcand.isolationMap.at("c30").at(0)
	    + pfcand.isolationMap.at("c30").at(2) 
	    + pfcand.isolationMap.at("c30").at(3) 
	    + pfcand.isolationMap.at("c30").at(4);
	  if (iso/pfcand.pt < 1.) leptonMatched = true;
	}
      } 
      if (!leptonMatched) continue;
      photonVec_.push_back(pfcand);
      
      // Now update vector of pair
      if (elindx > -1)
	looseElePhotonPairVec_.at(elindx).second.push_back(pfcand);
      else if (muindx > -1)
	looseMuPhotonPairVec_.at(muindx).second.push_back(pfcand);
    }
  }
  // now fill the vectors with tight leptons
  for (const auto& el: looseElePhotonPairVec_)
    if (electronBDT(el.first)) tightElePhotonPairVec_.push_back(el);

  for (const auto& el: looseMuPhotonPairVec_)
    if (el.first.isPFMuon) tightMuPhotonPairVec_.push_back(el);
}
void PhysicsObjSelector::leptonCrossCleaning() {
  vector<pair<vhtm::Electron, vector<vhtm::PackedPFCandidate> > > list;
  for (const auto& it: tightElePhotonPairVec_) 
    if (crossCleaned(it.first)) list.push_back(it);

  tightElePhotonPairVec_.clear();
  tightElePhotonPairVec_ = list;
}
bool PhysicsObjSelector::crossCleaned(const vhtm::Electron& electron) const {
  bool flag = true;
  for (const auto& mu: tightmuVec_) {
    if (getP4(electron).DeltaR(getP4(mu)) < 0.05) {
      flag = false;
      break;
    }
  }
  return flag;
}


void PhysicsObjSelector::jetSelector(double wt) {
  for (const auto& jet: *jetColl()) {
    AnaUtil::fillHist1D("jetCutFlow", 0, wt);
    if( jet.pt <= AnaUtil::cutValue(jetCutMap(), "pt") )                   continue;
    AnaUtil::fillHist1D("jetCutFlow", 1, wt);
    if( std::fabs(jet.eta) >= AnaUtil::cutValue(jetCutMap(), "eta") )      continue;
    AnaUtil::fillHist1D("jetCutFlow", 2, wt);
    if( !jetpuMVAid(jet) )                                                continue; 
    AnaUtil::fillHist1D("jetCutFlow", 3, wt);
    if( !jetLeptoncleaning( jet,AnaUtil::cutValue(jetCutMap(), "dRlep")) ) continue;
    AnaUtil::fillHist1D("jetCutFlow", 4, wt);
    AnaUtil::fillHist1D("btagfactor", jet.combinedInclusiveSecondaryVertexV2BJetTags, wt);
    if( jet.combinedInclusiveSecondaryVertexV2BJetTags > AnaUtil::cutValue(jetCutMap(), "btagFactor") )
      nbJets_++;
    jetVec_.push_back(jet);
  } 
  nJets_= jetVec_.size();
  if ( nJets_ > 1)
    std::sort(jetVec_.begin(), jetVec_.end(), PtComparator<vhtm::Jet>());
  AnaUtil::fillHist1D("nJets",nJets_,wt);
  AnaUtil::fillHist1D("nbJets",nbJets_,wt);
}


bool PhysicsObjSelector::jetpuMVAid(const vhtm::Jet& jet) {
  float jpumva = jet.jpumva;
  double jpt = jet.pt;
  double jeta = std::fabs(jet.eta);
  bool passPU = true;
  if(jpt>20){
    if(jeta>3.){
      if(jpumva<=-0.45) passPU=false;
    }else if(jeta>2.75){
      if(jpumva<=-0.55) passPU=false;
    }else if(jeta>2.5){
      if(jpumva<=-0.6) passPU=false;
    }else if(jpumva<=-0.63) passPU=false;
  }else{
    if(jeta>3.){
      if(jpumva<=-0.95) passPU=false;
    }else if(jeta>2.75){
      if(jpumva<=-0.94) passPU=false;
    }else if(jeta>2.5){
      if(jpumva<=-0.96) passPU=false;
    }else if(jpumva<=-0.95) passPU=false;
  }
  return passPU;
}

bool PhysicsObjSelector::jetLeptoncleaning( const vhtm::Jet& jet,const float dR){
  //bool jetCleaned = true;
  TLorentzVector jetP4 = getP4(jet);
  for( auto& mu: tightmuVec_ ) {
    if( ( mu.sumChargedHadronPt + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt - 0.5 * mu.sumPUPt))/mu.pt >= 0.4 )
      continue;
    if( jetP4.DeltaR(getP4(mu)) <= dR ) {
      //jetCleaned = false;
      //break;
      return false;
    }     
  }  

  for( auto& ele: tighteleVec_ ) {
    if( ( ele.chargedHadronIso + 
          std::max(0., ele.neutralHadronIso + ele.photonIso  -
                       getEleRhoEffectiveArea(std::fabs(ele.eta)) * fGridRhoFastjetAll_))/ele.pt >= 0.5)  continue;
    if( jetP4.DeltaR(getP4(ele)) <= dR ) {
      //jetCleaned = false;
      //break;
      return false;
    }     
  }  
 
  return true;
}

// clear lists
void PhysicsObjSelector::clear() {
  loosemuVec_.clear();
  tightmuVec_.clear();
  looseeleVec_.clear();
  tighteleVec_.clear();
  photonVec_.clear();
  jetVec_.clear();
  looseElePhotonPairVec_.clear();
  tightElePhotonPairVec_.clear();
  looseMuPhotonPairVec_.clear();
  tightMuPhotonPairVec_.clear();

  searchedMu_ = false; 
  searchedEle_ = false; 
  searchedPhoton_ = false;
  nJets_ = 0;
  nbJets_ = 0;
}

bool PhysicsObjSelector::electronBDT(const vhtm::Electron& electron) {
  double scEta = std::fabs(electron.scEta);
  double elePt = electron.pt;
  double BDT = electron.BDT;
  bool isBDT = (elePt <= 10 && ((scEta < 0.8 && BDT > -0.202)                     ||
				((scEta >= 0.8 && scEta < 1.479) && BDT > -0.444) ||
				(scEta >= 1.479 && BDT > 0.264)))                 ||
               (elePt >  10 && ((scEta < 0.8 && BDT > -0.110)                     ||
				((scEta >= 0.8 && scEta < 1.479) && BDT > -0.284) ||
				(scEta >= 1.479 && BDT > -0.212)));
  return isBDT;
}
bool PhysicsObjSelector::passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose) const {
  // Supercluster veto: remove all PF photons that match with any electron passing loose ID and SIP cuts; 
  // matching is according to (|deta| < 2, |dphi| < 0.05) OR (dR < 0.15), with respect to the electron's supercluster. 
  bool passedVeto = true;
  if (verbose)
  cout << "    pfPt   pfEta   pfPhi   elePt   scEta  elePhi    dEta    dPhi      dR" << endl;
  for (const auto& ele: looseeleVec_) {
    double deta = ele.scEta - pfcand.eta;
    double dphi = TVector2::Phi_mpi_pi(ele.scPhi - pfcand.phi);
    if (dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();  
    float dR = std::sqrt(deta * deta + dphi * dphi);
    if (verbose)
    cout << setw(8) << pfcand.pt
         << setw(8) << pfcand.eta
         << setw(8) << pfcand.phi
         << setw(8) << ele.pt
         << setw(8) << ele.scEta
         << setw(8) << ele.scPhi
         << setw(8) << deta
         << setw(8) << dphi
         << setw(8) << dR
         << endl;
    if ((std::fabs(deta) < 0.05 && std::fabs(dphi) < 2.0) || dR < 0.15) {
      passedVeto = false;
      break;
    }
  }
  return passedVeto;
}
double PhysicsObjSelector::findClosestLepton(const vhtm::PackedPFCandidate& pfPho, int& muindx, int& elindx) const {
  TLorentzVector phoP4 = getP4(pfPho);
  double dRmin = 999;
  muindx = -1;
  // First consider loose muons
  for (unsigned int i = 0; i < loosemuVec_.size(); ++i) {
    const vhtm::Muon& mu = loosemuVec_[i];
    TLorentzVector muP4 = getP4(mu);
    double dR = muP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      muindx = i;
    }
  }
  // Then consider loose electron
  elindx = -1;
  for (unsigned int i = 0; i < looseeleVec_.size(); ++i) {
    const vhtm::Electron& ele = looseeleVec_[i];
    TLorentzVector eleP4 = getP4(ele);
    double dR = eleP4.DeltaR(phoP4);
    if (dR < dRmin) {
      dRmin = dR;
      elindx = i;
    }
  }
  return dRmin;
}
void PhysicsObjSelector::findObjects(double evRho,double wt) {
  setEventGridRho(evRho);
  // muonSelector must precede electronSelector
  muonSelector(wt);
  electronSelector(wt);
  
  // after electrons and muons are found, find photons
  photonSelector();
  jetSelector(wt);
}
void PhysicsObjSelector::printZCandidate(const ZCandidate& Za, const std::string& tag, std::ostream& os) {
  os << std::setprecision(3);
  os << tag << endl;
  os << " -- Leptons: " << endl;
  os << "  indx      pT     eta     phi  energy  charge  reliso"
     << endl; 
  os << setw(6) << 1
     << setw(8) << Za.l1P4.Pt()
     << setw(8) << Za.l1P4.Eta()
     << setw(8) << Za.l1P4.Phi()
     << setw(8) << Za.l1P4.Energy()
     << setw(8) << Za.l1Charge
     << setw(8) << Za.l1Isolation
     << endl;
  os << setw(6) << 2
     << setw(8) << Za.l2P4.Pt()
     << setw(8) << Za.l2P4.Eta()
     << setw(8) << Za.l2P4.Phi()
     << setw(8) << Za.l2P4.Energy()
     << setw(8) << Za.l2Charge
     << setw(8) << Za.l2Isolation
     << endl;

  if (Za.fsrPhoP4.Et() > 2) {
    os << " -- FSR Photon: " << endl;
    os << "       eT     eta     phi   energy" 
       << endl; 
    os << setw(9) << Za.fsrPhoP4.Et()
       << setw(8) << Za.fsrPhoP4.Eta()
       << setw(8) << Za.fsrPhoP4.Phi()
       << setw(9) << Za.fsrPhoP4.Energy()
       << endl;
  }
  os << " -- Z Properties:" << endl;
  os << "  flavour    mass  massDiff" 
     << endl; 
  os << setw(9) << ((!Za.flavour) ? "muon" : "electron") 
     << setw(8) << Za.mass
     << setw(10) << Za.massDiff
     << endl << endl;
}
bool PhysicsObjSelector::sameFlavourZPair(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
  if ( (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::ee) ||
       (Z1Cand.flavour == ZType::ee   && Z2Cand.flavour == ZType::mumu) ) return false;
  return true;  
}
void PhysicsObjSelector::syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate Z2, int nJets,
                                    float jet1Pt, float jet2Pt, std::ostream& os) {
  //{runX}:{lumiX}:{eventX}:{mass4l:.2fX}:{mZ1:.2fX}:{mZ2:.2fX}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
  //{njets30:dX}:  {jet1pt:.2fX}:{jet2pt:.2fX}:{category} 

  os << std::fixed << setprecision(2);
  os << run << ":"
     << lumi << ":"
     << event << ":"
     << (Z1.l1P4 + Z1.l2P4 + Z1.fsrPhoP4 + Z2.l1P4 + Z2.l2P4 +  Z2.fsrPhoP4).M() << ":"
     << Z1.mass << ":"
     << Z2.mass << ":"
     << 0 << ":"
     << nJets <<  ":"
     << jet1Pt << ":"
     << jet2Pt 
     << endl;
}

void PhysicsObjSelector::syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate Z2, int nJets,
                                    float jet1Pt, float jet2Pt, const int category,std::ostream& os) {
  //{runX}:{lumiX}:{eventX}:{mass4l:.2fX}:{mZ1:.2fX}:{mZ2:.2fX}::{D_bkg^kin:.3f}:{D_bkg:.3f}:{D_gg:.3f}:{D_HJJ^VBF:.3f}:{D_0-:.3f}:
  //{njets30:dX}:  {jet1pt:.2fX}:{jet2pt:.2fX}:{category} 

  os << std::fixed << setprecision(2);
  os << run << ":"
     << lumi << ":"
     << event << ":"
     << (Z1.l1P4 + Z1.l2P4 + Z1.fsrPhoP4 + Z2.l1P4 + Z2.l2P4 +  Z2.fsrPhoP4).M() << ":"
     << Z1.mass << ":"
     << Z2.mass << ":"
     << 0 << ":"
     << 0 << ":"
     << 0 << ":"
     << 0 << ":"
     << 0 << ":"
     << nJets <<  ":"
     << jet1Pt << ":"
     << jet2Pt << ":"
     << category
     << endl;
}

double PhysicsObjSelector::getEleRhoEffectiveArea(double eta) {
  double area;
  if (eta >= 0.0 && eta < 0.8) area = 0.1830;
  else if (eta >= 0.8 && eta < 1.3) area = 0.1734;
  else if (eta >= 1.3 && eta < 2.0) area = 0.1077;
  else if (eta >= 2.0 && eta < 2.2) area = 0.1565;
  else if (eta >= 2.2) area = 0.2680;
  return area;
}
// -- from twiki --
// for each photon attached to a lepton of each Z, test if the ll(g) candidate would satisfy the two criteria:
// 4 < m(llg))) < 100 GeV
// |m(llg)))-91.1876| < |m(ll)-91.1876| 
void PhysicsObjSelector::selectFSRPhoforLepton(const TLorentzVector& lep1P4, const TLorentzVector& lep2P4,
					       const vector<vhtm::PackedPFCandidate>& lepprePhoVec,
					       vector<vhtm::PackedPFCandidate>& lepPhoVec) 
{
  for (const auto& v: lepprePhoVec) {
    const TLorentzVector& pfPhoP4 = getP4(v);
    double mass = (lep1P4 + lep2P4 + pfPhoP4).M();
    if ((mass > 4. && mass < 100.)
	&& std::fabs(mass - MZnominal) < std::fabs((lep1P4 + lep2P4).M() - MZnominal))
      lepPhoVec.push_back(v);
  }
  if (lepPhoVec.size() > 1) 
    std::sort(lepPhoVec.begin(), lepPhoVec.end(), PtComparator<vhtm::PackedPFCandidate>());
}
// -- from twiki --
// if within the same Z multiple photons satisfy these criteria, pick the best one according to this logic:
//    - if there's at least one photon with pT > 4 GeV, pick the one with highest pT
//    - if all photons have pT < 4 GeV, pick the one that has the smallest dR to its closest lepton. 
int PhysicsObjSelector::selectBestFSRforLepton(const TLorentzVector& lepP4, const vector<vhtm::PackedPFCandidate>& lepPhoVec) {
  int index = -1;
  if (lepPhoVec[0].pt > 4.) index = 0;
  else {
    double dRmin = 999.;
    for (unsigned int i = 0; i < lepPhoVec.size(); ++i) {
      const TLorentzVector& phoP4 = getP4(lepPhoVec[i]);
      double dR = phoP4.DeltaR(lepP4);
      if (dR < dRmin) {
	index = i;
	dRmin = dR;
      }
    }
  }
  return index;
}
void PhysicsObjSelector::addFSRtoAltZ(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, 
				      int l1indx, int l2indx, TLorentzVector& altZP4, const std::string& tag, bool verbose) {
  if (Z1Cand.fsrWithLep == l1indx) {
    if (std::fabs((altZP4+Z1Cand.fsrPhoP4).M() - MZnominal) < std::fabs(altZP4.M() - MZnominal)) {
      if (verbose) cout << "--- SmartCut: adding Photon, l1indx: " << l1indx 
			<< ", Photon pT: " << Z1Cand.fsrPhoP4.Pt() << tag 
			<< endl;
      altZP4 += Z1Cand.fsrPhoP4;
    }
  }
  if (Z2Cand.fsrWithLep == l2indx) {
    if (std::fabs((altZP4+Z2Cand.fsrPhoP4).M() - MZnominal) < std::fabs(altZP4.M() - MZnominal)) {
      if (verbose) cout << "--- SmartCut: adding Photon, l2indx: " << l2indx 
			<< ", Photon pT: " << Z2Cand.fsrPhoP4.Pt() << tag 
			<< endl;
      altZP4 += Z2Cand.fsrPhoP4;
    }
  }
}
// -- from twiki --
// For each FSR photon that was selected, exclude that photon from the photon isolation sum of all four leptons of both Zs
// (Note that in Run I we were instead doing this only for the Z the photon is attached to) if it was in the isolation cone
// and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons;
// note: these requirements should probably be rechecked for consistency with isolation algorithms)
double PhysicsObjSelector::computeElectronIso(const vhtm::Electron& ele, const ZCandidate& Za, const TLorentzVector& lP4,
					      const std::vector<TLorentzVector>& fsrVec, 
					      double vetoCone, double isoCone, bool verbose) const
{
  double phoEtSum = 0.;
  for (const auto& v: fsrVec) {
    if ((std::fabs(ele.scEta) < 1.479 || lP4.DeltaR(v) > vetoCone) && lP4.DeltaR(v) < isoCone)
      phoEtSum += v.Et();
  }
  const vhtm::Event& evt = eventColl()->at(0);
  double fGridRhoFastjetAll = evt.fGridRhoFastjetAll;
  double iso = ele.chargedHadronIso + 
    std::max(0., ele.neutralHadronIso + ele.photonIso - phoEtSum -
	     getEleRhoEffectiveArea(std::fabs(lP4.Eta())) * fGridRhoFastjetAll);

  if (verbose)
    cout << "lepton isolation: " << iso << ", " << lP4.Pt() << ", "
	 << ele.chargedHadronIso << ", "
	 << ele.neutralHadronIso << ", "
	 << ele.photonIso << ", "
	 << phoEtSum << ", " 
	 << getEleRhoEffectiveArea(std::fabs(lP4.Eta())) << ", " 
	 << fGridRhoFastjetAll
	 << endl;

  return iso/lP4.Pt();
}
bool PhysicsObjSelector::ZeeIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone, double isoCone, bool verbose) const {
  const auto& elePhotonPairVec = getTightElePhotonPairList();
  
  const vhtm::Electron& ele1 = elePhotonPairVec.at(Za.l1Index).first;
  double l1iso = computeElectronIso(ele1, Za, Za.l1P4, fsrVec, vetoCone, isoCone, verbose);
  Za.l1Isolation = l1iso;

  const vhtm::Electron& ele2 = elePhotonPairVec.at(Za.l2Index).first;
  double l2iso = computeElectronIso(ele2, Za, Za.l2P4, fsrVec, vetoCone, isoCone, verbose);
  Za.l2Isolation = l2iso;
  
  if (Za.l1Isolation < 0.5 && Za.l2Isolation < 0.5) return true;
  return false;
}
// -- from twiki --
// For each FSR photon that was selected, exclude that photon from the photon isolation sum of all four leptons of both Zs
// (Note that in Run I we were instead doing this only for the Z the photon is attached to) if it was in the isolation cone
// and outside the isolation veto (ΔR>0.01 for muons and (ele->supercluster()->eta() < 1.479 || dR > 0.08) for electrons;
// note: these requirements should probably be rechecked for consistency with isolation algorithms)
double PhysicsObjSelector::computeMuonIso(const vhtm::Muon& mu, const ZCandidate& Za, const TLorentzVector& lP4,
					  const std::vector<TLorentzVector>& fsrVec, 
					  double vetoCone, double isoCone, bool verbose) const 
{
  double phoEtSum = 0.;
  for (const auto& v: fsrVec) {
    if (lP4.DeltaR(v) > vetoCone && lP4.DeltaR(v) < isoCone)
      phoEtSum += v.Et();
  }
  double iso = mu.sumChargedHadronPt
    + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt - phoEtSum - 0.5 * mu.sumPUPt);
  if (verbose) 
    cout << "l1 isolation: " << iso << ", " << lP4.Pt() << ", "
	 << mu.sumChargedHadronPt << ", "
	 << mu.sumNeutralHadronEt << ", "
	 << mu.sumPhotonEt << ", "
	 << phoEtSum << ", " 
	 << mu.sumPUPt
	 << endl;

  return iso/lP4.Pt();
}  
bool PhysicsObjSelector::ZmumuIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone, double isoCone, bool verbose) const {
  const auto& muPhotonPairVec = getTightMuPhotonPairList();

  const vhtm::Muon& mu1 = muPhotonPairVec.at(Za.l1Index).first;
  double l1iso = computeMuonIso(mu1, Za, Za.l1P4, fsrVec, vetoCone, isoCone, verbose);
  Za.l1Isolation = l1iso;
  
  const vhtm::Muon& mu2 = muPhotonPairVec.at(Za.l2Index).first;
  double l2iso = computeMuonIso(mu2, Za, Za.l2P4, fsrVec, vetoCone, isoCone, verbose);
  Za.l2Isolation = l2iso;
  
  if (Za.l1Isolation < 0.4 && Za.l2Isolation < 0.4) return true;
  return false;
}
// Function to compute the isolation for the ZZ leptons and check the Z masses
// If they pass, form a candidate ZZ pair
void PhysicsObjSelector::ZZisoMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec) {
  std::vector<TLorentzVector> fsrVec;
  fsrVec.push_back(Za.fsrPhoP4);
  fsrVec.push_back(Zb.fsrPhoP4);
  // -- from twiki --
  // Require all four leptons to pass the isolation requirement

  // First Z Candidate
  bool ZalepIso = false; 
  if (Za.flavour == ZType::mumu)
    ZalepIso = ZmumuIso(Za, fsrVec);
  else if (Za.flavour == ZType::ee)
    ZalepIso = ZeeIso(Za, fsrVec);

  // Second Z candidate 
  bool ZblepIso = false;
  if (Zb.flavour == ZType::mumu)
    ZblepIso = ZmumuIso(Zb, fsrVec);
  else if (Zb.flavour == ZType::ee)
    ZblepIso = ZeeIso(Zb, fsrVec);
  
  // -- from twiki --
  // Require both Z candidate masses (computed including the FSR photons, if present) to be 12 < m(ll(g)) < 120 GeV
  bool ZZmasscond = (Za.mass > 12 && Za.mass < 120) && (Zb.mass > 12 && Zb.mass < 120);

  // if all the leptons pass isolation and Z pair pass the mass cut, we have a ZZ Candidate
  if (ZalepIso && ZblepIso && ZZmasscond) {
    // -- from twiki --
    // define the Z1 as the one with mass closest to the nominal mZ; require mZ1 > 40 GeV. The other Z is the Z2.
    if (Za.massDiff < Zb.massDiff) {
      if (Za.mass > 40.)
        ZZVec.push_back({Za, Zb});
    } 
    else {
      if (Zb.mass > 40.)
	ZZVec.push_back({Zb, Za});
    }
  }
}
void PhysicsObjSelector::printP4(const TLorentzVector& lv, const string& tag) {
  cout << setprecision(3);
  cout << tag << " = (" 
       << setw(7) << lv.Pt() 
       << "," << setw(7) << lv.Eta() 
       << "," << setw(7) << lv.Phi() 
       << "," << setw(7) << lv.Energy() << ")" 
       << endl;
}
void PhysicsObjSelector::showEfficiency(const string& hname, const string slist[], const string& tag) {
  cout << ">>> " << tag << " Efficiency" << endl;
  TH1 *h = AnaUtil::getHist1D(hname);
  if (h) {
    cout << setw(64) << "CutFlow"
         << setw(10) << "Events"
         << setw(10) << "AbsEff"
         << setw(10) << "RelEff"
         << endl;
    cout << setprecision(3);
    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; ++i)
      cout << setw(64) << slist[i-1]
           << setprecision(1) 
           << setw(10) << int(h->GetBinContent(i))
           << setprecision(5) 
           << setw(10) << ((h->GetBinContent(1)>0) ? h->GetBinContent(i)/h->GetBinContent(1) : 0.0)
           << setw(10) << ( i == 1?1.0:(h->GetBinContent(i-1)>0) ? h->GetBinContent(i)/h->GetBinContent(i-1) : 0.0)
           << endl;
  }
}
void PhysicsObjSelector::dumpEvent(bool dumpGen, bool showEvent, ostream& os) const {
  os << std::setprecision(3);
  const vhtm::Event& evt = eventColl()->at(0);
  double fGridRhoFastjetAll = evt.fGridRhoFastjetAll;

  // Event
  if (showEvent) showEventNumber(os);

  // Muons
  if (muonColl()->size()) {
    os << " -- # Muons: " << muonColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  global tracker      PF  nMatch  Type       SIP ghostCleaned  reliso"
       << endl;
    int indx = 0;
    for (auto const& muon: *muonColl()) {
      double iso = muon.sumChargedHadronPt + std::max(0., muon.sumNeutralHadronEt + muon.sumPhotonEt - 0.5 * muon.sumPUPt);
      os << setw(6) << ++indx 
	 << setw(8) << muon.pt 
	 << setw(8) << muon.eta
	 << setw(8) << muon.phi
         << setw(8) << muon.charge
	 << setw(8) << muon.dxyPV
	 << setw(8) << muon.dzPV
	 << setw(8) << muon.isGlobalMuon
	 << setw(8) << muon.isTrackerMuon
	 << setw(8) << muon.isPFMuon
	 << setw(8) << muon.matches
	 << setw(6) << muon.muonBestTrackType
	 << setw(10) << muon.dB3D/muon.edB3D
	 << setw(13) << muon.isghostCleaned
         << setw(8) << iso/muon.pt
	 << endl;
    }
  }
  // Electrons
  if (electronColl()->size()) {
    os << " -- # Electrons: " << electronColl()->size() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  misHit       SIP   isBDT  reliso crossCleaned"
       << endl; 
    int indx = 0;
    for (auto const& electron: *electronColl()) {
      double iso = electron.chargedHadronIso + std::max(0., electron.neutralHadronIso
	         + electron.photonIso - getEleRhoEffectiveArea(std::fabs(electron.eta)) * fGridRhoFastjetAll);
      os << setw(6)  << ++indx
	 << setw(8)  << electron.pt 
	 << setw(8)  << electron.eta
	 << setw(8)  << electron.phi
         << setw(8)  << electron.charge
	 << setw(8)  << electron.dxyPV
	 << setw(8)  << electron.dzPV
	 << setw(8)  << electron.missingHits
	 << setw(10) << electron.dB3D/electron.edB3D
	 << setw(8)  << electronBDT(electron)
         << setw(8)  << iso/electron.pt
         << setw(13)  << crossCleaned(electron)
	 << endl;
    }
  }
  // Photons
  int npho = 0;
  for (const auto& pfcand: *packedPFCandidateColl()) {
    // Preselection: pT > 2 GeV, |eta| < 2.4 
    if (pfcand.pdgId != 22 || pfcand.pt < 2. || std::fabs(pfcand.eta) >= 2.4) continue;
    ++npho;
  }
  if (npho) {
    os << " -- Photons: " << endl;
    os << "  indx      pT     eta     phi  reliso  passedSCVeto closestLepton   dRmin"
       << endl; 
    int indx = 0;
    for (const auto& pfcand: *packedPFCandidateColl()) {
      // Preselection: pT > 2 GeV, |eta| < 2.4 
      if (pfcand.pdgId != 22 || pfcand.pt < 2. || std::fabs(pfcand.eta) >= 2.4) continue;
      int muindx = -1, elindx = -1;
      string ltype = "none";
      double dRmin = findClosestLepton(pfcand, muindx, elindx);
      if (elindx > -1) ltype = "electron";
      else if (muindx > -1) ltype = "muon";
      double iso = pfcand.isolationMap.at("c30").at(0) 
	+ pfcand.isolationMap.at("c30").at(2) 
	+ pfcand.isolationMap.at("c30").at(3) 
	+ pfcand.isolationMap.at("c30").at(4);
      os << setw(6) << ++indx
	 << setw(8) << pfcand.pt 
	 << setw(8) << pfcand.eta
	 << setw(8) << pfcand.phi
	 << setw(8) << iso/pfcand.pt
	 << setw(14) << (passedSuperClusterVeto(pfcand) ? "yes" : "no")
	 << setw(14) << ltype
	 << setw(8) << dRmin
	 << endl;
    }
  }
  // Selected vertices
  if (vertexColl()->size()) {
    os << " -- # Vertices: " << vertexColl()->size() << endl;
    os << "  indx      ndf     rho    chi2     dxy       z  isfake" 
       << endl; 
    int indx = 0;
    for (auto const& vtx: *vertexColl()) {
      double dxy = std::sqrt(pow(vtx.x, 2) + pow(vtx.y, 2));
      os << setw(6) << ++indx
	 << setw(9) << vtx.ndf
	 << setw(8) << vtx.rho
	 << setw(8) << vtx.chi2
	 << setw(8) << dxy
	 << setw(8) << vtx.z
	 << setw(8) << (vtx.isfake ? "yes" : "no")
	 << endl;
    }
  }
  //Jets
   if (jetColl()->size()) {
    os << " -- # Jets: " << jetColl()->size() << endl;
    os << "  indx      pT     eta     phi  jpumva    bDiscriminator"//    leptonCleaned"
       << endl;
    int indx = 0;
    for (auto const& jet: *jetColl()) {
      os << setw(6)  << ++indx
	 << setw(8)  << jet.pt 
	 << setw(8)  << jet.eta
         << setw(8)  << jet.phi
	 << setw(8)  << jet.jpumva
         << setw(8)  << jet.combinedInclusiveSecondaryVertexV2BJetTags
         //<< setw(8)  << ( jetLeptoncleaning( jet,AnaUtil::cutValue(jetCutMap(), "dRlep"))  ? "yes" : "no" )
         << endl;
    }
   }
  if (isMC() && dumpGen) dumpGenInfo(os); 
}
