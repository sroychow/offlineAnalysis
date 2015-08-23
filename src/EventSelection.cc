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
#include "TH3.h"
#include "TProfile.h"

#include "AnaUtil.h"
#include "HZZ4lUtil.h"
#include "EventSelection.h"
#include "GenAnalysis.h"

//#include "MEMCalculators.h"
//using namespace MEMNames;

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

using namespace vhtm;

// -----------
// Constructor
// -----------
EventSelection::EventSelection()
  : PhysicsObjSelector(),
  checkGen_(true),
  dumpGenInfo_(false),
  useEventList_(false),
  selectEvType_(false),
  evtype_(-1),
  doKDcalc_(false),
  dumpFilename_("syncDumpFile.txt")
{
}
// ----------
// Destructor
// ----------
EventSelection::~EventSelection() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool EventSelection::beginJob() 
{ 
  // Open the output ROOT file (in AnaBase)
  PhysicsObjSelector::beginJob();

  histf()->cd();
  histf()->mkdir("EventSelection");
  
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void EventSelection::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("EventSelection");

  //------- Object PLots -----------------------------------------------
  new TH1D("nGoodmuon", "Number of Good muons(with selection cuts) per event", 20, -0.5, 19.5);
  new TH1D("nGoodelectron", "Number of Good electrons(with selection cuts) per event", 20, -0.5, 19.5);
  new TH1D("nGoodtau", "Number of tau(with selection cuts) per event", 20, -0.5, 19.5);

  // Z and Hplots
  new TH1D("nZcand", "Number of selected Zcandidates per event", 20, 0, 20);
  new TH1F("massnZcand", "Mass of selected Zcandidates", 100, 0., 200.);
  new TH1F("massZ1", "Mass of selected Z1", 200, 0., 200.);
  new TH1D("higgsMass","mass of the selected Higgs candidate",300,0,300);
	
  new TH1D("evtCutFlow", "Event CutFlow", 10, -0.5, 9.5);

  new TH1F("nEvntflavReco", "Flavour of ZZ decay at reco", 4, -0.5, 3.5);

  new TH1D("dRlepZal1Zal2", "dRlepZal1Zal2", 100, 0, 5);
  new TH1D("dRlepZbl1Zbl2", "dRlepZbl1Zbl2", 100, 0, 5);
  new TH1D("dRlepZal1Zbl1", "dRlepZal1Zbl1", 100, 0, 5);
  new TH1D("dRlepZal1Zbl2", "dRlepZal1Zbl2", 100, 0, 5);
  new TH1D("dRlepZal2Zbl1", "dRlepZal2Zbl1", 100, 0, 5);
  new TH1D("dRlepZal2Zbl2", "dRlepZal2Zbl2", 100, 0, 5);

  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);

  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void EventSelection::clearLists() {
  PhysicsObjSelector::clear();

  vtxList_.clear();
  genZList_.clear();
  ZCandList_.clear();
  leptauCandList_.clear();
  ZZPairVec_.clear();
  evtype_ = -1;
}
// -------------------
// The main event loop
// -------------------
void EventSelection::eventLoop()
{
  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = std::max(10000, nEvents()/1000);
  
  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  int fevt = (firstEvent() > -1) ? firstEvent() : 0;
  int levt = (lastEvent() > -1) ? lastEvent() : nEvents();
  cout << ">>> Event range: [" << fevt << ", " << levt -1 << "]" << endl;
  int eeetau =0, mumumutau=0, eemutau =0;
  for (int ev = fevt; ev < levt; ++ev) {
    clearEvent(); // reset tree variables
    clearLists();
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));
    const Event& evt = eventColl()->at(0);
    
    histf()->cd();
    
    // For data or for MC without pileup
    puevWt_ = 1;
    if (isMC() && usePUWt()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
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
    if (ev%nPrint == 0 || firstEvent() > -1)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << run
	   << " Lumis# " << lumis
	   << " Event# " << setw(8) << event << " >>> "
	   << " Events proc. " << setw(8) << ((firstEvent() > -1) ? ev - firstEvent() : ev)
	   << endl;
    
    if (useEventList_ && eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
      //if (eventIdMap().find(mkey.str()) != eventIdMap().end()) continue;
    }
    
    histf()->cd();
    histf()->cd("EventSelection");
    AnaUtil::fillHist1D("evtCutFlow", 0, puevWt_);
    
    bool bookGenHist = (ev > 0) ? false:true;
    GenAnalysis gAna(*genParticleColl(),histf(),bookGenHist);
    //if(  gAna.genOk() )
    //  cout << "GenPass = " << gAna.genOk() << std::endl;
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    double ngoodVtx = vtxList_.size();
    
    histf()->cd();
    histf()->cd("EventSelection");
    
    AnaUtil::fillHist1D("evtCutFlow", 1, puevWt_);
    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt_); 
    
    // is event triggered?
    if (useTrigger() && !isTriggered(true, false)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2, puevWt_);
    
    // at least 1 good PV
    if (ngoodVtx < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, puevWt_);
    
    // main analysis object selection
    double vz = (vtxList_.size() > 0) ? vtxList_.at(0).z : -999;
    
    findObjects(vz,puevWt_);
    
    // access selected objects 
    const auto& elePhotonPairVec = getTightElePhotonPairList();
    const auto& muPhotonPairVec  = getTightMuPhotonPairList();
    const auto& tauVec = getTauList();
    
    histf()->cd();
    histf()->cd("EventSelection");
    AnaUtil::fillHist1D("nGoodmuon", muPhotonPairVec.size(), puevWt_);
    AnaUtil::fillHist1D("nGoodelectron", elePhotonPairVec.size(), puevWt_);
    AnaUtil::fillHist1D("nGoodtau", tauVec.size(), puevWt_);
    
    //if ( elePhotonPairVec.size() >= 3 && tauVec.size() >= 1 ) eeetau++;
    if( elePhotonPairVec.size() < 2 && muPhotonPairVec.size() < 2  )      continue;
    AnaUtil::fillHist1D("evtCutFlow", 4, puevWt_);
    if( tauVec.empty() )                                                  continue;
    AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
    
    if(muPhotonPairVec.size() >= 2)    ZSelector<vhtm::Muon>(muPhotonPairVec,ZCandList_);
    if(elePhotonPairVec.size() >= 2)   ZSelector<vhtm::Electron>(elePhotonPairVec,ZCandList_);
    
    if( ZCandList_.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6, puevWt_);
  
    std::sort(ZCandList_.begin(), ZCandList_.end(),dmComparator);
    ZCandidate selectedZ = ZCandList_.at(0);

    if( muPhotonPairVec.size() > 2 )     leptauSelector<vhtm::Muon>(muPhotonPairVec,tauVec,selectedZ,leptauCandList_);
    if( elePhotonPairVec.size() > 2 )    leptauSelector<vhtm::Electron>(elePhotonPairVec,tauVec,selectedZ,leptauCandList_);
    
    if( leptauCandList_.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7, puevWt_);
    
    std::sort( leptauCandList_.begin(), leptauCandList_.end(),massComparator);
    
    AnaUtil::fillHist1D("massZ1",selectedZ.mass, puevWt_);
    AnaUtil::fillHist1D("higgsMass", leptauCandList_[0].mass, puevWt_);
    
  }


// Analysis is over
std::cout <<"Neeetau=" << eeetau << "\t Nmumumutau=" << mumumutau << "\t Neemutau=" << eemutau << std::endl;
endJob();
}

//find a opp sign lepton tau pair
template <typename T>
void EventSelection::leptauSelector(const std::vector<std::pair<T, std::vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec, 
                                    const std::vector<vhtm::Tau> tauVec,const ZCandidate& eventZ,std::vector<ZCandidate>& leptauList) {
  for( auto& lep: lepPhotonPairVec) {
    const TLorentzVector& lepP4 = HZZ4lUtil::getP4(lep.first);
    if( lepP4 == eventZ.l1P4 || lepP4 == eventZ.l2P4 )   continue;
    for( auto& tau: tauVec ) {
      if (lep.first.charge + tau.charge != 0) continue;
      const TLorentzVector& tauP4 = HZZ4lUtil::getP4(tau);
      if( lepP4.DeltaR(tauP4) < 0.02 )        continue;
      ZCandidate ltau;
      ltau.flavour = -1;

      ltau.l1Index = -1;
      ltau.l1P4 = lepP4;
      ltau.l1Charge = lep.first.charge;

      ltau.l2Index = -1;
      ltau.l2P4 = tauP4;
      ltau.l2Charge = tau.charge;
      ltau.fsrWithLep = -1;
      ltau.mass = ( lepP4 + tauP4 ).M();
      ltau.massDiff = -99999999;
      leptauList.push_back(ltau);
    }
  } 
}


// create unique lepton pair combination giving a Z statisfying Z2 mass cuts
// and push them into a vector
template <typename T>
void EventSelection::ZSelector(const vector<pair<T, vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec, std::vector<ZCandidate>& candList) {
  for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
    const auto& ip = lepPhotonPairVec[i];
    const TLorentzVector& lep1P4 = HZZ4lUtil::getP4(ip.first);
    for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
      const auto& jp = lepPhotonPairVec[j];

      // opposite charge
      if (ip.first.charge + jp.first.charge != 0) continue; 
      const TLorentzVector& lep2P4 = HZZ4lUtil::getP4(jp.first);

      // Select FSR photon(s) for lepton[i]
      // both the leptons are needed to form Z mass
      vector<vhtm::PackedPFCandidate> lep1PhoList;
      if (!ip.second.empty())
        HZZ4lUtil::selectFSRPhoforLepton(lep1P4, lep2P4, ip.second, lep1PhoList);

      // Select FSR photon(s) for lepton[j]
      vector<vhtm::PackedPFCandidate> lep2PhoList;
      if (!jp.second.empty())
        HZZ4lUtil::selectFSRPhoforLepton(lep2P4, lep1P4, jp.second, lep2PhoList);

      int l1phoidx = (!lep1PhoList.empty()) ? HZZ4lUtil::selectBestFSRforLepton(lep1P4, lep1PhoList) : -1;
      int l2phoidx = (!lep2PhoList.empty()) ? HZZ4lUtil::selectBestFSRforLepton(lep2P4, lep2PhoList) : -1;

      // Now choose the correct FSR Photon for this ll pair
      TLorentzVector phoP4;    // LorentzVector of the finally selected FSR
      int withLep = -1;
      if (l1phoidx == -1 && l2phoidx == -1) {
        phoP4.SetPtEtaPhiE(0., 0., 0., 0);  
      }
      else if (l1phoidx > -1 && l2phoidx == -1) {
        phoP4 = HZZ4lUtil::getP4(lep1PhoList[l1phoidx]);
        withLep = 1;
      }
      else if (l1phoidx == -1 && l2phoidx > -1) {
        phoP4 = HZZ4lUtil::getP4(lep2PhoList[l2phoidx]);
        withLep = 2;
      }
      else {
        const vhtm::PackedPFCandidate& lep1Photon = lep1PhoList[l1phoidx];
        const vhtm::PackedPFCandidate& lep2Photon = lep2PhoList[l2phoidx];
        if (lep1Photon.pt > 4. || lep2Photon.pt > 4.) {
          if (lep1Photon.pt > lep2Photon.pt) {
            phoP4 = HZZ4lUtil::getP4(lep1Photon);
            withLep = 1;
          }
          else {
            phoP4 = HZZ4lUtil::getP4(lep2Photon);
            withLep = 2;
          }
        }
        else {
          const TLorentzVector& l1PhoP4 = HZZ4lUtil::getP4(lep1Photon);
          const TLorentzVector& l2PhoP4 = HZZ4lUtil::getP4(lep2Photon);
          bool dRcond = l1PhoP4.DeltaR(lep1P4) < l2PhoP4.DeltaR(lep2P4);
          if (dRcond) {
            phoP4 = l1PhoP4;
            withLep = 1;
          }
          else {
            phoP4 = l2PhoP4;
            withLep = 2;
          }
        }
      }
      ZCandidate ztmp;
      if (typeid(jp.first) == typeid(vhtm::Muon))
        ztmp.flavour = HZZ4lUtil::ZType::mumu;
      else if (typeid(jp.first) == typeid(vhtm::Electron))
        ztmp.flavour = HZZ4lUtil::ZType::ee;
      else 
        ztmp.flavour = -1;

      ztmp.l1Index = i;
      ztmp.l1P4 = lep1P4;
      ztmp.l1Charge = ip.first.charge;

      ztmp.l2Index = j;
      ztmp.l2P4 = lep2P4;
      ztmp.l2Charge = jp.first.charge;

      ztmp.fsrWithLep = withLep; 
      ztmp.fsrPhoP4 = phoP4;
      double Zmass = (lep1P4 + lep2P4 + phoP4).M();
      ztmp.mass = Zmass;
      ztmp.massDiff = std::fabs(Zmass - HZZ4lUtil::MZnominal);

      candList.push_back(ztmp);
    }
  }
}
// For ZZ pairs
void EventSelection::ZZselector() {
  histf()->cd();
  histf()->cd("EventSelection");

  // Note. ZZPairVec_ will be overwritten!
  std::vector<std::pair<ZCandidate, ZCandidate> > list = ZZPairVec_;
  ZZPairVec_.clear();

  int zzcounter[] = {0,0,0,0,0,0,0};
  for (const auto& v: list) {
    ++zzcounter[0];
    const ZCandidate& Z1Cand = v.first;
    const ZCandidate& Z2Cand = v.second;

    // -- from twiki --
    // dR(eta,phi)>0.02 between each of the four leptons (to remove ghosts)
    double dra1a2 = Z1Cand.l1P4.DeltaR(Z1Cand.l2P4);
    double drb1b2 = Z2Cand.l1P4.DeltaR(Z2Cand.l2P4);  
    double dra1b1 = Z1Cand.l1P4.DeltaR(Z2Cand.l1P4);
    double dra1b2 = Z1Cand.l1P4.DeltaR(Z2Cand.l2P4);
    double dra2b1 = Z1Cand.l2P4.DeltaR(Z2Cand.l1P4);
    double dra2b2 = Z1Cand.l2P4.DeltaR(Z2Cand.l2P4);
    AnaUtil::fillHist1D("dRlepZal1Zal2", dra1a2, puevWt_);
    AnaUtil::fillHist1D("dRlepZbl1Zbl2", drb1b2, puevWt_);
    AnaUtil::fillHist1D("dRlepZal1Zbl1", dra1b1, puevWt_);
    AnaUtil::fillHist1D("dRlepZal1Zbl2", dra1b2, puevWt_);
    AnaUtil::fillHist1D("dRlepZal2Zbl1", dra2b1, puevWt_);
    AnaUtil::fillHist1D("dRlepZal2Zbl2", dra2b2, puevWt_);
    bool dRlep = dra1a2 > 0.02 && 
      drb1b2 > 0.02 && 
      dra1b1 > 0.02 && 
      dra1b2 > 0.02 && 
      dra2b1 > 0.02 && 
      dra2b2 > 0.02;
    if (!dRlep) continue;
    ++zzcounter[1];

    vector<TLorentzVector> lepP4List;
    lepP4List.push_back(Z1Cand.l1P4);
    lepP4List.push_back(Z1Cand.l2P4);
    lepP4List.push_back(Z2Cand.l1P4);
    lepP4List.push_back(Z2Cand.l2P4);
    std::sort(lepP4List.begin(), lepP4List.end(), PtComparatorTL<TLorentzVector>());

    // -- from twiki --
    // any two leptons of the four have pt > 20/10
    if (lepP4List[0].Pt() <= 20) continue;
    ++zzcounter[2];

    if (lepP4List[1].Pt() <= 10) continue;
    ++zzcounter[3];

    // -- from twiki --
    // QCD suppression: m(ll)>4 GeV on all four Opposite Signed pairs that can be made with the four leptons (regardless of flavour).
    // FSR photons are not used in computing m(ll).
    TLorentzVector ZaP4, ZbP4;
    if (Z1Cand.l1Charge + Z2Cand.l1Charge == 0) {
      ZaP4 = Z1Cand.l1P4 + Z2Cand.l1P4;
      ZbP4 = Z1Cand.l2P4 + Z2Cand.l2P4;
    }
    else {
      ZaP4 = Z1Cand.l1P4 + Z2Cand.l2P4;
      ZbP4 = Z1Cand.l2P4 + Z2Cand.l1P4;
    }
    if (Z1Cand.mass <= 4 || Z2Cand.mass <= 4 || ZaP4.M() <= 4. || ZbP4.M() <= 4.) continue;
    ++zzcounter[4];

    // -- from twiki --
    // additional "smart cut" to reject 4mu/4e pairs where the alternative pairing looks like an on-shell Z+low-mass ll. 
    // (This is required to avoid a background increase when selecting the "best candidate" after all cuts; cf. Simon's slides. 
    // The requirement is: !( |mZa-mZ| < |mZ1-mZ| && mZb<12), where Za and Zb are the mass-sorted alternative pairing Z candidates 
    // (Za being the one closest to the nominal Z mass). FSR photons that are present in the candidate under study are considered as 
    // FSR for Za and Zb leptons with the same FSR-attaching requirements (see above). This means that essentially the same FSR 
    // recovery criteria used for ZZ are applied to Za and Zb, but starting from the list of photons already attached to the ZZ 
    // candidate, and without recomputing isolation. 
    if (HZZ4lUtil::sameFlavourZPair(Z1Cand, Z2Cand)) {
      TLorentzVector altZ1P4, altZ2P4;
      if (Z1Cand.l1Charge + Z2Cand.l1Charge == 0) {
        altZ1P4 = Z1Cand.l1P4 + Z2Cand.l1P4;
        HZZ4lUtil::addFSRtoAltZ(Z1Cand, Z2Cand, 1, 1, altZ1P4, "(lep1,lep1)"); 

        altZ2P4 = Z1Cand.l2P4 + Z2Cand.l2P4;
        HZZ4lUtil::addFSRtoAltZ(Z1Cand, Z2Cand, 2, 2, altZ2P4, "(lep2,lep2)"); 
      }
      else {
        altZ1P4 = Z1Cand.l1P4 + Z2Cand.l2P4;
        HZZ4lUtil::addFSRtoAltZ(Z1Cand, Z2Cand, 1, 2, altZ1P4, "(lep1,lep2)"); 

        altZ2P4 = Z1Cand.l2P4 + Z2Cand.l1P4;
        HZZ4lUtil::addFSRtoAltZ(Z1Cand, Z2Cand, 2, 1, altZ2P4, "(lep2,lep1)"); 
      }
      if (std::fabs(altZ2P4.M() - HZZ4lUtil::MZnominal) < std::fabs(altZ1P4.M() - HZZ4lUtil::MZnominal)) {
        TLorentzVector lv = altZ1P4;
        altZ1P4 = altZ2P4;
        altZ2P4 = lv;
      }
      if (0) 
        cout << "--- SmartCut Leptons: Z1: (" << Z1Cand.l1Index << "," << Z1Cand.l2Index
          << "), Z2: (" << Z2Cand.l1Index << "," << Z2Cand.l2Index << ")" << endl
          << "--- SmartCut: (" << Z1Cand.flavour << ", " << Z2Cand.flavour
          << ", " << Z1Cand.mass << ", " << altZ1P4.M() 
          << ", " << Z2Cand.mass << ", " << altZ2P4.M() << ")"
          << endl;
      if (std::fabs(altZ1P4.M() - HZZ4lUtil::MZnominal) < Z1Cand.massDiff && altZ2P4.M() < 12) {
        if (0) cout << " -- SmartCut: skip ZZ Candidate" << endl;
        continue;
      }
    }
    ++zzcounter[5];

    // --- from twiki --
    // m(4l) > 70 GeV
    double mass4l = (Z1Cand.l1P4 + Z1Cand.l2P4 + Z1Cand.fsrPhoP4 + Z2Cand.l1P4 + Z2Cand.l2P4 + Z2Cand.fsrPhoP4).M();
    if (mass4l <= 70.) continue;
    ++zzcounter[6];

    // Z1,Z2 satisfy the conditions, hence fill them in ZZPairVec_
    ZZPairVec_.push_back({Z1Cand, Z2Cand});
  }
  for (unsigned int i = 0; i < NEL(zzcounter); ++i)
    if (zzcounter[i] > 0) AnaUtil::fillHist1D("zzSelCutFlow", i, puevWt_);
}
void EventSelection::finalZZSelector(int run, int lumi, int event) {
  histf()->cd();
  histf()->cd("EventSelection");

  // -- from twiki --
  // If more than one ZZ candidate survives, choose the one with the Z1 closest in mass to nominal Z. If two or more
  // candidates include the same Z1 and differ just by the Z2, choose the one with the highest-pT Z2 leptons.
  // (Other Z pairing choices,e.g. by best Dbkg, could be considered later).
  ZCandidate Z1Cand = ZZPairVec_[0].first;
  ZCandidate Z2Cand = ZZPairVec_[0].second;
  for (unsigned int i = 1; i < ZZPairVec_.size(); ++i)  {
    // choose the one with the Z1 closest in mass to nominal Z
    double massDiff = ZZPairVec_[i].first.massDiff;
    if (massDiff < Z1Cand.massDiff) {
      Z1Cand = ZZPairVec_[i].first;
      Z2Cand = ZZPairVec_[i].second;
    } 
    else if (massDiff == Z1Cand.massDiff) {
      // If two or more candidates include the same Z1 and differ just by the Z2, choose
      // the one with the highest-pT Z2 leptons.
      if ((ZZPairVec_[i].second.l1P4.Pt() + ZZPairVec_[i].second.l2P4.Pt())/2. > (Z2Cand.l1P4.Pt() + Z2Cand.l2P4.Pt())/2.)
        Z2Cand = ZZPairVec_[i].second;
    }
  }
  double mass4l = (Z1Cand.l1P4 + Z1Cand.l2P4 + Z1Cand.fsrPhoP4 + Z2Cand.l1P4 + Z2Cand.l2P4 + Z2Cand.fsrPhoP4).M();
  AnaUtil::fillHist1D("massZ1", Z1Cand.mass, puevWt_);
  AnaUtil::fillHist1D("massZ2", Z2Cand.mass, puevWt_);
  AnaUtil::fillHist1D("mass4l", mass4l, puevWt_);

  cout << "---- New Event" << endl;
  showEventNumber(true);

  const auto& jetVec = getLeptonCleanedLooseJetList();
  int nJets = jetVec.size();

  double jet1Pt = -1.0, jet2Pt= -1.0;
  TLorentzVector jet1P4, jet2P4;
  jet1P4.SetPtEtaPhiE(0.,0.,0.,0.);
  jet2P4.SetPtEtaPhiE(0.,0.,0.,0.);
  if (jetVec.size()) {
    jet1Pt = jetVec[0].pt;
    jet1P4 = HZZ4lUtil::getP4(jetVec[0]);
  }
  if (jetVec.size() > 1) {
    jet2Pt = jetVec[1].pt;
    jet2P4 = HZZ4lUtil::getP4(jetVec[1]);
  }

  int nextlep = findExtraLeptons(Z1Cand, Z2Cand);
  AnaUtil::fillHist1D("nextraLepton", nextlep, puevWt_);

  int cat = findEventCategory(4+nextlep, jetVec, getNbJets(), Z1Cand, Z2Cand); 
  AnaUtil::fillHist1D("evntcatReco", cat, puevWt_);

  // synchronization exercise dump

  HZZ4lUtil::syncDumper(run, lumi, event, Z1Cand, Z2Cand, nJets, jet1Pt, jet2Pt, cat, syncDumpf_);

  // Event Flavour
  if (Z1Cand.flavour == HZZ4lUtil::ZType::mumu && Z2Cand.flavour == HZZ4lUtil::ZType::mumu)
    AnaUtil::fillHist1D("nEvntflavReco", 0, puevWt_);
  else if (Z1Cand.flavour == HZZ4lUtil::ZType::ee && Z2Cand.flavour == HZZ4lUtil::ZType::ee)
    AnaUtil::fillHist1D("nEvntflavReco", 1, puevWt_);
  else if ( (Z1Cand.flavour == HZZ4lUtil::ZType::mumu && Z2Cand.flavour == HZZ4lUtil::ZType::ee) || 
      (Z1Cand.flavour == HZZ4lUtil::ZType::ee   && Z2Cand.flavour == HZZ4lUtil::ZType::mumu) )
    AnaUtil::fillHist1D("nEvntflavReco", 2, puevWt_);
  else
    AnaUtil::fillHist1D("nEvntflavReco", 3, puevWt_);

  // dump event
  cout << " --- #Z Cadidates: " << ZCandList_.size() 
    << ", #ZZ Candidates(final): " << ZZPairVec_.size() 
                                      << ", flavour: "  << (HZZ4lUtil::sameFlavourZPair(Z1Cand, Z2Cand) ? "same" : "different") 
                                        << endl;
  cout << "--- mass4l: " << mass4l << " GeV" 
    << endl;
  cout << "--- fGridRhoFastjetAll: " << getEventGridRho() << endl;
  HZZ4lUtil::printZCandidate(Z1Cand, " >> Z Candidate 1");
  HZZ4lUtil::printZCandidate(Z2Cand, " >> Z Candidate 2");
  dumpEvent();
}
// Function to calculate number of extra tight leptons passing isolation apart from leptons in
// ZZ candidate
int EventSelection::findExtraLeptons(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
  int nextlep = 0;
  for (const auto& mu: getTightMuList()) {
    TLorentzVector muP4 = HZZ4lUtil::getP4(mu);
    if (muP4 == Z1Cand.l1P4 || muP4 == Z1Cand.l2P4 || muP4 == Z2Cand.l1P4 || muP4 == Z2Cand.l2P4) continue;
    if (HZZ4lUtil::pfiso(mu)/mu.pt < 0.4) nextlep++;
  }
  for (const auto& ele: getTightEleList()) {
    TLorentzVector eleP4 = HZZ4lUtil::getP4(ele);
    if (eleP4 == Z1Cand.l1P4 || eleP4 == Z1Cand.l2P4 || eleP4 == Z2Cand.l1P4 || eleP4 == Z2Cand.l2P4) continue;
    if (HZZ4lUtil::pfiso(ele, getEventGridRho())/ele.pt < 0.5) nextlep++;
  }
  return nextlep;
}
//
//  exactly 4 leptons + at least 2 jets + at most 1 b-tag jet in the event +Djet>0.5
//   VBF-tagged category (#2)
//
//  exactly 4 leptons + at least one pair of jets passing { |eta|<2.4, pT>40, 60<mjj<120 } + pT,4l>m4l,
//  OR: exactly 4 leptons + exactly 2 jets, both with b-tag
//    VH-hadronic-tagged category (#4)
//
//  no more than 2 jets + no b-tagged jet in the event; + at least 5 leptons
//    VH-leptonic-tagged category (#3)
//
//  at least 3 jets, of which at least 1 b-tagged,
//  OR: at least 5 leptons
//    ttH-tagged category (#5)
//
//  at least 1 jet
//    1-jet-tagged category (#1)
//
//  remaining events
//    untagged category (#0) 
//
//  Djet = 0.18*|eta_j1 - eta_j2| + 1.92e-04 * mj1j2
//
int EventSelection::findEventCategory(int nleptons, const std::vector<vhtm::Jet>& jetList, int nbjets,
    const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, bool verbose) {
  int njets = jetList.size();

  TLorentzVector j1P4, j2P4;
  j1P4.SetPtEtaPhiE(0.,0.,0.,0.);
  j2P4.SetPtEtaPhiE(0.,0.,0.,0.);
  if (njets) j1P4 = HZZ4lUtil::getP4(jetList[0]);
  if (njets > 1) j2P4 = HZZ4lUtil::getP4(jetList[1]);

  double djet = -1., mjj = 0.;
  TLorentzVector final4lP4 = Z1Cand.l1P4 + Z1Cand.l2P4 + Z1Cand.fsrPhoP4 + Z2Cand.l1P4 + Z2Cand.l2P4 + Z2Cand.fsrPhoP4;
  if (njets >= 2) {
    mjj = (j1P4 + j2P4).M();
    djet = 0.18 * std::fabs(j1P4.Eta() - j2P4.Eta()) + 1.92E-04 * mjj;
  }
  histf()->cd();
  histf()->cd("EventSelection");
  AnaUtil::fillHist1D("djet", djet, 1);

  int cat = 0;
  if (nleptons == 4 && njets >= 2 && nbjets <= 1 && djet > 0.5) 
    cat = 2;
  else if ( (nleptons == 4 && 
        njets >= 2 && EventSelection::hasJetPair(jetList) &&
        final4lP4.Pt() > final4lP4.M() ) ||
      (nleptons == 4 && njets == 2 && nbjets == 2) )
    cat = 4;
  else if (njets <= 2 && nbjets == 0 && nleptons >= 5)
    cat = 3;
  else if ( (njets >= 3 && nbjets >= 1) || nleptons >= 5 )
    cat = 5;
  else if (njets >= 1)
    cat = 1;

  if (verbose) {
    cout << "---- Event Category" << endl;
    cout << "  nlep  njet nbjet   jet1Pt  jet1Eta   jet2Pt  jet2Eta      mjj     4lPt      4lM     djet category"
      << endl;
    cout << setw(6) << nleptons
      << setw(6) << njets
      << setw(6) << nbjets
      << setw(9) << j1P4.Pt()
      << setw(9) << (j1P4.Pt() > 0 ? j1P4.Eta() : 99)
      << setw(9) << j2P4.Pt()
      << setw(9) << (j2P4.Pt() > 0 ? j2P4.Eta() : 99)
      << setw(9) << mjj
      << setw(9) << final4lP4.Pt() 
      << setw(9) << final4lP4.M()
      << setw(9) << djet
      << setw(9) << cat
      << endl << endl;
  }
  return cat;
}
bool EventSelection::hasJetPair(const std::vector<vhtm::Jet>& jetList) {
  for (unsigned int i = 0; i < jetList.size(); ++i) {
    auto const& j1 = jetList[i];
    auto const& j1P4 = HZZ4lUtil::getP4(j1);
    for (unsigned int j = i+1; j < jetList.size(); ++j) {
      auto const& j2 = jetList[j];
      auto const& j2P4 = HZZ4lUtil::getP4(j2);
      double mjj = (j1P4+j2P4).M();      
      if (0) cout << "mjj[" << i << ", " << j << "] = " << mjj << endl;
      if ( (std::fabs(j1.eta) < 2.4 && j1.pt > 40.) && 
          (std::fabs(j2.eta) < 2.4 && j2.pt > 40.) && 
          (mjj > 60. && mjj < 120.) ) return true;
    }
  }
  return false;
}

  void EventSelection::endJob() {
    syncDumpf_.close();
    closeFiles();

    histf()->cd();
    histf()->cd("EventSelection");
    string evlabels[] = {
      "Events processed",
      "Gen Filter",
      "Events passing trigger",
      "Events with > 0 good vertex",
      "Correct Channel",
      "4 leptons",
      "# of Z Candidates > 0",
      "# of Z Candidates > 1",
      "Candidate (mass + isolation)",
      "Candidate final cuts"
    };
    HZZ4lUtil::showEfficiency("evtCutFlow", evlabels, "Event Selection");  
    string zzlabels[] = {
      "All",
      "dR(Lep1,Lep2) Overlap",
      "Lep1 pT >= 20 GeV",
      "Lep2 pT >= 10 GeV",
      "Lepton Cross Mass",
      "Smart Cut",
      "Z Mass"
    };
    HZZ4lUtil::showEfficiency("zzSelCutFlow", zzlabels, "ZZ Selection");  

    histf()->cd();
    histf()->Write();
    histf()->Close();
    delete histf();
  }
  // -------------------------------------------------------------------------------
  // Poor man's way of a datacard. Each line between the 'START' and 'END' tags
  // is read in turn, split into words, where the first element is the 'key' and
  // the rest the value(s). If more than one values are present they are 
  // stored in a vector. No safety mechanism is in place. Any line with an unknown 
  // key is skipped. Comments lines should start with either '#' or '//', preferably
  // in the first column. Empty lines are skipped. The file containing the datacards 
  // is passed as the only argument of the program, there is no default
  // -------------------------------------------------------------------------------
  bool EventSelection::readJob(const string& jobFile, int& nFiles)
  {
    if (!AnaBase::readJob(jobFile, nFiles)) return false;

    static const int BUF_SIZE = 256;

    // Open the file containing the datacards
    ifstream fin(jobFile.c_str(), std::ios::in);    
    if (!fin) {
      cerr << "Input File: " << jobFile << " could not be opened!" << endl;
      return false;
    }

    char buf[BUF_SIZE];
    vector<string> tokens;
    while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
      string line(buf);
      if (line.empty() || line == "START") continue;   

      // enable '#' and '//' style comments
      if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
      if (line == "END") break;

      // Split the line into words
      AnaUtil::tokenize(line, tokens);
      assert(tokens.size() > 1);
      string key = tokens[0];
      string value = tokens[1];
      if (key == "useEventList")
        useEventList_ = (std::stoi(value.c_str()) > 0) ? true : false;
      else if (key == "dumpGenInfo")
        dumpGenInfo_ = (std::stoi(value.c_str()) > 0) ? true : false;
      else if (key == "checkGenInfo")
        checkGen_ = (std::stoi(value.c_str()) > 0) ? true : false;
      else if (key == "doKD")
        doKDcalc_ = (std::stoi(value.c_str()) > 0) ? true : false;
      else if (key == "syncDumpFile")
        dumpFilename_ = value.c_str();

      tokens.clear();
    }
    // Close the file
    fin.close();

    syncDumpf_.open(dumpFilename_, std::ios::out);
    if (!syncDumpf_) {
      cerr << "Output File: " << dumpFilename_ << " could not be opened!" << endl;
      return false;
    }

    selectEvType_ = (static_cast<int>(AnaUtil::cutValue(evselCutMap(), "selectEvType")) > 0) ? true : false;
    printJob();

    return true;
  }
  void EventSelection::printJob(ostream& os) const
  {
    AnaBase::printJob(os);
  }
