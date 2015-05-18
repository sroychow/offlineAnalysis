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
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
EventSelection::EventSelection()
  : PhysicsObjSelector(),
    dumpGenInfo_(false),
    useEventList_(false),
    selectEvType_(false),
    evtype_(-1)
{
  syncDumpf_.open("syncDumpFile.txt");
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
  // Open the output ROOT file (in AnaBase)
  PhysicsObjSelector::beginJob();

  histf()->cd();
  histf()->mkdir("EventSelection");
  histf()->mkdir("GeneratorLevel");
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
  new TH1D("nGoodmuon", "Number of Good muons(with selection cuts) per event", 20, 0, 20);
  new TH1D("nGoodelectron", "Number of Good electrons(with selection cuts) per event", 20, 0, 20);
  new TH1D("nextraLepton", "Number of extra tight leptons passing isolation in selected events(apart from 4l)", 20, 0, 20);
  // Z and Hplots
  new TH1D("nZcand", "Number of selected Zcandidates per event", 20, 0, 20);
  new TH1D("nZZcand", "Number of selected ZZ candidates per event", 20, 0, 20);
  new TH1F("massnZcand", "Mass of selected Zcandidates", 200, 0., 200.);
  new TH1F("massZ1", "Mass of selected Z1", 200, 0., 200.);
  new TH1F("massZ2", "Mass of selected Z2", 200, 0., 200.);
  new TH1F("mass4l", "4lepton mass", 200, 0., 200.);
  
  new TH1D("evtCutFlow", "Event CutFlow", 10, -0.5, 9.5);
  new TH1D("zzSelCutFlow", "ZZ Selection CutFlow", 7, -0.5, 6.5);
  
  new TH1D("dRlepZal1Zal2", "dRlepZal1Zal2", 100, 0, 5);
  new TH1D("dRlepZbl1Zbl2", "dRlepZbl1Zbl2", 100, 0, 5);
  new TH1D("dRlepZal1Zbl1", "dRlepZal1Zbl1", 100, 0, 5);
  new TH1D("dRlepZal1Zbl2", "dRlepZal1Zbl2", 100, 0, 5);
  new TH1D("dRlepZal2Zbl1", "dRlepZal2Zbl1", 100, 0, 5);
  new TH1D("dRlepZal2Zbl2", "dRlepZal2Zbl2", 100, 0, 5);
  
  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);
  
  // GEN Level histograms
  histf()->cd();
  histf()->cd("GeneratorLevel");
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
  for (int ev = 0; ev < nEvents(); ++ev) {
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
    if (ev%nPrint == 0)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << run
	   << " Lumis# " << lumis
	   << " Event# " << setw(8) << event << " >>> "
	   << " Events proc. " << setw(8) << ev
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

    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    double ngoodVtx = vtxList_.size();
    fGridRhoFastjetAll_ = evt.fGridRhoFastjetAll;
    if (AnaUtil::cutValue(evselCutMap(), "Isbkg")) { 
    }
    else {
      bool genPass = genOk();
      if (!genPass) {
	if (dumpGenInfo_) {
	  showEventNumber();
	  dumpGenInfo();
        }
        continue;
      }

      // crucial
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
      findObjects(fGridRhoFastjetAll_,puevWt_);

      
      // access selected objects 
      const auto& elePhotonPairVec = getTightElePhotonPairList();
      const auto& muPhotonPairVec  = getTightMuPhotonPairList();


      
      histf()->cd();
      histf()->cd("EventSelection");
      
      AnaUtil::fillHist1D("nGoodmuon", muPhotonPairVec.size(), puevWt_);
      AnaUtil::fillHist1D("nGoodelectron", elePhotonPairVec.size(), puevWt_);
      
      // one can select a particular chanel (4mu,4e,2e2mu)
      if (selectEvType_ && evtype_ != static_cast<int>(AnaUtil::cutValue(evselCutMap(), "evType"))) continue;
      AnaUtil::fillHist1D("evtCutFlow", 4, puevWt_);
      
      if (evtype_ == EventType::eemm) {
	if (elePhotonPairVec.size() >= 2 && muPhotonPairVec.size() >= 2) {
	  AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
	  ZSelector<vhtm::Electron>(elePhotonPairVec, ZCandList_);
	  ZSelector<vhtm::Muon>(muPhotonPairVec, ZCandList_);
	}
      }
      else if (evtype_ == EventType::eeee) {
	if (elePhotonPairVec.size() >= 4) {
	  AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
	  ZSelector<vhtm::Electron>(elePhotonPairVec, ZCandList_);
	}
      }
      else if (evtype_ == EventType::mmmm) {
	if (muPhotonPairVec.size() >= 4) {
	  AnaUtil::fillHist1D("evtCutFlow", 5, puevWt_);
	  ZSelector<vhtm::Muon>(muPhotonPairVec, ZCandList_);
	}
      }
      else {
	cerr << "eventLoop: Wrong evtype: " << evtype_ << endl;
        continue;
      }
      AnaUtil::fillHist1D("nZcand", ZCandList_.size(), puevWt_);
      // At least one Z candidate found 
      if (ZCandList_.size() >= 1) {
	AnaUtil::fillHist1D("evtCutFlow", 6, puevWt_);
	for (unsigned int i = 0; i < ZCandList_.size(); ++i)
	  AnaUtil::fillHist1D("massnZcand", ZCandList_[i].mass, puevWt_);
	
	// At least two Z candidates found 
	if (ZCandList_.size() >= 2) {
	  AnaUtil::fillHist1D("evtCutFlow", 7, puevWt_);
	  for (unsigned int i = 0; i < ZCandList_.size(); ++i) {
	    for (unsigned int j = i+1; j < ZCandList_.size(); ++j)
	      ZZisoMass(ZCandList_[i], ZCandList_[j], ZZPairVec_);
	  }
	  AnaUtil::fillHist1D("nZZcand", ZZPairVec_.size(), puevWt_);
	  if (ZZPairVec_.size() > 0) {
	    AnaUtil::fillHist1D("evtCutFlow", 8, puevWt_);
	    ZZselector(); // Note. ZZPairVec_ gets overwritten
	    if (ZZPairVec_.size() > 0) {
	      AnaUtil::fillHist1D("evtCutFlow", 9, puevWt_);
	      finalZZSelector(run, lumis, event);
	    }
	  }
	}
      }
    }
  }
  // Analysis is over
  endJob();
}
// create unique lepton pair combination giving a Z statisfying Z2 mass cuts
// and push them into a vector
template <typename T>
void EventSelection::ZSelector(const vector<pair<T, vector<vhtm::PackedPFCandidate> > >& lepPhotonPairVec, std::vector<ZCandidate>& candList) {
  for (unsigned int i = 0; i < lepPhotonPairVec.size(); ++i) {
    const auto& ip = lepPhotonPairVec[i];
    const TLorentzVector& lep1P4 = PhysicsObjSelector::getP4(ip.first);
    for (unsigned int j = i+1; j < lepPhotonPairVec.size(); ++j) {
      const auto& jp = lepPhotonPairVec[j];

      // opposite charge
      if (ip.first.charge + jp.first.charge != 0) continue; 
      const TLorentzVector& lep2P4 = PhysicsObjSelector::getP4(jp.first);
      
      // Select FSR photon(s) for lepton[i]
      // both the leptons are needed to form Z mass
      vector<vhtm::PackedPFCandidate> lep1PhoList;
      if (!ip.second.empty())
	PhysicsObjSelector::selectFSRPhoforLepton(lep1P4, lep2P4, ip.second, lep1PhoList);
      
      // Select FSR photon(s) for lepton[j]
      vector<vhtm::PackedPFCandidate> lep2PhoList;
      if (!jp.second.empty())
	PhysicsObjSelector::selectFSRPhoforLepton(lep2P4, lep1P4, jp.second, lep2PhoList);
      
      int l1phoidx = (!lep1PhoList.empty()) ? PhysicsObjSelector::selectBestFSRforLepton(lep1P4, lep1PhoList) : -1;
      int l2phoidx = (!lep2PhoList.empty()) ? PhysicsObjSelector::selectBestFSRforLepton(lep2P4, lep2PhoList) : -1;

      // Now choose the correct FSR Photon for this ll pair
      TLorentzVector phoP4;    // LorentzVector of the finally selected FSR
      int withLep = -1;
      if (l1phoidx == -1 && l2phoidx == -1) {
        phoP4.SetPtEtaPhiE(0., 0., 0., 0);  
      }
      else if (l1phoidx > -1 && l2phoidx == -1) {
	phoP4 = PhysicsObjSelector::getP4(lep1PhoList[l1phoidx]);
        withLep = 1;
      }
      else if (l1phoidx == -1 && l2phoidx > -1) {
	phoP4 = PhysicsObjSelector::getP4(lep2PhoList[l2phoidx]);
        withLep = 2;
      }
      else {
        const vhtm::PackedPFCandidate& lep1Photon = lep1PhoList[l1phoidx];
	const vhtm::PackedPFCandidate& lep2Photon = lep2PhoList[l2phoidx];
        if (lep1Photon.pt > 4. || lep2Photon.pt > 4.) {
          if (lep1Photon.pt > lep2Photon.pt) {
            phoP4 = PhysicsObjSelector::getP4(lep1Photon);
            withLep = 1;
          }
          else {
            phoP4 = PhysicsObjSelector::getP4(lep2Photon);
            withLep = 2;
          }
	}
	else {
	  const TLorentzVector& l1PhoP4 = PhysicsObjSelector::getP4(lep1Photon);
	  const TLorentzVector& l2PhoP4 = PhysicsObjSelector::getP4(lep2Photon);
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
        ztmp.flavour = PhysicsObjSelector::ZType::mumu;
      else if (typeid(jp.first) == typeid(vhtm::Electron))
        ztmp.flavour = PhysicsObjSelector::ZType::ee;
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
      ztmp.massDiff = std::fabs(Zmass - PhysicsObjSelector::MZnominal);

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
    if (PhysicsObjSelector::sameFlavourZPair(Z1Cand, Z2Cand)) {
      TLorentzVector altZ1P4, altZ2P4;
      if (Z1Cand.l1Charge + Z2Cand.l1Charge == 0) {
	altZ1P4 = Z1Cand.l1P4 + Z2Cand.l1P4;
	PhysicsObjSelector::addFSRtoAltZ(Z1Cand, Z2Cand, 1, 1, altZ1P4, "(lep1,lep1)"); 

	altZ2P4 = Z1Cand.l2P4 + Z2Cand.l2P4;
        PhysicsObjSelector::addFSRtoAltZ(Z1Cand, Z2Cand, 2, 2, altZ2P4, "(lep2,lep2)"); 
      }
      else {
	altZ1P4 = Z1Cand.l1P4 + Z2Cand.l2P4;
	PhysicsObjSelector::addFSRtoAltZ(Z1Cand, Z2Cand, 1, 2, altZ1P4, "(lep1,lep2)"); 

	altZ2P4 = Z1Cand.l2P4 + Z2Cand.l1P4;
	PhysicsObjSelector::addFSRtoAltZ(Z1Cand, Z2Cand, 2, 1, altZ2P4, "(lep2,lep1)"); 
      }
      if (std::fabs(altZ2P4.M() - PhysicsObjSelector::MZnominal) < std::fabs(altZ1P4.M() - PhysicsObjSelector::MZnominal)) {
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
      if (std::fabs(altZ1P4.M() - PhysicsObjSelector::MZnominal) < Z1Cand.massDiff && altZ2P4.M() < 12) {
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
  const auto& jetVec = getLeptonCleanedJetList();
  
  double jet1Pt,jet2Pt;
  if( !jetVec.empty() )  {
    jet1Pt = jetVec.at(0).pt;
    if( jetVec.size() > 1 ) jet2Pt = jetVec.at(1).pt;
    else                    jet2Pt = -1.;
  } else {
    jet1Pt = -1.; 
    jet2Pt = -1.;
  }        
  int nextlep = findNExtralepton(Z1Cand,Z2Cand);
  AnaUtil::fillHist1D("nextraLepton", nextlep, puevWt_);
  PhysicsObjSelector::syncDumper(run, lumi, event, Z1Cand, Z2Cand, jetVec.size(),jet1Pt,jet2Pt,syncDumpf_);
  dumpEvent();
#if 1
  // dump event
  cout << "---- New Event" << endl;
  showEventNumber(false);
  cout << " --- #Z Cadidates: " << ZCandList_.size() 
       << ", #ZZ Candidates(final): " << ZZPairVec_.size() 
       << ", flavour: "  << (PhysicsObjSelector::sameFlavourZPair(Z1Cand, Z2Cand) ? "same" : "different") 
       << endl;
  cout << "--- mass4l: " << mass4l << " GeV" 
       << endl;
  cout << "--- fGridRhoFastjetAll: " << fGridRhoFastjetAll_ << endl;
  PhysicsObjSelector::printZCandidate(Z1Cand, " >> Z Candidate 1");
  PhysicsObjSelector::printZCandidate(Z2Cand, " >> Z Candidate 2");
  dumpEvent();
#endif
}

//Function to calculate number of extra tight leptons passing isolation apart from leptons in
//ZZ candidate
int EventSelection::findNExtralepton(const ZCandidate& Z1,const ZCandidate& Z2) {
  int nextlep=0;
  for( auto& mu: getTightMuList() ) {
    TLorentzVector muVec = getP4(mu);
    if( muVec == Z1.l1P4 || muVec == Z1.l2P4 || muVec == Z2.l1P4 || muVec == Z2.l2P4 ) 
      continue;
    if( ( mu.sumChargedHadronPt + std::max(0., mu.sumNeutralHadronEt + mu.sumPhotonEt - 0.5 * mu.sumPUPt))/mu.pt < 0.4 )
      nextlep++;
  }
  for( auto& ele: getTightEleList() ) {
    TLorentzVector eleVec = getP4(ele);
    if( eleVec == Z1.l1P4 || eleVec == Z1.l2P4 || eleVec == Z2.l1P4 || eleVec == Z2.l2P4 ) 
      continue;
    if( ( ele.chargedHadronIso + 
          std::max(0., ele.neutralHadronIso + ele.photonIso  -
          getEleRhoEffectiveArea(std::fabs(ele.eta)) * fGridRhoFastjetAll_))/ele.pt < 0.5 )
      nextlep++;
  }
  return nextlep; 
}
//---------------------
// Gen Level Functions
//---------------------
int EventSelection::getGenDauPgd(const GenParticle& gp) {
  vector<int> d = gp.daughtIndices;
  vector<int> dpdgid;
  for (auto di: d) {
    if (di >= ngenparticle()) continue;
    const GenParticle& dgp = genParticleColl()->at(di);
    int pid = std::abs(dgp.pdgId);
    if (pid == 23) continue;
    dpdgid.push_back(pid);
  }
  if (dpdgid.empty()) return -1;
  else if (dpdgid[0] == 15) return dpdgid[1];
  else return dpdgid[0];
}
double EventSelection::getHmassfromZdau(const GenParticle& Z1, const GenParticle& Z2) {
  TLorentzVector lZ1, lZ2;
  lZ1.SetPtEtaPhiE(0., 0., 0., 0.);
  lZ2.SetPtEtaPhiE(0., 0., 0., 0.);

  const vector<int>& d1 = Z1.daughtIndices;
  for (auto di: d1) {
    if (di >= ngenparticle()) continue;
    const GenParticle& dgp = genParticleColl()->at(di);
    if (std::abs(dgp.pdgId) == 23) continue;
    lZ1 += getP4(dgp);
  }

  const vector<int>& d2 = Z2.daughtIndices;
  for (auto di: d2) {
    if (di >= ngenparticle()) continue;
    const GenParticle& dgp = genParticleColl()->at(di);
    if (std::abs(dgp.pdgId) == 23) continue;
    lZ2 += getP4(dgp);
  }
  return (lZ1+lZ2).M();
}
bool EventSelection::genOk() {
  histf()->cd();
  histf()->cd("GeneratorLevel");
  bool genPass = false;
  for (const auto& gp: *genParticleColl()) {
    int pdgid = std::abs(gp.pdgId);
    int status = gp.status;
    
    if (pdgid == 25) {
      TLorentzVector l = getP4(gp);
      AnaUtil::fillHist1D("genHmass", l.M(), puevWt_);
    }
    // if (pdgid == 23 && status==2 ) {
    else if (pdgid == 23 && status == 3 ) { // Z
      int mmid = -1;
      int index = getMotherId(gp, mmid);
      AnaUtil::fillHist1D("genZmother", mmid, puevWt_);
      if (index < 0 || std::abs(mmid) != 25) continue; // H -> ZZ
      genZList_.push_back(gp);
      TLorentzVector l1 = getP4(gp);
      AnaUtil::fillHist1D("genZmass", l1.M(), puevWt_);
    }
  }
  if (0) cout << "#Gen Z = " << genZList_.size() << endl;
  if (genZList_.size() < 2) return genPass;
  //  std::sort(genZList_.begin(), genZList_.end(), MassComparator<vhtm::GenParticle>());
  AnaUtil::fillHist1D("genZ1mass", getP4(genZList_[0]).M(), puevWt_); 
  AnaUtil::fillHist1D("genZ2mass", getP4(genZList_[1]).M(), puevWt_);
  AnaUtil::fillHist1D("genHmassfromZ", (getP4(genZList_[0])+getP4(genZList_[1])).M(), puevWt_);
  
  int nZ = genZList_.size();
  int* nMu  = new int[nZ];
  int* nEle = new int[nZ];
  for (int i = 0; i < nZ; ++i) {nMu[i] = 0; nEle[i] = 0;}
  int iz = 0;
  for (auto const& gp: genZList_) {
    // check e mu decay of Z
    vector<int> d = gp.daughtIndices;
    for (auto di: d) {
      if (di >= ngenparticle()) continue;
      const GenParticle& dgp = genParticleColl()->at(di);
      int pid = std::abs(dgp.pdgId);
      if (pid == 11) {
	++nEle[iz];
	AnaUtil::fillHist1D("nZdauflav", 3, puevWt_);
      }
      else if (pid == 13) {
	++nMu[iz];
	AnaUtil::fillHist1D("nZdauflav", 1, puevWt_); 
      }
      else if (pid == 15) {
	AnaUtil::fillHist1D("nZdauflav", 5, puevWt_);
        GenParticle xgp = (dgp.status == 3) ? genParticleColl()->at(dgp.daughtIndices[0]) : dgp;
	const GenParticle& ygp = genParticleColl()->at(xgp.daughtIndices[0]);
	if (std::abs(ygp.pdgId) == 15 && ygp.status == 2) xgp = ygp;
	const GenParticle& zgp = genParticleColl()->at(xgp.daughtIndices[0]);
	if (std::abs(zgp.pdgId) == 15 && zgp.status == 2) xgp = zgp;

	vector<int> dtauIndices = xgp.daughtIndices;
	for (auto in: dtauIndices) {
	  if (in >= ngenparticle()) continue;
	  const GenParticle& dtaugp = genParticleColl()->at(in);
	  int dtaupid = std::abs(dtaugp.pdgId);
	  if (dtaupid == 13) {
	    ++nMu[iz];
	    AnaUtil::fillHist1D("nZdauflav", 1, puevWt_); 
	  }
	  else if (dtaupid == 11 || dtaupid == 213) {
	    ++nEle[iz];
	    AnaUtil::fillHist1D("nZdauflav", 3, puevWt_);
	  }
	}
      }
      else AnaUtil::fillHist1D("nZdauflav", 7, puevWt_);
    }                
    iz++;
  }
  if (0) {
    for (int i = 0; i < nZ; ++i) {
      cout << "Z[" << i << "] daughters: " << nMu[i] << ", " << nEle[i] << endl;
    }
  }  
  // Event Flavour
  if (nMu[0] == 2 && nMu[1] == 2) {
    AnaUtil::fillHist1D("nEvntflav", 1, puevWt_);
    AnaUtil::fillHist1D("genHmassfromZdau", getHmassfromZdau(genZList_[0], genZList_[1]), puevWt_);
    genPass = true;
    evtype_ = EventType::mmmm;
  } 
  else if (nEle[0] == 2 && nEle[1] == 2) {
    AnaUtil::fillHist1D("nEvntflav", 3, puevWt_);
    AnaUtil::fillHist1D("genHmassfromZdau", getHmassfromZdau(genZList_[0], genZList_[1]), puevWt_);
    genPass = true;
    evtype_ = EventType::eeee;
  }
  else if ((nMu[0] == 2 && nEle[1] == 2) || (nEle[0] == 2 && nMu[1] == 2)) {
    AnaUtil::fillHist1D("nEvntflav", 5, puevWt_);
    AnaUtil::fillHist1D("genHmassfromZdau", getHmassfromZdau(genZList_[0], genZList_[1]), puevWt_);
    genPass = true;
    evtype_ = EventType::eemm;
  } 
  else {
    AnaUtil::fillHist1D("nEvntflav", 7, puevWt_);
  }
  delete [] nMu;
  delete [] nEle;

  if (genPass) {   
    if (0) cout << "Event Type = " << evtype_ << endl;

    // Z1 daughters
    if (getGenDauPgd(genZList_[0]) == 13 ) AnaUtil::fillHist1D("nZ1dauflav", 1, puevWt_);
    else if (getGenDauPgd(genZList_[0]) == 11) AnaUtil::fillHist1D("nZ1dauflav", 3, puevWt_);
    else if (getGenDauPgd(genZList_[0]) == 15) AnaUtil::fillHist1D("nZ1dauflav", 5, puevWt_);
    else AnaUtil::fillHist1D("nZ1dauflav", 7, puevWt_);
 
    // Z2 daughter
    if (getGenDauPgd(genZList_[1]) == 13 ) AnaUtil::fillHist1D("nZ2dauflav", 1, puevWt_);
    else if (getGenDauPgd(genZList_[1]) == 11) AnaUtil::fillHist1D("nZ2dauflav", 3, puevWt_);
    else if (getGenDauPgd(genZList_[1]) == 15) AnaUtil::fillHist1D("nZ2dauflav", 5, puevWt_);
    else AnaUtil::fillHist1D("nZ2dauflav", 7, puevWt_);
  }
  return genPass;
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
  PhysicsObjSelector::showEfficiency("evtCutFlow", evlabels, "Event Selection");  
  string zzlabels[] = {
    "All",
    "dR(Lep1,Lep2) Overlap",
    "Lep1 pT >= 20 GeV",
    "Lep2 pT >= 10 GeV",
    "Lepton Cross Mass",
    "Smart Cut",
    "Z Mass"
  };
  PhysicsObjSelector::showEfficiency("zzSelCutFlow", zzlabels, "ZZ Selection");  
  
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
  ifstream fin(jobFile.c_str(), ios::in);    
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
      useEventList_ = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = (atoi(value.c_str()) > 0) ? true : false;

    tokens.clear();
  }
  // Close the file
  fin.close();

  selectEvType_ = (static_cast<int>(AnaUtil::cutValue(evselCutMap(), "selectEvType")) > 0) ? true : false;
  printJob();

  return true;
}
void EventSelection::printJob(ostream& os) const
{
  AnaBase::printJob(os);
}
