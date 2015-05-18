#ifndef __PhysicsObjSelector__hh
#define __PhysicsObjSelector__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "AnaBase.h"
#include "ZCandidate.h"

class PhysicsObjSelector: public AnaBase {
 public:
  PhysicsObjSelector();
  virtual ~PhysicsObjSelector() {}
  virtual void bookHistograms();
  virtual bool beginJob();
  
  const std::vector<vhtm::Muon>& getLooseMuList() const {return loosemuVec_;}
  const std::vector<vhtm::Muon>& getTightMuList() const {return tightmuVec_;}
  const std::vector<vhtm::Electron>& getLooseEleList() const {return looseeleVec_;}
  const std::vector<vhtm::Electron>& getTightEleList() const {return tighteleVec_;}
  const std::vector<vhtm::Jet>& getLeptonCleanedJetList() const {return jetVec_;}
  const std::vector<vhtm::PackedPFCandidate>& getFSRPhotonList() const {return photonVec_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getLooseElePhotonPairList() const {return looseElePhotonPairVec_;}
  const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& getTightElePhotonPairList() const {return tightElePhotonPairVec_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getLooseMuPhotonPairList() const {return looseMuPhotonPairVec_;}
  const std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& getTightMuPhotonPairList() const {return tightMuPhotonPairVec_;}
  void setTightElePhotonPairList(const std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > >& list) {tightElePhotonPairVec_ = list;}
  const int getNJets() const { return nJets_; }  
  const int getNbJets() const { return nbJets_; }  

  void findObjects(double evRho,double wt=1);
  void muonSelector(double wt=1);
  void electronSelector(double wt=1);
  void photonSelector();
  void jetSelector(double wt=1);
  static bool jetpuMVAid(const vhtm::Jet& jet);
  bool jetLeptoncleaning( const vhtm::Jet& jet,const float dR); 
  void setEventGridRho( double evRho ) { fGridRhoFastjetAll_ = evRho; }
  void clear();
  // must be defined inside the header to be effective
  template <class T>
  static TLorentzVector getP4(const T& obj) {
    TLorentzVector lv;
    lv.SetPtEtaPhiE(obj.pt, obj.eta, obj.phi, obj.energy);
    return lv;
  }
  bool passedSuperClusterVeto(const vhtm::PackedPFCandidate& pfcand, bool verbose=false) const;
  double findClosestLepton(const vhtm::PackedPFCandidate& photon, int& muindx, int& elindx) const;
  void leptonCrossCleaning();
  bool crossCleaned(const vhtm::Electron& electron) const;
  bool ZmumuIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.01, double isoCone=0.4, bool verbose=false) const;
  bool ZeeIso(ZCandidate& Za, const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.08, double isoCone=0.4, bool verbose=false) const;
  void ZZisoMass(ZCandidate& Za, ZCandidate& Zb, std::vector<std::pair<ZCandidate, ZCandidate> >& ZZVec);
  void dumpEvent(bool dumpGen=false, bool showEvent=false, ostream& os=std::cout) const;
  double computeMuonIso(const vhtm::Muon& mu, const ZCandidate& Za, const TLorentzVector& lP4,
			const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.01, double isoCone=0.4, bool verbose=false) const;
  double computeElectronIso(const vhtm::Electron& ele, const ZCandidate& Za, const TLorentzVector& lP4,
			    const std::vector<TLorentzVector>& fsrVec, double vetoCone=0.08, double isoCone=0.4, bool verbose=false) const;
  // static functions
  static void printP4(const TLorentzVector& lv, const std::string& tag);
  static double getEleRhoEffectiveArea(double eta);
  static bool sameFlavourZPair(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand);
  static bool electronBDT(const vhtm::Electron& electron);
  static void printZCandidate(const ZCandidate& Za, const std::string& tag, std::ostream& os=std::cout);
  static void syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate Z2, int nJets,
                         float jet1Pt, float jet2Pt, std::ostream& os=std::cout);
 
  static void syncDumper(int run, int lumi, int event, const ZCandidate& Z1, const ZCandidate Z2, int nJets,
                         float jet1Pt, float jet2Pt, const int category,std::ostream& os); 
  static void selectFSRPhoforLepton(const TLorentzVector& lep1P4, const TLorentzVector& lep2P4,
				    const std::vector<vhtm::PackedPFCandidate>& lep1prePhoVec,
				    std::vector<vhtm::PackedPFCandidate>& lep1PhoVec);
  static int selectBestFSRforLepton(const TLorentzVector& lepP4, const std::vector<vhtm::PackedPFCandidate>& lepPhoVec);
  static void addFSRtoAltZ(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, 
			   int l1indx, int l2indx, TLorentzVector& altZP4, const std::string& tag, bool verbose=false);
  static void showEfficiency(const std::string& hname, const std::string slist[], const std::string& tag);

 public:
  static const double MZnominal;
  enum ZType {
    mumu = 0, ee, wrong
  };
 
 private:
  bool dumpEvent_;
  std::vector<vhtm::Muon> loosemuVec_, tightmuVec_;
  std::vector<vhtm::Electron> looseeleVec_, tighteleVec_;
  std::vector<vhtm::PackedPFCandidate> photonVec_;
  std::vector<vhtm::Jet> jetVec_;
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > looseElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > tightElePhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > looseMuPhotonPairVec_; 
  std::vector<std::pair<vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > tightMuPhotonPairVec_; 
  bool searchedEle_, 
    searchedMu_, 
    searchedPhoton_;
  double fGridRhoFastjetAll_;
  int nJets_;
  int nbJets_;
};
#endif
