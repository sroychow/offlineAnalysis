#ifndef __LeptonSync__hh
#define __LeptonSync__hh

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
#include "PhysicsObjSelector.h"

class LeptonSync : public PhysicsObjSelector {
    
public:
  enum EventType
  {
   mmmm=0,eeee,eemm
  };
  LeptonSync();
  virtual ~LeptonSync();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  
  void clearLists();

  void selectEvent(){};
  virtual void bookHistograms();
    
  void muonSelector();
  void electronSelector();
  void pfphotonSelector();

//Functions for gen Level
  template <typename T>
  TLorentzVector getP4(const T& obj);
  int getGenDauPgd(const vhtm::GenParticle gp);
  double getHmassfromZdau(const vhtm::GenParticle Z1,const vhtm::GenParticle Z2);
  bool genOk();
  void dumpMuon(std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> >);
  void dumpElectron(std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> >);
  double calcRhoLepIso(vhtm::Electron lep, double rho);
  double getEleRhoEffectiveArea( double fSCeta);
public:

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::GenParticle> genZ;
  std::vector<vhtm::Muon> loosemuVec_,tightmuVec_;
  std::vector<vhtm::Electron> looseeleVec_,tighteleVec_;
  std::vector<vhtm::PackedPFCandidate> pfphotonVec_;
  std::vector< std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> >> looseElePhoFSRpair_;
  std::vector< std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > looseMuPhoFSRpair_;
  
public:
  bool _dumpEvent;
  int evtype,run_,lumi_,event_;
  double rho_;
  ofstream lepsyncDumper_;
};
#endif
