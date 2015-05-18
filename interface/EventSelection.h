#ifndef __EventSelection__hh
#define __EventSelection__hh

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


class EventSelection : public PhysicsObjSelector {
    
public:
  enum EventType
  {
   mmmm=0,eeee,eemm
  };
  enum Ztype
  {
   mumu=0,ee,wrong
  };
  struct Zcand {
    TLorentzVector l1P;
    TLorentzVector l2P;
    TLorentzVector fsrPhoP;
    int chargel1;
    int chargel2;
    int flavour;//0for mumu//1 for ee//2 for others
    double mass;
    double massDiff;
  };

  EventSelection();
  virtual ~EventSelection();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  void selectEvent(){};
  virtual void bookHistograms();
 // template <typename T>
 // bool Z1Selector();
 // template <typename T1,typename T2>
 // bool Z2Selector();
 // template <typename T>
 // void plotZmass(std::string ztype,std::pair<T,T> Z);
 // template <typename T1,typename T2>
 // void plotHmass(std::pair<T1,T1> Z1,std::pair<T2,T2> Z2);
 // void plotHmass(std::string hname,int val);//plot Higgs mass from struct method
 // template <typename T>
 // void ZSelector( std::vector< std::pair< T, std::vector<vhtm::PackedPFCandidate> > > lepFSRVec);
    
  
    void ZSelector( std::vector< std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > lepFSRVec);
    
    void ZSelector( std::vector< std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> > > lepFSRVec);
 // template <typename T>
 // void Z2Selector(const std::vector<T>& lepVec);  
  static bool ZSorter(const Zcand &a, const Zcand &b) {
    return a.massDiff < b.massDiff;
  }
    template <typename T>
    double calcDBLepIso(T,double phoEnergy);
    template <typename T>
    double calcRhoLepIso(T lep,double phoEnergy, double effArea, double rho);
    double getEleRhoEffectiveArea( double fSCeta);
    int ZZselector(Zcand, Zcand);
    void final4lSelector(int run, int lumi,int event);
//Functions for gen Level
  template <typename T>
  TLorentzVector getP4(const T& obj);
  int getGenDauPgd(const vhtm::GenParticle gp);
  double getHmassfromZdau(const vhtm::GenParticle Z1,const vhtm::GenParticle Z2);
  bool genOk();
  //void createDiscriminant(TLorentzVector,TLorentzVector,TLorentzVector,TLorentzVector);
  void syncDumper(int run, int lumi,int event,Zcand Z1,Zcand Z2);
  
    void selectFSRPhoforlepton( TLorentzVector lep1P4,TLorentzVector lep2P4,
                                std::vector<vhtm::PackedPFCandidate> lep1prePhoVec,
                                std::vector<vhtm::PackedPFCandidate>& lep1PhoVec);
    int selectBestFSRforlepton( TLorentzVector lepP4, std::vector<vhtm::PackedPFCandidate> lepPhoVec);
public:

    std::vector<vhtm::Vertex> vtxList;
    std::vector<vhtm::Electron> looseeleVec,tighteleVec;
    std::vector<vhtm::Muon> loosemuVec,tightmuVec;
    std::vector<vhtm::PackedPFCandidate> pfPhoFSR;
    std::vector< std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> >> eleFSRPhopair;
    std::vector< std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > muFSRPhopair;
    std::vector< std::pair< Zcand, Zcand> > ZZvec;
  bool _dumpEvent;
  const double MZnominal=91.188;
  double Mzdiff;//|difference of Z1-Z mass|
  std::vector< std::pair<vhtm::Muon,vhtm::Muon> > Z1mu,Z2mu;
  std::vector< std::pair<vhtm::Electron,vhtm::Electron> > Z1ele,Z2ele;
  bool Z1present,Z2present,Z1tomu,Z1toele;
  std::vector<Zcand> myZ1,myZ2;
  std::vector<Zcand> myZ;
  Zcand Z2actual;
  int n4mu,n2e2mu,n4e,nZ1clash,ntauev=0;
  std::vector<vhtm::GenParticle> genZ;
  int evtype;
   double fGridRhoFastjetAll_;
  ofstream syncDumpf;
  
};
#endif

