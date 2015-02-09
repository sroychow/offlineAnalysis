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


class EventSelection : public AnaBase {
    
public:
  enum eventType
  {
   mmmm,eeee,eemm
  };
  enum Ztype
  {
   mumu,ee
  };
  struct Zcand {
    TLorentzVector l1P;
    TLorentzVector l2P;
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
  void muonSelector(); 
  void tauSelector(); 
  void electronSelector();  
 // template <typename T>
 // bool Z1Selector(const std::vector<T>& lepVec,std::vector< std::pair<T,T> >& Z1);
 // template <typename T1,typename T2>
 // bool Z2Selector(std::vector< std::pair<T1,T1> >& Z1,std::vector< std::pair<T2,T2> >& Z2,const std::vector<T2>& lepVec);
 // template <typename T>
 // void plotZmass(std::string ztype,std::pair<T,T> Z);
 // template <typename T1,typename T2>
 // void plotHmass(std::pair<T1,T1> Z1,std::pair<T2,T2> Z2);
  void plotHmass();//plot Higgs mass from struct method
  template <typename T>
  void ZSelector(const std::vector<T>& lepVec);
 // template <typename T>
 // void Z2Selector(const std::vector<T>& lepVec);  
  static bool ZSorter(const Zcand &a, const Zcand &b) {
    return a.massDiff < b.massDiff;
  }
//Functions for gen Level
  template <typename T>
  TLorentzVector getP4(const T& obj);
  int getGenDauPgd(const vhtm::GenParticle gp);
  double getHmassfromZdau(const vhtm::GenParticle Z1,const vhtm::GenParticle Z2);
  bool genOk();
public:
  //int nProbe[15];
  //int nSingleCut[15];

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Electron> eleVec;
  std::vector<vhtm::Muon> muVec;
  std::vector<vhtm::Tau> tauVec;
  //std::vector<vhtm::TriggerObject> trigObjList;
public:
  //std::map<std::string, double> _evselCutMap;
  //std::vector<std::string> _triggerPathTagList;
  bool _dumpEvent;
  const double MZnominal=91.1876;
  double Mzdiff;//|difference of Z1-Z mass|
  std::vector< std::pair<vhtm::Muon,vhtm::Muon> > Z1mu,Z2mu;
  std::vector< std::pair<vhtm::Electron,vhtm::Electron> > Z1ele,Z2ele;
  bool Z1present,Z2present,Z1tomu,Z1toele;
  std::vector<Zcand> myZ1,myZ2;
  Zcand Z2actual;
  int n4mu,n2e2mu,n4e,nZ1clash,ntauev=0;
  std::vector<vhtm::GenParticle> genZ;
};
#endif

