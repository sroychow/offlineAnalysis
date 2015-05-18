//
//  FSRStudy.h
//
//
//  Created by Suvankar Roy Chowdhury on 13/03/15.
//
//

#ifndef __FSRStudy__
#define __FSRStudy__

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


class FSRStudy : public PhysicsObjSelector {
    
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
    int chargel1;
    int chargel2;
    int flavour;//0for mumu//1 for ee//2 for others
    double mass;
    double massDiff;
  };

  FSRStudy();
  virtual ~FSRStudy();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  void selectEvent(){};
  virtual void bookHistograms();


  //std::vector<vhtm::Tau> tauVec;
  //std::vector<vhtm::TriggerObject> trigObjList;
public:
  //std::map<std::string, double> _evselCutMap;
  //std::vector<std::string> _triggerPathTagList;
  bool _dumpEvent;
  const double MZnominal=91.1876;
  double Mzdiff;//|difference of Z1-Z mass|

  
  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Electron> looseeleVec,tighteleVec;
  std::vector<vhtm::Muon> loosemuVec,tightmuVec;
  std::vector<vhtm::PackedPFCandidate> pfPhoFSR;
  std::vector< std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> >> eleFSRPhopair;
  std::vector< std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > muFSRPhopair;

  
};
#endif

