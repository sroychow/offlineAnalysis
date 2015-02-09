#ifndef __ElectronIsolation__hh
#define __ElectronIsolation__hh

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


class ElectronIsolation : public AnaBase {
    
public:

  ElectronIsolation();
  virtual ~ElectronIsolation();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  //bool readJob(const std::string& jobFile, int& nFiles);
  //void printJob(std::ostream& os=std::cout) const;
  
  void clearLists();
 // void findMuonIDInfo();
 // void computeMuonEff();

  void selectEvent(){};
  virtual void bookHistograms();
  void readRhofitSlope();
public:
  //int nProbe[15];
  //int nSingleCut[15];

  std::vector<vhtm::Vertex> vtxList;
  //std::vector<vhtm::Electron> eleList;
  //std::vector<vhtm::Tau> tauList;
  //std::vector<vhtm::Jet> bjetList;
  //std::vector<vhtm::TriggerObject> trigObjList;
  //std::vector<vhtm::Muon> probeMuonList;
  //std::vector<vhtm::Muon> tagMuonList;
  std::vector<vhtm::Electron> eventElectrons;
public:
  //std::map<std::string, double> _evselCutMap;
  //std::vector<std::string> _triggerPathTagList;
  bool _dumpEvent;
  std::map<std::string,double> rhisoSlope;
};
#endif
