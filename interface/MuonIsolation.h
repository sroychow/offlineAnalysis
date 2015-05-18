#ifndef __MuonIsolation__hh
#define __MuonIsolation__hh

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


class MuonIsolation : public AnaBase {
    
public:

  MuonIsolation();
  virtual ~MuonIsolation();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  
  void clearLists();
  void selectEvent(){};
  virtual void bookHistograms();
  void readRhofitSlope();
  void fillHistoforCone(std::string c,double mupt,double ch,double nu,double pu,int ngoodVtx,double fGridRhoFastjetAll);
  
public:
  //int nProbe[15];
  //int nSingleCut[15];

  std::vector<vhtm::Vertex> vtxList;
public:
  bool _dumpEvent;
  std::map<std::string,std::string> cone_;
};
#endif
