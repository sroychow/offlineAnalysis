#ifndef __LeptonIsolation__hh
#define __LeptonIsolation__hh

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
#include "PhysicsObjSelector.h"

class LeptonIsolation : public PhysicsObjSelector {
    
public:

  LeptonIsolation();
  virtual ~LeptonIsolation();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  
  void clearLists();
  void selectEvent(){};
  virtual void bookHistograms();
  void bookHistograms(TString folder);
  void readRhofitSlope();
  
  void fillHistoforCone(std::string c,double mupt,double ch_had,double ch_lep,double nu_had,double nu_photon,double pu,
                        int ngoodVtx,double fGridRhoFastjetAll,const std::map<std::string, double>& lepCutMap,TString& leptype);
  template <typename T>
    void getLeptonIsolationInfo( std::vector<T> lepvec, double rho,int ngoodVtx,const std::map<std::string, double>& lepCutMap,
                                 TString leptype );

  std::vector<vhtm::Vertex> vtxList;  
  bool _dumpEvent;
  std::map<std::string,std::string> cone_;
    
};
#endif
