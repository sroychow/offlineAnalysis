#ifndef __GenAnalysis__hh
#define __GenAnalysis__hh

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


class GenAnalysis : public AnaBase {
    
public:
  GenAnalysis();
  virtual ~GenAnalysis();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void clearLists();

  void selectEvent(){};
  virtual void bookHistograms();
  template <typename T>
  TLorentzVector getP4(const T& obj);
  int getGenDauPgd(const vhtm::GenParticle gp);
  double getHmassfromZdau(const vhtm::GenParticle Z1,const vhtm::GenParticle Z2);
private:
  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::GenParticle> genZ;  
  bool _dumpEvent;
};
#endif
