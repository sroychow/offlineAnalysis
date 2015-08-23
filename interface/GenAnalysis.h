#ifndef __GenAnalysis__hh
#define __GenAnalysis__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>

#include<TFile.h>

#include "PhysicsObjects.h"

class GenAnalysis  {
    
public:
  enum EventType {
    mmmt = 0, mmet, eemt, eeet
  };
  GenAnalysis(std::vector<vhtm::GenParticle> genPartList_, TFile* fout, bool bookHisto);
  ~GenAnalysis();

  void bookHistograms();
  
  bool genOk(){ return genPass_; }
  void checkGenLevel( bool bookH );
    
public:
  std::vector<vhtm::GenParticle> genPList_;
  TFile* histf_;
  bool genPass_;
  std::string dumpGenFilename_;
};
#endif
