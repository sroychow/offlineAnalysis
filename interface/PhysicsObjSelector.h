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

class PhysicsObjSelector : public AnaBase{
  public :
    PhysicsObjSelector();
    virtual void bookHistograms();
    
    std::vector<vhtm::Muon>& getLooseMu(){ return loosemuVec_;}
    std::vector<vhtm::Muon>& getTightMu(){ return tightmuVec_;}
    std::vector<vhtm::Electron>& getLooseEle(){ return looseeleVec_;}
    std::vector<vhtm::Electron>& getTightEle(){ return tighteleVec_;}
    std::vector<vhtm::PackedPFCandidate>& getPFPhotonFSR(){ return pfphotonVec_; }
    std::vector< std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> >>& 
                                                        getlooseElePhoFSR() { return looseElePhoFSRpair_; }
    std::vector< std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > >& 
                                                        getlooseMuPhoFSR() { return looseMuPhoFSRpair_;}
    void muonSelector();
    void electronSelector();
    void pfphotonSelector();
    void clear();
  private :
    bool _dumpEvent;
    std::vector<vhtm::Muon> loosemuVec_,tightmuVec_;
    std::vector<vhtm::Electron> looseeleVec_,tighteleVec_;
    std::vector<vhtm::PackedPFCandidate> pfphotonVec_;
    std::vector< std::pair< vhtm::Electron, std::vector<vhtm::PackedPFCandidate> >> looseElePhoFSRpair_; 
    std::vector< std::pair< vhtm::Muon, std::vector<vhtm::PackedPFCandidate> > > looseMuPhoFSRpair_; 
    bool searchedEle_,searchedMu_,searchedPhoton_;
 };
#endif
