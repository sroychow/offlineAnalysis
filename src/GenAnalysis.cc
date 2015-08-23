#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <functional>
#include <numeric>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility> 
#include <typeinfo>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "AnaUtil.h"
#include "HZZ4lUtil.h"
#include "GenAnalysis.h"

GenAnalysis::GenAnalysis( std::vector<vhtm::GenParticle> genPartList_, TFile* fout, bool bookHisto ) {
  genPList_ = genPartList_;
  histf_ = fout;  
  genPass_ =  false;
  checkGenLevel(bookHisto);
}

GenAnalysis::~GenAnalysis() {
}

void GenAnalysis::bookHistograms() {
  histf_->cd();
  histf_->mkdir("GenLevel");
  histf_->cd("GenLevel");
  new TH1D("h1decayFlav", "h1 decay flavour", 4, 0.5, 4.5);
  new TH1D("genHmass","mass os the H in generator level",200,0,200);
  histf_->cd();
}

void GenAnalysis::checkGenLevel(bool bookH) {
  if(bookH) bookHistograms();
  histf_->cd();
  histf_->cd("GenLevel");
  genPass_ = false;
  //vector<vhtm::GenParticle> genZList_;
  //vector<vhtm::GenParticle> genHList_;
  for ( const auto& gp: genPList_ ) {
    int pdgid = std::abs(gp.pdgId);
    int status = gp.status;
    if (pdgid == 4000035) {
      std::vector<int> d = gp.daughtIndices;
      if(d.size()!=2) continue;
      const auto& dgp1 = genPList_.at(d.at(0));
      const auto& dgp2 = genPList_.at(d.at(1));
      int d_pdgid1 = std::abs(dgp1.pdgId);
      int d_pdgid2 = std::abs(dgp2.pdgId);
      std::vector<int>  h1dauPdg;
      h1dauPdg.push_back(d_pdgid1);
      h1dauPdg.push_back(d_pdgid2);
      if ( !( ( (d_pdgid1 == 11||d_pdgid1 == 13) && d_pdgid2 == 15) || 
                (d_pdgid1 == 15 && (d_pdgid2 == 11||d_pdgid2 == 13) ) ) ) 
      continue;
      //EventType eflav;
      if( std::find (h1dauPdg.begin(), h1dauPdg.end(), 15) != h1dauPdg.end() &&
          std::find (h1dauPdg.begin(), h1dauPdg.end(), 13) != h1dauPdg.end() )
        AnaUtil::fillHist1D("h1decayFlav", 1,1); 
      else if( std::find (h1dauPdg.begin(), h1dauPdg.end(), 15) != h1dauPdg.end() &&
          std::find (h1dauPdg.begin(), h1dauPdg.end(), 11) != h1dauPdg.end() )
        AnaUtil::fillHist1D("h1decayFlav", 2,1); 
      genPass_ = true;
      //genHList_.push_bak(gp); 
      TLorentzVector l = HZZ4lUtil::getP4(gp);
      AnaUtil::fillHist1D("genHmass", l.M(),1);
    }
  }
}







