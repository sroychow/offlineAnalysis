#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "HZZ4lUtil.h"
#include "AnaUtil.h"
#include "LeptonIsolation.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
LeptonIsolation::LeptonIsolation()
  : PhysicsObjSelector(),
    _dumpEvent(false)
{
 cone_["c15"]="015";
 cone_["c20"]="020";
 cone_["c25"]="025";
 cone_["c30"]="030";
 cone_["c35"]="035";
 cone_["c40"]="040";
 cone_["c45"]="045";
}
// ----------
// Destructor
// ----------
LeptonIsolation::~LeptonIsolation() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool LeptonIsolation::beginJob() 
{ 
  PhysicsObjSelector::beginJob();

  // Open the output ROOT file
  histf()->cd();
  bookHistograms();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void LeptonIsolation::bookHistograms() 
{
 histf()->cd();
 histf()->mkdir("Event");
 new TH1D("nGoodvtx","Number of Vertices",50,0,50);
//Rho Values
 new TH1F("fGridRhoAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAllCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralChargedPileUp","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralNeutral","Grid Rho for event",100,0.,100.);
 bookHistograms("electron");
 bookHistograms("muon");
}
//-------------------------------------
//Book isolation histograms for leptons
//-------------------------------------
void LeptonIsolation::bookHistograms(TString folder) 
{
 histf()->cd();
 histf()->mkdir(folder);
//Kinematic info
 new TH1F("pt", "Lepton Pt", 150, 0., 150.);
 new TProfile("ptvsnVtx", "Lepton Pt vs nVertex", 80, 0, 80, 0.,150.);
 new TH1F("eta", "Lepton Eta", 30, -3., 3.);
//Charged Isolation
 new TH1F("ChIso_015","Charged Isolation for cone size 0.15",100,0.,10.);
 new TH1F("ChIso_020","Charged Isolation for cone size 0.20",100,0.,10.);
 new TH1F("ChIso_025","Charged Isolation for cone size 0.25",100,0.,10.);
 new TH1F("ChIso_030","Charged Isolation for cone size 0.30",100,0.,10.);
 new TH1F("ChIso_035","Charged Isolation for cone size 0.35",100,0.,10.);
 new TH1F("ChIso_040","Charged Isolation for cone size 0.40",100,0.,10.);
 new TH1F("ChIso_045","Charged Isolation for cone size 0.45",100,0.,10.);
//Neutral Isolation
 new TH1F("NIso_015","Neutral Isolation for cone size 0.15",100,0.,10.);
 new TH1F("NIso_020","Neutral Isolation for cone size 0.20",100,0.,10.);
 new TH1F("NIso_025","Neutral Isolation for cone size 0.25",100,0.,10.);
 new TH1F("NIso_030","Neutral Isolation for cone size 0.30",100,0.,10.);
 new TH1F("NIso_035","Neutral Isolation for cone size 0.35",100,0.,10.);
 new TH1F("NIso_040","Neutral Isolation for cone size 0.40",100,0.,10.);
 new TH1F("NIso_045","Neutral Isolation for cone size 0.45",100,0.,10.);
//Pileup 
 new TH1F("PuIso_015","Isolation from Pileup for cone size 0.15",100,0.,10.);
 new TH1F("PuIso_020","Isolation from Pileup for cone size 0.20",100,0.,10.);
 new TH1F("PuIso_025","Isolation from Pileup for cone size 0.25",100,0.,10.);
 new TH1F("PuIso_030","Isolation from Pileup for cone size 0.30",100,0.,10.);
 new TH1F("PuIso_035","Isolation from Pileup for cone size 0.35",100,0.,10.);
 new TH1F("PuIso_040","Isolation from Pileup for cone size 0.40",100,0.,10.);
 new TH1F("PuIso_045","Isolation from Pileup for cone size 0.45",100,0.,10.);
//Charged Isolation vs Rho
 new TProfile("chIso015vsrho", "Charged Isolation 015 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso020vsrho", "Charged Isolation 020 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso025vsrho", "Charged Isolation 025 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso030vsrho", "Charged Isolation 030 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso035vsrho", "Charged Isolation 035 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso040vsrho", "Charged Isolation 040 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso045vsrho", "Charged Isolation 045 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
//Isolation from PU vs Rho
 new TProfile("puIso015vsrho", "Pileup Isolation 015 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("puIso020vsrho", "Pileup Isolation 020 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("puIso025vsrho", "Pileup Isolation 025 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("puIso030vsrho", "Pileup Isolation 030 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("puIso035vsrho", "Pileup Isolation 035 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("puIso040vsrho", "Pileup Isolation 040 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.); 
 new TProfile("puIso045vsrho", "Pileup Isolation 045 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
//Neutral Isolation vs Rho
 new TProfile("NIso015vsrho", "Neutral Isolation 015 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso020vsrho", "Neutral Isolation 020 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso025vsrho", "Neutral Isolation 025 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso030vsrho", "Neutral Isolation 030 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso035vsrho", "Neutral Isolation 035 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso040vsrho", "Neutral Isolation 040 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("NIso045vsrho", "Neutral Isolation 045 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
/////////////////////////////////////////////////Isolation vs NVertex//////////////////////////////////////
//Charged Isolation vs NVertex
 new TProfile("chIso015vsnVtx", "Charged Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("chIso020vsnVtx", "Charged Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("chIso025vsnVtx", "Charged Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("chIso030vsnVtx", "Charged Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("chIso035vsnVtx", "Charged Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("chIso040vsnVtx", "Charged Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("chIso045vsnVtx", "Charged Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
//Isolation from PU vs NVertex
 new TProfile("puIso015vsnVtx", "Pileup Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("puIso020vsnVtx", "Pileup Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("puIso025vsnVtx", "Pileup Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("puIso030vsnVtx", "Pileup Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("puIso035vsnVtx", "Pileup Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("puIso040vsnVtx", "Pileup Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("puIso045vsnVtx", "Pileup Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
//Neutral Isolation vs NVertex
 new TProfile("NIso015vsnVtx", "Neutral Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("NIso020vsnVtx", "Neutral Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("NIso025vsnVtx", "Neutral Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("NIso030vsnVtx", "Neutral Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("NIso035vsnVtx", "Neutral Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("NIso040vsnVtx", "Neutral Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("NIso045vsnVtx", "Neutral Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
//Rho Corrected Raw Isolation vs NVertex
 new TProfile("rhIsoRaw015vsnVtx", "Rho Corrected Raw Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRaw020vsnVtx", "Rho Corrected Raw Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRaw025vsnVtx", "Rho Corrected Raw Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRaw030vsnVtx", "Rho Corrected Raw Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRaw035vsnVtx", "Rho Corrected Raw Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRaw040vsnVtx", "Rho Corrected Raw Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRaw045vsnVtx", "Rho Corrected Raw Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);

//Rho Corrected Relative Isolation vs NVertex
 new TProfile("rhIsoRel015vsnVtx", "Rho Corrected Relative Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRel020vsnVtx", "Rho Corrected Relative Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRel025vsnVtx", "Rho Corrected Relative Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRel030vsnVtx", "Rho Corrected Relative Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRel035vsnVtx", "Rho Corrected Relative Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRel040vsnVtx", "Rho Corrected Relative Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIsoRel045vsnVtx", "Rho Corrected Relative Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);

//Raw Delta Beta Corrected Isolation vs NVertex
 new TProfile("dbIsoRaw015vsnVtx", "Delta-beta Corrected Raw Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRaw020vsnVtx", "Delta-beta Corrected Raw Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRaw025vsnVtx", "Delta-beta Corrected Raw Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRaw030vsnVtx", "Delta-beta Corrected Raw Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRaw035vsnVtx", "Delta-beta Corrected Raw Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRaw040vsnVtx", "Delta-beta Corrected Raw Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRaw045vsnVtx", "Delta-beta Corrected Raw Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);

//Delta Beta Rel Corrected Isolation vs NVertex
 new TProfile("dbIsoRel015vsnVtx", "Delta-beta Corrected Relative Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRel020vsnVtx", "Delta-beta Corrected Relative Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRel025vsnVtx", "Delta-beta Corrected Relative Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRel030vsnVtx", "Delta-beta Corrected Relative Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRel035vsnVtx", "Delta-beta Corrected Relative Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRel040vsnVtx", "Delta-beta Corrected Relative Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoRel045vsnVtx", "Delta-beta Corrected Relative Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);

 //rho counter
 new TH1F("lep_rho_015", "no of leptons passing rho cut for cone size 0.15", 7, -0.5, 6.5);
 new TH1F("lep_rho_020", "no of leptons passing rho cut for cone size 0.20", 7, -0.5, 6.5);
 new TH1F("lep_rho_025", "no of leptons passing rho cut for cone size 0.25", 7, -0.5, 6.5);
 new TH1F("lep_rho_030", "no of leptons passing rho cut for cone size 0.30", 7, -0.5, 6.5);
 new TH1F("lep_rho_035", "no of leptons passing rho cut for cone size 0.35", 7, -0.5, 6.5);
 new TH1F("lep_rho_040", "no of leptons passing rho cut for cone size 0.40", 7, -0.5, 6.5);
 new TH1F("lep_rho_045", "no of leptons passing rho cut for cone size 0.45", 7, -0.5, 6.5);
 //Delta-beta counter
 new TH1F("lep_db_015", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
 new TH1F("lep_db_020", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
 new TH1F("lep_db_025", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
 new TH1F("lep_db_030", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
 new TH1F("lep_db_035", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
 new TH1F("lep_db_040", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
 new TH1F("lep_db_045", "no of leptons passing db cut for cone size 0.15",7,-0.5,6.5);
}
//---------------------------------
// Read Slope from NIso vs Rho fits
//---------------------------------
void LeptonIsolation::readRhofitSlope() {
}
//----------------------------------------------
//Fill Histo for a cone
//----------------------------------------------
void LeptonIsolation::fillHistoforCone(std::string c,double mupt,double ch_had,double ch_lep,double nu_had,double nu_photon,double pu,int ngoodVtx,double fGridRhoFastjetAll, TString& leptype) {
 histf()->cd();
 histf()->cd(leptype);
 if(cone_.find(c)!=cone_.end()) {
   AnaUtil::fillHist1D("ChIso_"+cone_[c],ch_had,1);
   AnaUtil::fillHist1D("NIso_"+cone_[c],nu_had+nu_photon,1);
   AnaUtil::fillHist1D("PuIso_"+cone_[c],pu,1);

   // do deltaBeta
   double dbiso = ch_had + std::max(0.0, nu_had+nu_photon-AnaUtil::cutValue(muonCutMap(), "db"+cone_[c])*pu-AnaUtil::cutValue(muonCutMap(), "dbInt"+cone_[c]));

   double rhoSlope=AnaUtil::cutValue(muonCutMap(),"brho"+cone_[c]); 

   double rhoIso=ch_had + std::max(0.0, nu_had+nu_photon - rhoSlope*fGridRhoFastjetAll);

   //isolation vs rho
   AnaUtil::fillProfile("chIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,ch_had);
   AnaUtil::fillProfile("NIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,nu_had+nu_photon);
   AnaUtil::fillProfile("puIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,pu);

   //Vs NVertex
   AnaUtil::fillProfile("chIso"+cone_[c]+"vsnVtx",ngoodVtx,ch_had);
   AnaUtil::fillProfile("NIso"+cone_[c]+"vsnVtx",ngoodVtx,nu_had+nu_photon);
   AnaUtil::fillProfile("puIso"+cone_[c]+"vsnVtx",ngoodVtx,pu);    
   AnaUtil::fillProfile("rhIsoRaw"+cone_[c]+"vsnVtx",ngoodVtx,rhoIso);
   AnaUtil::fillProfile("rhIsoRel"+cone_[c]+"vsnVtx",ngoodVtx,rhoIso/mupt);
   AnaUtil::fillProfile("dbIsoRaw"+cone_[c]+"vsnVtx",ngoodVtx,dbiso); 
   AnaUtil::fillProfile("dbIsoRel"+cone_[c]+"vsnVtx",ngoodVtx,dbiso/mupt); 

   AnaUtil::fillHist1D("lep_rho_"+cone_[c],0);
   ///Added parts for iso efficiency
   if( rhoIso/mupt < 0.5 ) {
     AnaUtil::fillHist1D("lep_rho_"+cone_[c],1);
     if( rhoIso/mupt < 0.45 ){
       AnaUtil::fillHist1D("lep_rho_"+cone_[c],2);
       if( rhoIso/mupt < 0.4 ){
         AnaUtil::fillHist1D("lep_rho_"+cone_[c],3);
         if( rhoIso/mupt < 0.35 ){
           AnaUtil::fillHist1D("lep_rho_"+cone_[c],4);
           if( rhoIso/mupt < 0.3 ){
             AnaUtil::fillHist1D("lep_rho_"+cone_[c],5);
             if( rhoIso/mupt < 0.25 ){ 
               AnaUtil::fillHist1D("lep_rho_"+cone_[c],6);
             }
           }  
         }
       }
     }
   }
   AnaUtil::fillHist1D("lep_db_"+cone_[c],0);
   if( dbiso/mupt < 0.5 ){
     AnaUtil::fillHist1D("lep_db_"+cone_[c],1);
     if( dbiso/mupt < 0.45 ){
       AnaUtil::fillHist1D("lep_db_"+cone_[c],2);
       if( dbiso/mupt < 0.4 ){
         AnaUtil::fillHist1D("lep_db_"+cone_[c],3);
         if( dbiso/mupt < 0.35 ){
           AnaUtil::fillHist1D("lep_db_"+cone_[c],4);
           if( dbiso/mupt < 0.3 ){
             AnaUtil::fillHist1D("lep_db_"+cone_[c],5);
             if( dbiso/mupt < 0.25 ){ 
               AnaUtil::fillHist1D("lep_db_"+cone_[c],6);
             }
           }  
         }
       }
     }
   }
 }
 else
  std::cout<<"Wrong cone value"<<std::endl;
}
//--------------------
//--------------------
template <typename T>
void LeptonIsolation::getLeptonIsolationInfo( std::vector<T> lepvec, double rho, int ngoodVtx, TString leptype ) {
  for (const auto& l : lepvec ) {
      double lpt = l.pt;
      AnaUtil::fillHist1D("pt",lpt,1);
        AnaUtil::fillHist1D("eta",l.eta,1);
        AnaUtil::fillProfile("ptvsnVtx",ngoodVtx,l.pt);
        if(l.isolationMap.find("c15") != l.isolationMap.end() ) {
            fillHistoforCone("c15",lpt,l.isolationMap.at("c15").at(0),
                             l.isolationMap.at("c15").at(1),
                             l.isolationMap.at("c15").at(2),
                             l.isolationMap.at("c15").at(3),
                             l.isolationMap.at("c15").at(4),
                             ngoodVtx,
                             rho,
                             leptype);
        }
        if(l.isolationMap.find("c20") != l.isolationMap.end() ) {
            fillHistoforCone("c20",lpt,l.isolationMap.at("c20").at(0),
                             l.isolationMap.at("c20").at(1),
                             l.isolationMap.at("c20").at(2),
                             l.isolationMap.at("c20").at(3),
                             l.isolationMap.at("c20").at(4),
                             ngoodVtx,
                             rho,
                             leptype);
        }
        if(l.isolationMap.find("c25") != l.isolationMap.end() ) {
            fillHistoforCone("c25",lpt,l.isolationMap.at("c25").at(0),
                             l.isolationMap.at("c25").at(1),
                             l.isolationMap.at("c25").at(2),
                             l.isolationMap.at("c25").at(3),
                             l.isolationMap.at("c25").at(4),
                             ngoodVtx,
                             rho,
                             leptype);
            
        }
        if(l.isolationMap.find("c30") != l.isolationMap.end() ) {
            fillHistoforCone("c30",lpt,l.isolationMap.at("c30").at(0),
                             l.isolationMap.at("c30").at(1),
                             l.isolationMap.at("c30").at(2),
                             l.isolationMap.at("c30").at(3),
                             l.isolationMap.at("c30").at(4),
                             ngoodVtx,
                             rho,
                             leptype);

        }
        if(l.isolationMap.find("c35") != l.isolationMap.end() ) {
            fillHistoforCone("c35",lpt,l.isolationMap.at("c35").at(0),
                             l.isolationMap.at("c35").at(1),
                             l.isolationMap.at("c35").at(2),
                             l.isolationMap.at("c35").at(3),
                             l.isolationMap.at("c35").at(4),
                             ngoodVtx,
                             rho,
                             leptype);

        }
        if(l.isolationMap.find("c40") != l.isolationMap.end() ) {
            fillHistoforCone("c40",lpt,l.isolationMap.at("c40").at(0),
                             l.isolationMap.at("c40").at(1),
                             l.isolationMap.at("c40").at(2),
                             l.isolationMap.at("c40").at(3),
                             l.isolationMap.at("c40").at(4),
                             ngoodVtx,
                             rho,
                             leptype);

        }
        if(l.isolationMap.find("c45") != l.isolationMap.end() ) {
            fillHistoforCone("c45",lpt,l.isolationMap.at("c45").at(0),
                             l.isolationMap.at("c45").at(1),
                             l.isolationMap.at("c45").at(2),
                             l.isolationMap.at("c45").at(3),
                             l.isolationMap.at("c45").at(4),
                             ngoodVtx,
                             rho,
                             leptype);

        }
    }
}



// -------------------
// The main event loop
// -------------------
void LeptonIsolation::clearLists() {
 vtxList.clear();
}
void LeptonIsolation::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = max(10000, nEvents()/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  //long int nRecoMuons=0,nGoodRecoMuons=0;
  std::cout<<"Bunch Crossing>>>>"<<bunchCrossing()<<std::endl;
  for (int ev = 0; ev < nEvents(); ++ev) {
    clearEvent();
    clearLists();
    int lflag = chain()->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName())); 

    const Event& evt = eventColl()->at(0);

    histf()->cd();

    //For data or for MC without pileup
    puevWt_ = 1;
    /*
    if (isMC()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
    */
    // Show status of the run
    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;

    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << chain()->GetTreeNumber()  
         << " ==> " << chain()->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;

   op.verbose = (logOption() >> 1 & 0x1); 
   findVtxInfo(vtxList, op, fLog());
   
   if( vtxList.empty() )             continue;
  
   int ngoodVtx=vtxList.size();
   histf()->cd();
   histf()->cd("Event"); 
   AnaUtil::fillHist1D("nGoodvtx",ngoodVtx,1);
   AnaUtil::fillHist1D("fGridRhoAll",evt.fGridRhoAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAll",evt.fGridRhoFastjetAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAllCalo",evt.fGridRhoFastjetAllCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralCalo",evt.fGridRhoFastjetCentralCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralChargedPileUp",evt.fGridRhoFastjetCentralChargedPileUp,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralNeutral",evt.fGridRhoFastjetCentralNeutral,1);

   const auto& loosemuVec_ = getLooseMuList();
   const auto& tightmuVec_ = getTightMuList();

   const auto& looseeleVec_ = getLooseEleList();
   const auto& tighteleVec_ = getTightEleList();

      if( !loosemuVec_.empty() ) getLeptonIsolationInfo<vhtm::Muon>(loosemuVec_, evt.fGridRhoFastjetAll,ngoodVtx,"muon");
      if( !looseeleVec_.empty() ) getLeptonIsolationInfo<vhtm::Electron>(looseeleVec_, evt.fGridRhoFastjetAll,ngoodVtx,"electron");
  }
  // Analysis is over
  endJob();
}

void LeptonIsolation::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}
