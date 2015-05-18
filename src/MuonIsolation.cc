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

#include "MuonIsolation.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

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
MuonIsolation::MuonIsolation()
  : AnaBase(),
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
MuonIsolation::~MuonIsolation() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MuonIsolation::beginJob() 
{ 
  AnaBase::beginJob();

  // Open the output ROOT file
  histf()->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void MuonIsolation::bookHistograms() 
{
 new TH1D("nGoodvtx","Number of Vertices",50,0,50);
 new TH1D("nRecomuon","Number of Reconstructed Muons(without cut)",20,0,20);
 new TH1D("nGoodmuon","Number of Good Muons(with selection cuts)",20,0,20);
 new TH1F("muon1Pt", "Pt of muon", 150, 0., 150.);
 new TProfile("muPtvsnVtx", "Muon Pt vs nVertex", 80, 0, 80, 0.,150.);
// new TH1F("muon1Pt", "Pt of sub-leading muon", 150, 0., 150.);
 new TH1F("muon1Eta", "Eta of muon", 30, -3., 3.);
//Rho Values
 new TH1F("fGridRhoAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAllCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralChargedPileUp","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralNeutral","Grid Rho for event",100,0.,100.);
//Charged Isolation
 new TH1F("muonChIso_015","Muon Charged Isolation for cone size 0.15",100,0.,10.);
 new TH1F("muonChIso_020","Muon Charged Isolation for cone size 0.20",100,0.,10.);
 new TH1F("muonChIso_025","Muon Charged Isolation for cone size 0.25",100,0.,10.);
 new TH1F("muonChIso_030","Muon Charged Isolation for cone size 0.30",100,0.,10.);
 new TH1F("muonChIso_035","Muon Charged Isolation for cone size 0.35",100,0.,10.);
 new TH1F("muonChIso_040","Muon Charged Isolation for cone size 0.40",100,0.,10.);
 new TH1F("muonChIso_045","Muon Charged Isolation for cone size 0.45",100,0.,10.);
//Neutral Isolation
 new TH1F("muonNIso_015","Muon Neutral Isolation for cone size 0.15",100,0.,10.);
 new TH1F("muonNIso_020","Muon Neutral Isolation for cone size 0.20",100,0.,10.);
 new TH1F("muonNIso_025","Muon Neutral Isolation for cone size 0.25",100,0.,10.);
 new TH1F("muonNIso_030","Muon Neutral Isolation for cone size 0.30",100,0.,10.);
 new TH1F("muonNIso_035","Muon Neutral Isolation for cone size 0.35",100,0.,10.);
 new TH1F("muonNIso_040","Muon Neutral Isolation for cone size 0.40",100,0.,10.);
 new TH1F("muonNIso_045","Muon Neutral Isolation for cone size 0.45",100,0.,10.);
//Pileup 
 new TH1F("muonPuIso_015","Muon Isolation from Pileup for cone size 0.15",100,0.,10.);
 new TH1F("muonPuIso_020","Muon Isolation from Pileup for cone size 0.20",100,0.,10.);
 new TH1F("muonPuIso_025","Muon Isolation from Pileup for cone size 0.25",100,0.,10.);
 new TH1F("muonPuIso_030","Muon Isolation from Pileup for cone size 0.30",100,0.,10.);
 new TH1F("muonPuIso_035","Muon Isolation from Pileup for cone size 0.35",100,0.,10.);
 new TH1F("muonPuIso_040","Muon Isolation from Pileup for cone size 0.40",100,0.,10.);
 new TH1F("muonPuIso_045","Muon Isolation from Pileup for cone size 0.45",100,0.,10.);
//Charged Isolation vs Rho
 new TProfile("chIso015vsrho", "Charged Isolation 015 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso020vsrho", "Charged Isolation 020 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso025vsrho", "Charged Isolation 025 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso030vsrho", "Charged Isolation 030 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso035vsrho", "Charged Isolation 035 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso040vsrho", "Charged Isolation 040 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
 new TProfile("chIso045vsrho", "Charged Isolation 045 vs fixedGridRhoFastjetAll", 100, 0., 100., 0.,10.);
//Isolation from vs Rho
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
//Rho Corrected Isolation vs NVertex
 new TProfile("rhIso015vsnVtx", "Rho Corrected Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIso020vsnVtx", "Rho Corrected Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIso025vsnVtx", "Rho Corrected Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIso030vsnVtx", "Rho Corrected Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIso035vsnVtx", "Rho Corrected Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIso040vsnVtx", "Rho Corrected Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("rhIso045vsnVtx", "Rho Corrected Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
//Delta Beta Corrected Isolation vs NVertex
 new TProfile("dbIsoraw015vsnVtx", "Delta-beta Corrected Raw Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoraw020vsnVtx", "Delta-beta Corrected Raw Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoraw025vsnVtx", "Delta-beta Corrected Raw Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoraw030vsnVtx", "Delta-beta Corrected Raw Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoraw035vsnVtx", "Delta-beta Corrected Raw Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoraw040vsnVtx", "Delta-beta Corrected Raw Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsoraw045vsnVtx", "Delta-beta Corrected Raw Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
//Quadratic Delta-Beta Corrected Isolation vs NVertex
 new TProfile("QdbIsoraw015vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("QdbIsoraw020vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("QdbIsoraw025vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("QdbIsoraw030vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("QdbIsoraw035vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("QdbIsoraw040vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("QdbIsoraw045vsnVtx", "Quadratic Delta-beta Corrected Raw Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
//Quadratic Delta-Beta Corrected Isolation vs muon pt
 new TProfile("QdbIsoraw015vspt", "Quadratic Delta-beta Corrected Relative Isolation 015 vs muonPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw020vspt", "Quadratic Delta-beta Corrected Relative Isolation 020 vs muonPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw025vspt", "Quadratic Delta-beta Corrected Relative Isolation 025 vs muonPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw030vspt", "Quadratic Delta-beta Corrected Relative Isolation 030 vs muonPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw035vspt", "Quadratic Delta-beta Corrected Relative Isolation 035 vs muonPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw040vspt", "Quadratic Delta-beta Corrected Relative Isolation 040 vs muonPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw045vspt", "Quadratic Delta-beta Corrected Relative Isolation 045 vs muonPt", 150, 0, 150, 0.,5.);
    //Default Delta Beta Rel Corrected Isolation vs NVertex
    new TProfile("OdbIsorel015vsnVtx", "Delta-beta Corrected Relative Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("OdbIsorel020vsnVtx", "Delta-beta Corrected Relative Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("OdbIsorel025vsnVtx", "Delta-beta Corrected Relative Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("OdbIsorel030vsnVtx", "Delta-beta Corrected Relative Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("OdbIsorel035vsnVtx", "Delta-beta Corrected Relative Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("OdbIsorel040vsnVtx", "Delta-beta Corrected Relative Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("OdbIsorel045vsnVtx", "Delta-beta Corrected Relative Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);

    //Delta Beta Rel Corrected Isolation vs NVertex
 new TProfile("dbIsorel015vsnVtx", "Delta-beta Corrected Relative Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsorel020vsnVtx", "Delta-beta Corrected Relative Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsorel025vsnVtx", "Delta-beta Corrected Relative Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsorel030vsnVtx", "Delta-beta Corrected Relative Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsorel035vsnVtx", "Delta-beta Corrected Relative Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsorel040vsnVtx", "Delta-beta Corrected Relative Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
 new TProfile("dbIsorel045vsnVtx", "Delta-beta Corrected Relative Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
    //Quadratic Delta Beta Rel Corrected Isolation vs NVertex
    new TProfile("QdbIsorel015vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("QdbIsorel020vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("QdbIsorel025vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("QdbIsorel030vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("QdbIsorel035vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("QdbIsorel040vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("QdbIsorel045vsnVtx", "Quadratic Delta-beta Corrected Relative Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
    //Rho Corrected Relative Isolation vs NVertex
    new TProfile("rhIsorel015vsnVtx", "Rho Corrected Relative Isolation 015 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("rhIsorel020vsnVtx", "Rho Corrected Relative Isolation 020 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("rhIsorel025vsnVtx", "Rho Corrected Relative Isolation 025 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("rhIsorel030vsnVtx", "Rho Corrected Relative Isolation 030 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("rhIsorel035vsnVtx", "Rho Corrected Relative Isolation 035 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("rhIsorel040vsnVtx", "Rho Corrected Relative Isolation 040 vs nVertex", 80, 0, 80, 0.,5.);
    new TProfile("rhIsorel045vsnVtx", "Rho Corrected Relative Isolation 045 vs nVertex", 80, 0, 80, 0.,5.);
    //Relative DB Isolation original
     new TH1F("OdbmuonrelIso_015","Muon Relative Delta-Beta Isolation for cone size 0.15",100,0.,1.);
     new TH1F("OdbmuonrelIso_020","Muon Relative Delta-Beta Isolation for cone size 0.20",100,0.,1.);
     new TH1F("OdbmuonrelIso_025","Muon Relative Delta-Beta Isolation for cone size 0.25",100,0.,1.);
     new TH1F("OdbmuonrelIso_030","Muon Relative Delta-Beta Isolation for cone size 0.30",100,0.,1.);
     new TH1F("OdbmuonrelIso_035","Muon Relative Delta-Beta Isolation for cone size 0.35",100,0.,1.);
     new TH1F("OdbmuonrelIso_040","Muon Relative Delta-Beta Isolation for cone size 0.40",100,0.,1.);
     new TH1F("OdbmuonrelIso_045","Muon Relative Delta-Beta Isolation for cone size 0.45",100,0.,1.);
    //Relative DB Isolation original
    new TH1F("QdbmuonrelIso_015","Muon Relative Delta-Beta Isolation for cone size 0.15",100,0.,1.);
    new TH1F("QdbmuonrelIso_020","Muon Relative Delta-Beta Isolation for cone size 0.20",100,0.,1.);
    new TH1F("QdbmuonrelIso_025","Muon Relative Delta-Beta Isolation for cone size 0.25",100,0.,1.);
    new TH1F("QdbmuonrelIso_030","Muon Relative Delta-Beta Isolation for cone size 0.30",100,0.,1.);
    new TH1F("QdbmuonrelIso_035","Muon Relative Delta-Beta Isolation for cone size 0.35",100,0.,1.);
    new TH1F("QdbmuonrelIso_040","Muon Relative Delta-Beta Isolation for cone size 0.40",100,0.,1.);
    new TH1F("QdbmuonrelIso_045","Muon Relative Delta-Beta Isolation for cone size 0.45",100,0.,1.);
    //Relative DB Isolation original
    new TH1F("RhmuonrelIso_015","Muon Relative Rho Isolation for cone size 0.15",100,0.,1.);
    new TH1F("RhmuonrelIso_020","Muon Relative Rho Isolation for cone size 0.20",100,0.,1.);
    new TH1F("RhmuonrelIso_025","Muon Relative Rho Isolation for cone size 0.25",100,0.,1.);
    new TH1F("RhmuonrelIso_030","Muon Relative Rho Isolation for cone size 0.30",100,0.,1.);
    new TH1F("RhmuonrelIso_035","Muon Relative Rho Isolation for cone size 0.35",100,0.,1.);
    new TH1F("RhmuonrelIso_040","Muon Relative Rho Isolation for cone size 0.40",100,0.,1.);
    new TH1F("RhmuonrelIso_045","Muon Relative Rho Isolation for cone size 0.45",100,0.,1.);
    //DB factor
    new TH1F("QdbFac_015","Delta-Beta for cone size 0.15",200,0.15,.6);
    new TH1F("QdbFac_020","Delta-Beta for cone size 0.20",200,0.15,.6);
    new TH1F("QdbFac_025","Delta-Beta for cone size 0.25",200,0.15,.6);
    new TH1F("QdbFac_030","Delta-Beta for cone size 0.30",200,0.15,.6);
    new TH1F("QdbFac_035","Delta-Beta for cone size 0.35",200,0.15,.6);
    new TH1F("QdbFac_040","Delta-Beta for cone size 0.40",200,0.15,.6);
    new TH1F("QdbFac_045","Delta-Beta for cone size 0.45",200,0.15,.6);
    //Rho factor
    new TH1F("QrhoFac_015","Rho for cone size 0.15",200,0.,.6);
    new TH1F("QrhoFac_020","Rho for cone size 0.20",200,0.,.6);
    new TH1F("QrhoFac_025","Rho for cone size 0.25",200,0.,.6);
    new TH1F("QrhoFac_030","Rho for cone size 0.30",200,0.,.6);
    new TH1F("QrhoFac_035","Rho for cone size 0.35",200,0.,.6);
    new TH1F("QrhoFac_040","Rho for cone size 0.40",200,0.,.6);
    new TH1F("QrhoFac_045","Rho for cone size 0.45",200,0.,.6);

}
//----------------------------------------------
//Fill Histo for a cone
//----------------------------------------------
void MuonIsolation::fillHistoforCone(std::string c,double mupt,double ch,double nu,double pu,int ngoodVtx,double fGridRhoFastjetAll) {
 if(cone_.find(c)!=cone_.end()) {
   AnaUtil::fillHist1D("muonChIso_"+cone_[c],ch,1);
   AnaUtil::fillHist1D("muonNIso_"+cone_[c],nu,1);
   AnaUtil::fillHist1D("muonPuIso_"+cone_[c],pu,1);
   //do default delta-beta
     double Odbiso = ch + std::max(0.0, nu-0.5*pu);
     AnaUtil::fillHist1D("OdbmuonrelIso_"+cone_[c],Odbiso/mupt,1);
   // do deltaBeta
   double dbiso = ch + std::max(0.0, nu-AnaUtil::cutValue(muonCutMap(), "db"+cone_[c])*pu-AnaUtil::cutValue(muonCutMap(), "dbInt"+cone_[c]));
   //do deltabeta from quadratic fit
   double apu=AnaUtil::cutValue(muonCutMap(), "apuIso"+cone_[c]);
   double aN=AnaUtil::cutValue(muonCutMap(), "aNIso"+cone_[c]);
   double bpu=AnaUtil::cutValue(muonCutMap(), "bpuIso"+cone_[c]);
   double bN=AnaUtil::cutValue(muonCutMap(), "bNIso"+cone_[c]);
   double cpu=AnaUtil::cutValue(muonCutMap(), "cpuIso"+cone_[c]);
   double cN=AnaUtil::cutValue(muonCutMap(), "cNIso"+cone_[c]);

   double dbFac= (bN + 2.*cN*ngoodVtx)/(bpu + 2.*cpu*ngoodVtx);
   AnaUtil::fillHist1D("QdbFac_"+cone_[c],dbFac,1);
     double Qdbiso= ch + std::max(0.0, nu - dbFac*pu - aN);
   AnaUtil::fillHist1D("QdbmuonrelIso_"+cone_[c],Qdbiso/mupt,1);
   //do rho Correction
   double rhoSlope=AnaUtil::cutValue(muonCutMap(),"brho"+cone_[c]) + 2.*AnaUtil::cutValue(muonCutMap(),"crho"+cone_[c])*fGridRhoFastjetAll;
   AnaUtil::fillHist1D("QrhoFac_"+cone_[c],rhoSlope,1);
//   double rhoIso=ch + std::max(0.0, nu - rhoSlope*fGridRhoFastjetAll - AnaUtil::cutValue(muonCutMap(),"arho"+cone_[c]));
   double rhoIso=ch + std::max(0.0, nu - rhoSlope*fGridRhoFastjetAll);
   AnaUtil::fillHist1D("RhmuonrelIso_"+cone_[c],rhoIso/mupt,1);
   //AnaUtil::fillProfile("dbIso"+cone_[c]+"vsnVtx",ngoodVtx,dbiso);
   //AnaUtil::fillHist1D("muonrelIso_"+cone_[c],iso/mupt,1);
   //vs Rho
   AnaUtil::fillProfile("chIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,ch);
   AnaUtil::fillProfile("NIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,nu);
   AnaUtil::fillProfile("puIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,pu);

   //AnaUtil::fillProfile("rhIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,rhoIso);
   //AnaUtil::fillHist1D("muonRhorelIso_"+cone_[c],rhoIso/mupt,1);   
   //Vs NVertex
   AnaUtil::fillProfile("chIso"+cone_[c]+"vsnVtx",ngoodVtx,ch);
   AnaUtil::fillProfile("NIso"+cone_[c]+"vsnVtx",ngoodVtx,nu);
   AnaUtil::fillProfile("puIso"+cone_[c]+"vsnVtx",ngoodVtx,pu);    
   AnaUtil::fillProfile("rhIso"+cone_[c]+"vsnVtx",ngoodVtx,rhoIso);
   AnaUtil::fillProfile("dbIsoraw"+cone_[c]+"vsnVtx",ngoodVtx,dbiso); 
   AnaUtil::fillProfile("QdbIsoraw"+cone_[c]+"vsnVtx",ngoodVtx,Qdbiso); 
   AnaUtil::fillProfile("dbIsorel"+cone_[c]+"vsnVtx",ngoodVtx,dbiso/mupt);
   AnaUtil::fillProfile("QdbIsorel"+cone_[c]+"vsnVtx",ngoodVtx,Qdbiso/mupt);
   AnaUtil::fillProfile("OdbIsorel"+cone_[c]+"vsnVtx",ngoodVtx,Odbiso/mupt);
   AnaUtil::fillProfile("rhIsorel"+cone_[c]+"vsnVtx",ngoodVtx,rhoIso/mupt);

   //vsPt
   AnaUtil::fillProfile("QdbIsoraw"+cone_[c]+"vspt",mupt,Qdbiso/mupt); 
 }
 else
  std::cout<<"Wrong cone value"<<std::endl;
}

// -------------------
// The main event loop
// -------------------
void MuonIsolation::clearLists() {
 vtxList.clear();
}
void MuonIsolation::eventLoop() 
{
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
   double ngoodVtx=vtxList.size();
   AnaUtil::fillHist1D("nGoodvtx",ngoodVtx,1);
   // AnaUtil::fillHist1D("evcounter", 0, _puevWt);
   AnaUtil::fillHist1D("fGridRhoAll",evt.fGridRhoAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAll",evt.fGridRhoFastjetAll,1);
   AnaUtil::fillHist1D("fGridRhoFastjetAllCalo",evt.fGridRhoFastjetAllCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralCalo",evt.fGridRhoFastjetCentralCalo,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralChargedPileUp",evt.fGridRhoFastjetCentralChargedPileUp,1);
   AnaUtil::fillHist1D("fGridRhoFastjetCentralNeutral",evt.fGridRhoFastjetCentralNeutral,1);
   if(!muonColl()->empty())
     AnaUtil::fillHist1D("nRecomuon",muonColl()->size(),1);
   int ngoodmu=0;
   for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
     const Muon& muon = (*it);
     if (abs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta"))                      continue;
     if (muon.pt < AnaUtil::cutValue(muonCutMap(), "pt"))                            continue;
     if (!muon.isTrackerMuon)                                                         continue;
     if (!muon.isGlobalMuonPromptTight)                                               continue;
     if (muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers"))                continue;
     if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches"))                  continue;
     if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations"))  continue; 
     if (muon.pixHits < AnaUtil::cutValue(muonCutMap(), "pixHits"))                    continue;
     if (muon.trkHits < AnaUtil::cutValue(muonCutMap(), "trkHits"))                    continue;
     if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2"))              continue;
     if (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0"))                   continue;
     ngoodmu++;
     //for(auto &mm : muon.isolationMap) {
     // std::cout<<mm.first << std::endl;
     // }

     AnaUtil::fillHist1D("muon1Pt",muon.pt,1);
     AnaUtil::fillHist1D("muon1Eta",muon.eta,1);
     
     AnaUtil::fillProfile("muPtvsnVtx",ngoodVtx,muon.pt);
     if(muon.isolationMap.find("c15") != muon.isolationMap.end() ) {
       fillHistoforCone("c15",muon.pt,muon.isolationMap.at("c15").at(0),muon.isolationMap.at("c15").at(1),muon.isolationMap.at("c15").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
     }
     if(muon.isolationMap.find("c20") != muon.isolationMap.end() ) {
       fillHistoforCone("c20",muon.pt,muon.isolationMap.at("c20").at(0),muon.isolationMap.at("c20").at(1),
                         muon.isolationMap.at("c20").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
     }
     if(muon.isolationMap.find("c25") != muon.isolationMap.end() ) {
       fillHistoforCone("c25",muon.pt,muon.isolationMap.at("c25").at(0),muon.isolationMap.at("c25").at(1),
                         muon.isolationMap.at("c25").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
     }
     if(muon.isolationMap.find("c30") != muon.isolationMap.end() ) {
       fillHistoforCone("c30",muon.pt,muon.isolationMap.at("c30").at(0),muon.isolationMap.at("c30").at(1),
                         muon.isolationMap.at("c30").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
     }
     if(muon.isolationMap.find("c35") != muon.isolationMap.end() ) {
       fillHistoforCone("c35",muon.pt,muon.isolationMap.at("c35").at(0),muon.isolationMap.at("c35").at(1),
                         muon.isolationMap.at("c35").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
     }
     if(muon.isolationMap.find("c40") != muon.isolationMap.end() ) {
       fillHistoforCone("c40",muon.pt,muon.isolationMap.at("c40").at(0),muon.isolationMap.at("c40").at(1),
                         muon.isolationMap.at("c40").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
    }
    if(muon.isolationMap.find("c45") != muon.isolationMap.end() ) {
       fillHistoforCone("c45",muon.pt,muon.isolationMap.at("c45").at(0),muon.isolationMap.at("c45").at(1),
                         muon.isolationMap.at("c45").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
    }
   }
   if(ngoodmu > 0)
     AnaUtil::fillHist1D("nGoodmuon",ngoodmu,1);
  }  
  // Analysis is over
  endJob();
}

void MuonIsolation::endJob() {

  closeFiles();

  histf()->cd();
  histf()->Write();
  histf()->Close();
  delete histf();
}
