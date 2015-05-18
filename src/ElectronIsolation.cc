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

#include "ElectronIsolation.h"
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
 ElectronIsolation::ElectronIsolation()
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
ElectronIsolation::~ElectronIsolation() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ElectronIsolation::beginJob() 
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
void ElectronIsolation::bookHistograms()
{
 new TH1D("nGoodvtx","Number of Vertices",50,0,50);
 new TH1D("nRecoelectron","Number of Reconstructed Electrons(without cut)",20,0,20);
 new TH1D("nGoodelectron","Number of Good Electrons(with selection cuts)",20,0,20);
 new TH1F("electron1Pt", "Pt of electron", 150, 0., 150.);
 new TProfile("elecPtvsnVtx", "Electron Pt vs nVertex", 80, 0, 80, 0.,150.);
 // new TH1F("electron1Pt", "Pt of sub-leading electron", 150, 0., 150.);
 new TH1F("electron1Eta", "Eta of electron", 30, -3., 3.);
 //Rho Values
 new TH1F("fGridRhoAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAll","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetAllCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralCalo","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralChargedPileUp","Grid Rho for event",100,0.,100.);
 new TH1F("fGridRhoFastjetCentralNeutral","Grid Rho for event",100,0.,100.);
 //Charged Isolation
 new TH1F("electronChIso_015","Electron Charged Isolation for cone size 0.15",100,0.,10.);
 new TH1F("electronChIso_020","Electron Charged Isolation for cone size 0.20",100,0.,10.);
 new TH1F("electronChIso_025","Electron Charged Isolation for cone size 0.25",100,0.,10.);
 new TH1F("electronChIso_030","Electron Charged Isolation for cone size 0.30",100,0.,10.);
 new TH1F("electronChIso_035","Electron Charged Isolation for cone size 0.35",100,0.,10.);
 new TH1F("electronChIso_040","Electron Charged Isolation for cone size 0.40",100,0.,10.);
 new TH1F("electronChIso_045","Electron Charged Isolation for cone size 0.45",100,0.,10.);
 //Neutral Isolation
 new TH1F("electronNIso_015","Electron Neutral Isolation for cone size 0.15",100,0.,10.);
 new TH1F("electronNIso_020","Electron Neutral Isolation for cone size 0.20",100,0.,10.);
 new TH1F("electronNIso_025","Electron Neutral Isolation for cone size 0.25",100,0.,10.);
 new TH1F("electronNIso_030","Electron Neutral Isolation for cone size 0.30",100,0.,10.);
 new TH1F("electronNIso_035","Electron Neutral Isolation for cone size 0.35",100,0.,10.);
 new TH1F("electronNIso_040","Electron Neutral Isolation for cone size 0.40",100,0.,10.);
 new TH1F("electronNIso_045","Electron Neutral Isolation for cone size 0.45",100,0.,10.);
 //Pileup 
 new TH1F("electronPuIso_015","Electron Isolation from Pileup for cone size 0.15",100,0.,10.);
 new TH1F("electronPuIso_020","Electron Isolation from Pileup for cone size 0.20",100,0.,10.);
 new TH1F("electronPuIso_025","Electron Isolation from Pileup for cone size 0.25",100,0.,10.);
 new TH1F("electronPuIso_030","Electron Isolation from Pileup for cone size 0.30",100,0.,10.);
 new TH1F("electronPuIso_035","Electron Isolation from Pileup for cone size 0.35",100,0.,10.);
 new TH1F("electronPuIso_040","Electron Isolation from Pileup for cone size 0.40",100,0.,10.);
 new TH1F("electronPuIso_045","Electron Isolation from Pileup for cone size 0.45",100,0.,10.);
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
 //Quadratic Delta-Beta Corrected Isolation vs electron pt
 new TProfile("QdbIsoraw015vspt", "Quadratic Delta-beta Corrected Relative Isolation 015 vs electronPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw020vspt", "Quadratic Delta-beta Corrected Relative Isolation 020 vs electronPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw025vspt", "Quadratic Delta-beta Corrected Relative Isolation 025 vs electronPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw030vspt", "Quadratic Delta-beta Corrected Relative Isolation 030 vs electronPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw035vspt", "Quadratic Delta-beta Corrected Relative Isolation 035 vs electronPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw040vspt", "Quadratic Delta-beta Corrected Relative Isolation 040 vs electronPt", 150, 0, 150, 0.,5.);
 new TProfile("QdbIsoraw045vspt", "Quadratic Delta-beta Corrected Relative Isolation 045 vs electronPt", 150, 0, 150, 0.,5.);
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
 new TH1F("OdbelectronrelIso_015","Electron Relative Delta-Beta Isolation for cone size 0.15",100,0.,1.);
 new TH1F("OdbelectronrelIso_020","Electron Relative Delta-Beta Isolation for cone size 0.20",100,0.,1.);
 new TH1F("OdbelectronrelIso_025","Electron Relative Delta-Beta Isolation for cone size 0.25",100,0.,1.);
 new TH1F("OdbelectronrelIso_030","Electron Relative Delta-Beta Isolation for cone size 0.30",100,0.,1.);
 new TH1F("OdbelectronrelIso_035","Electron Relative Delta-Beta Isolation for cone size 0.35",100,0.,1.);
 new TH1F("OdbelectronrelIso_040","Electron Relative Delta-Beta Isolation for cone size 0.40",100,0.,1.);
 new TH1F("OdbelectronrelIso_045","Electron Relative Delta-Beta Isolation for cone size 0.45",100,0.,1.);
 //Relative DB Isolation original
 new TH1F("QdbelectronrelIso_015","Electron Relative Delta-Beta Isolation for cone size 0.15",100,0.,1.);
 new TH1F("QdbelectronrelIso_020","Electron Relative Delta-Beta Isolation for cone size 0.20",100,0.,1.);
 new TH1F("QdbelectronrelIso_025","Electron Relative Delta-Beta Isolation for cone size 0.25",100,0.,1.);
 new TH1F("QdbelectronrelIso_030","Electron Relative Delta-Beta Isolation for cone size 0.30",100,0.,1.);
 new TH1F("QdbelectronrelIso_035","Electron Relative Delta-Beta Isolation for cone size 0.35",100,0.,1.);
 new TH1F("QdbelectronrelIso_040","Electron Relative Delta-Beta Isolation for cone size 0.40",100,0.,1.);
 new TH1F("QdbelectronrelIso_045","Electron Relative Delta-Beta Isolation for cone size 0.45",100,0.,1.);
 //Relative DB Isolation original
 new TH1F("RhelectronrelIso_015","Electron Relative Rho Isolation for cone size 0.15",100,0.,1.);
 new TH1F("RhelectronrelIso_020","Electron Relative Rho Isolation for cone size 0.20",100,0.,1.);
 new TH1F("RhelectronrelIso_025","Electron Relative Rho Isolation for cone size 0.25",100,0.,1.);
 new TH1F("RhelectronrelIso_030","Electron Relative Rho Isolation for cone size 0.30",100,0.,1.);
 new TH1F("RhelectronrelIso_035","Electron Relative Rho Isolation for cone size 0.35",100,0.,1.);
 new TH1F("RhelectronrelIso_040","Electron Relative Rho Isolation for cone size 0.40",100,0.,1.);
 new TH1F("RhelectronrelIso_045","Electron Relative Rho Isolation for cone size 0.45",100,0.,1.);
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
//////Added these histograms for electron Isolation efficiency
    //db
    new TH1F("ePassdbIsorel_015","db corrected Iso(for e pass)for cone size 0.15",200,0.,.6);
    new TH1F("ePassdbIsorel_020","db corrected Iso(for e pass)for cone size 0.20",200,0.,.6);
    new TH1F("ePassdbIsorel_025","db corrected Iso(for e pass)for cone size 0.25",200,0.,.6);
    new TH1F("ePassdbIsorel_030","db corrected Iso(for e pass)for cone size 0.30",200,0.,.6);
    new TH1F("ePassdbIsorel_035","db corrected Iso(for e pass)for cone size 0.35",200,0.,.6);
    new TH1F("ePassdbIsorel_040","db corrected Iso(for e pass)for cone size 0.40",200,0.,.6);
    new TH1F("ePassdbIsorel_045","db corrected Iso(for e pass)for cone size 0.45",200,0.,.6);
    //rho
    new TH1F("ePassrhIsorel_015","rho corrected Iso(for e pass)for cone size 0.15",200,0.,.6);
    new TH1F("ePassrhIsorel_020","rho corrected Iso(for e pass)for cone size 0.20",200,0.,.6);
    new TH1F("ePassrhIsorel_025","rho corrected Iso(for e pass)for cone size 0.25",200,0.,.6);
    new TH1F("ePassrhIsorel_030","rho corrected Iso(for e pass)for cone size 0.30",200,0.,.6);
    new TH1F("ePassrhIsorel_035","rho corrected Iso(for e pass)for cone size 0.35",200,0.,.6);
    new TH1F("ePassrhIsorel_040","rho corrected Iso(for e pass)for cone size 0.40",200,0.,.6);
    new TH1F("ePassrhIsorel_045","rho corrected Iso(for e pass)for cone size 0.45",200,0.,.6);
    
    
}
//---------------------------------
// Read Slope from NIso vs Rho fits
//---------------------------------
void ElectronIsolation::readRhofitSlope() {
}
//----------------------------------------------
//Fill Histo for a cone
//----------------------------------------------
void ElectronIsolation::fillHistoforCone(std::string c,double elecpt,double ch,double nu,double pu,int ngoodVtx,double fGridRhoFastjetAll) {
 if(cone_.find(c)!=cone_.end()) {
  AnaUtil::fillHist1D("electronChIso_"+cone_[c],ch,1);
  AnaUtil::fillHist1D("electronNIso_"+cone_[c],nu,1);
  AnaUtil::fillHist1D("electronPuIso_"+cone_[c],pu,1);
  //do default delta-beta
  double Odbiso = ch + std::max(0.0, nu-0.5*pu);
  AnaUtil::fillHist1D("OdbelectronrelIso_"+cone_[c],Odbiso/elecpt,1);
  // do deltaBeta
  double dbiso = ch + std::max(0.0, nu-AnaUtil::cutValue(electronCutMap(), "db"+cone_[c])*pu-AnaUtil::cutValue(electronCutMap(), "dbInt"+cone_[c]));
  //do deltabeta from quadratic fit
  double apu=AnaUtil::cutValue(electronCutMap(), "apuIso"+cone_[c]);
  double aN=AnaUtil::cutValue(electronCutMap(), "aNIso"+cone_[c]);
  double bpu=AnaUtil::cutValue(electronCutMap(), "bpuIso"+cone_[c]);
  double bN=AnaUtil::cutValue(electronCutMap(), "bNIso"+cone_[c]);
  double cpu=AnaUtil::cutValue(electronCutMap(), "cpuIso"+cone_[c]);
  double cN=AnaUtil::cutValue(electronCutMap(), "cNIso"+cone_[c]);

  double dbFac= (bN + 2.*cN*ngoodVtx)/(bpu + 2.*cpu*ngoodVtx);
  AnaUtil::fillHist1D("QdbFac_"+cone_[c],dbFac,1);
  double Qdbiso= ch + std::max(0.0, nu - dbFac*pu - aN);
  AnaUtil::fillHist1D("QdbelectronrelIso_"+cone_[c],Qdbiso/elecpt,1);
  //do rho Correction
  //double rhoSlope=AnaUtil::cutValue(electronCutMap(),"brho"+cone_[c]) + 2.*AnaUtil::cutValue(electronCutMap(),"crho"+cone_[c])*fGridRhoFastjetAll;
  double rhoSlope=AnaUtil::cutValue(electronCutMap(),"brho"+cone_[c]);
  AnaUtil::fillHist1D("QrhoFac_"+cone_[c],rhoSlope,1);
  //   double rhoIso=ch + std::max(0.0, nu - rhoSlope*fGridRhoFastjetAll - AnaUtil::cutValue(electronCutMap(),"arho"+cone_[c]));
  double rhoIso=ch + std::max(0.0, nu - rhoSlope*fGridRhoFastjetAll);
  AnaUtil::fillHist1D("RhelectronrelIso_"+cone_[c],rhoIso/elecpt,1);
  //AnaUtil::fillProfile("dbIso"+cone_[c]+"vsnVtx",ngoodVtx,dbiso);
  //AnaUtil::fillHist1D("electronrelIso_"+cone_[c],iso/elecpt,1);
  //vs Rho
  AnaUtil::fillProfile("chIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,ch);
  AnaUtil::fillProfile("NIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,nu);
  AnaUtil::fillProfile("puIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,pu);

  //AnaUtil::fillProfile("rhIso"+cone_[c]+"vsrho",fGridRhoFastjetAll,rhoIso);
  //AnaUtil::fillHist1D("electronRhorelIso_"+cone_[c],rhoIso/elecpt,1);   
  //Vs NVertex
  AnaUtil::fillProfile("chIso"+cone_[c]+"vsnVtx",ngoodVtx,ch);
  AnaUtil::fillProfile("NIso"+cone_[c]+"vsnVtx",ngoodVtx,nu);
  AnaUtil::fillProfile("puIso"+cone_[c]+"vsnVtx",ngoodVtx,pu);
  AnaUtil::fillProfile("rhIso"+cone_[c]+"vsnVtx",ngoodVtx,rhoIso);
  AnaUtil::fillProfile("dbIsoraw"+cone_[c]+"vsnVtx",ngoodVtx,dbiso);
  AnaUtil::fillProfile("QdbIsoraw"+cone_[c]+"vsnVtx",ngoodVtx,Qdbiso);
  AnaUtil::fillProfile("dbIsorel"+cone_[c]+"vsnVtx",ngoodVtx,dbiso/elecpt);
  AnaUtil::fillProfile("QdbIsorel"+cone_[c]+"vsnVtx",ngoodVtx,Qdbiso/elecpt);
  AnaUtil::fillProfile("OdbIsorel"+cone_[c]+"vsnVtx",ngoodVtx,Odbiso/elecpt);
  AnaUtil::fillProfile("rhIsorel"+cone_[c]+"vsnVtx",ngoodVtx,rhoIso/elecpt);

  //vsPt
  AnaUtil::fillProfile("QdbIsoraw"+cone_[c]+"vspt",elecpt,Qdbiso/elecpt);
  ///Added parts for iso efficiency
     if( rhoIso/elecpt < AnaUtil::cutValue(electronCutMap(),"eleRelIso") )
       AnaUtil::fillHist1D("ePassrhIsorel_"+cone_[c],rhoIso/elecpt,1);
     if( dbiso/elecpt < AnaUtil::cutValue(electronCutMap(),"eleRelIso") )
         AnaUtil::fillHist1D("ePassdbIsorel_"+cone_[c],dbiso/elecpt,1);
 }
 else
  std::cout<<"Wrong cone value"<<std::endl;
}
void ElectronIsolation::electronSelector() {
 for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it) {
  const vhtm::Electron& electron = (*it);
  double elePt = electron.pt;
  if (elePt <= AnaUtil::cutValue(electronCutMap(), "pt")) continue;
  double eleta  = std::fabs(electron.eta);
  bool   etaCut = (eleta >= AnaUtil::cutValue(electronCutMap(), "etaLow") &&
    eleta <= AnaUtil::cutValue(electronCutMap(), "etaUp")) ||
   eleta  >= AnaUtil::cutValue(electronCutMap(), "eta");
  if ( elePt <= AnaUtil::cutValue(electronCutMap(), "pt")                         ||
    etaCut                                                                     ||
    electron.dxyPV > AnaUtil::cutValue(electronCutMap(), "dxyPV")              ||
    electron.dzPV > AnaUtil::cutValue(electronCutMap(), "dzPV")                ||
    electron.missingHits > AnaUtil::cutValue(electronCutMap(), "missingHits")  ||
    std::fabs(electron.dB3D/electron.edB3D)  > 4.     
    ) continue;
     
  looseeleVec_.push_back(electron);

  double fSCeta = std::fabs(electron.scEta);
  double BDT=electron.BDT;
  bool isBDT = (elePt <= 10 && (( fSCeta < 0.8 && BDT > 0.47)   ||
     (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > 0.004) ||
     (fSCeta >= 1.479               && BDT > 0.295))) ||
     (elePt >  10 && ((fSCeta < 0.8 && BDT > -0.34)      ||
                    (fSCeta >= 0.8 && fSCeta < 1.479 && BDT > -0.65) ||
                    (fSCeta >= 1.479               && BDT > 0.6)));
  if( isBDT ) tighteleVec_.push_back(electron);

 }
 if( !looseeleVec_.empty() )
  std::sort(looseeleVec_.begin(),looseeleVec_.end(),PtComparator<vhtm::Electron>());
 if( !tighteleVec_.empty() )
  std::sort(tighteleVec_.begin(),tighteleVec_.end(),PtComparator<vhtm::Electron>());
}



// -------------------
// The main event loop
// -------------------
void ElectronIsolation::clearLists() {
 eventElectrons.clear();
 vtxList.clear();
 looseeleVec_.clear();
 tighteleVec_.clear();
 extraeleVec_.clear();
}
void ElectronIsolation::eventLoop() 
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
 //long int nRecoelectrons=0,nGoodRecoelectrons=0;
 std::cout<<"Bunch Crossing>>>>"<<bunchCrossing()<<std::endl;
 for (int ev = 0; ev < nEvents(); ++ev){
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
  op.verbose =(logOption() >> 1 & 0x1);
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
  if(!electronColl()->empty())
     AnaUtil::fillHist1D("nRecoelectron",electronColl()->size(),1);
  int ngoodelectron=0;
  electronSelector();
     std::cout << "Loose Electron size" << looseeleVec_.size();
     findExtraelectron();
     std::cout << "Extra Electron size" << extraeleVec_.size();
     if(ngoodelectron > 0)
         AnaUtil::fillHist1D("nGoodelectron",ngoodelectron,1);
  /*
  for (size_t  it =0; it< looseeleVec_.size(); ++it) {
   const vhtm::Electron& electron = looseeleVec_.at(it);    
   ngoodelectron++;
   AnaUtil::fillHist1D("electron1Pt",electron.pt,1);
   AnaUtil::fillHist1D("electron1Eta",electron.eta,1);

   AnaUtil::fillProfile("elecPtvsnVtx",ngoodVtx,electron.pt);
   if(electron.isolationMap.find("c15") != electron.isolationMap.end() ) {
    fillHistoforCone("c15",electron.pt,electron.isolationMap.at("c15").at(0),electron.isolationMap.at("c15").at(1),electron.isolationMap.at("c15").at(2),ngoodVtx,evt.fGridRhoFastjetAll);   
   }
   if(electron.isolationMap.find("c20") != electron.isolationMap.end() ) {
    fillHistoforCone("c20",electron.pt,electron.isolationMap.at("c20").at(0),electron.isolationMap.at("c20").at(1),
      electron.isolationMap.at("c20").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
   }
   if(electron.isolationMap.find("c25") != electron.isolationMap.end() ) {
    fillHistoforCone("c25",electron.pt,electron.isolationMap.at("c25").at(0),electron.isolationMap.at("c25").at(1),
      electron.isolationMap.at("c25").at(2),ngoodVtx,evt.fGridRhoFastjetAll);

   }
   if(electron.isolationMap.find("c30") != electron.isolationMap.end() ) {
    fillHistoforCone("c30",electron.pt,electron.isolationMap.at("c30").at(0),electron.isolationMap.at("c30").at(1),
      electron.isolationMap.at("c30").at(2),ngoodVtx,evt.fGridRhoFastjetAll);

   }
   if(electron.isolationMap.find("c35") != electron.isolationMap.end() ) {
    fillHistoforCone("c35",electron.pt,electron.isolationMap.at("c35").at(0),electron.isolationMap.at("c35").at(1),
      electron.isolationMap.at("c35").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
   }
   if(electron.isolationMap.find("c40") != electron.isolationMap.end() ) {
    fillHistoforCone("c40",electron.pt,electron.isolationMap.at("c40").at(0),electron.isolationMap.at("c40").at(1),
      electron.isolationMap.at("c40").at(2),ngoodVtx,evt.fGridRhoFastjetAll);

   }
   if(electron.isolationMap.find("c45") != electron.isolationMap.end() ) {
    fillHistoforCone("c45",electron.pt,electron.isolationMap.at("c45").at(0),electron.isolationMap.at("c45").at(1),
      electron.isolationMap.at("c45").at(2),ngoodVtx,evt.fGridRhoFastjetAll);
   }

  }*/
     
 }  
 // Analysis is over
 endJob();
}
//
void ElectronIsolation::findExtraelectron() {
    //std::vector<vhtm::Electron> Zele;
    int ze1idx=-1,ze2idx=-1;
    const float MZnominal=91.1876;
    float MZdiff = 9999.;
    if(looseeleVec_.size() > 2) {
      for (unsigned int  k1 =1; k1 < looseeleVec_.size(); ++k1) {
         TLorentzVector e1P4;
         e1P4.SetPtEtaPhiE(looseeleVec_[k1].pt,
                           looseeleVec_[k1].eta,
                           looseeleVec_[k1].phi,
                           looseeleVec_[k1].energy);
 
        for (unsigned int k2 =0; k2< looseeleVec_.size(); ++k2) {
            if(looseeleVec_[k1].charge + looseeleVec_[k2].charge !=0 ) continue;
            TLorentzVector e2P4;
            e2P4.SetPtEtaPhiE(looseeleVec_[k2].pt,
                              looseeleVec_[k2].eta,
                              looseeleVec_[k2].phi,
                              looseeleVec_[k2].energy);
            double MZtemp=(e1P4+e2P4).M();
            if( MZtemp >= AnaUtil::cutValue(evselCutMap(),"mZ1lowbound") &&
                MZtemp <= AnaUtil::cutValue(evselCutMap(),"mZ1upperbound") &&
               std::fabs(MZtemp - MZnominal) < MZdiff ) {
                MZdiff = std::fabs(MZtemp - MZnominal);
                ze1idx = k1;
                ze2idx = k2;
            }

        }
      }
      if( ze1idx != -1  && ze1idx != -1 ) {
        for( unsigned int  k1 =1; k1 < looseeleVec_.size(); ++k1 ) {
            if( k1 == ze1idx || k1 == ze1idx ) continue;
            extraeleVec_.push_back(looseeleVec_[k1]);
        }
      }
    }
}

void ElectronIsolation::endJob() {

 closeFiles();

 histf()->cd();
 histf()->Write();
 histf()->Close();
 delete histf();
}
