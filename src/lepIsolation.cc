#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStopwatch.h"

#include "LeptonIsolation.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " jobFile " << endl;
    exit(0);
  }     
  string jobFile(argv[1]);

  // Create  analysis object 
  LeptonIsolation myeff;

  // Read job input
  int nFiles = 0;
  bool succeed = myeff.readJob(jobFile, nFiles);
  if (!succeed) exit(1);
  if (myeff.getEntries() <= 0) {
    cerr << "No events present in the input chain, exiting ...!" << endl;
    exit(2);
  }
  
  gROOT->SetBatch(kTRUE);
  
  // Now go
  TStopwatch timer;
  cout << "Start event loop now with " << nFiles << " file(s)" << endl;
  timer.Start();
  myeff.eventLoop();
  timer.Stop();
  cout << "Realtime/CpuTime = " << timer.RealTime() << "/" << timer.CpuTime() << endl;
  
  return 0;
}
