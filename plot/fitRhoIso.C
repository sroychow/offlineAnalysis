#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<fstream>
#include<sstream>
#include<TROOT.h>
#include<TH1.h>
#include<TF1.h>
#include<TFile.h>
#include<TMath.h>
#include<TKey.h>
#include<TClass.h>
//#include<TIter.h>
#include<TProfile.h>
#include<TPaveStats.h>

//cone-arrray
TString cone[]={"015","020","025","030","035","040","045"};

//function to bind a stat box to a histogram
template<typename T>
void bindstat(T* h,double x1ndc,double y1ndc, double x2ndc,double y2ndc)
{
 TPaveStats *ptstats = new TPaveStats(x1ndc,y1ndc,x2ndc,y2ndc,"brNDC");
 h->GetListOfFunctions()->Add(ptstats);
 //ptstats->SetOptStat(1111);
 ptstats->SetOptFit(1111);
 ptstats->SetParent(h);
}


// define the linear fit function 
   Double_t fitf(Double_t *x,Double_t *par) {
      return par[0] + par[1]*x[0];
   }
// define the quadratic fit function 
   Double_t fitQf(Double_t *x,Double_t *par) {
      return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
   }


void dofitRhoIso(TFile* fin,TFile* fout, int rLow,int rUp,std::ofstream& ff) {

      
      // open a file and get a histogram
      //TFile *fin = TFile::Open("test_HZZ4lmuIso.root");
      //TFile *fout = TFile::Open("test_HZZ4lmuRhoIsofit.root","RECREATE");
      TIter next(fin->GetListOfKeys());
      TKey *key;
      
      while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if( cl->InheritsFrom("TProfile") ) {
          if(std::string(key->GetName()).find("rho")==std::string::npos ) continue;
          std::cout<<key->GetName()<<std::endl;
          TProfile *hpx = dynamic_cast<TProfile*>(fin->Get(key->GetName()));
          TF1 *func = new TF1("fit",fitf,0,100,2);
          //func->SetParameters(0.0768942,0.0215659);
          // give the parameters meaningful names
          func->SetParNames ("Constant","Slope");
          // call TH1::Fit with the name of the TF1 object
          hpx->Fit("fit","", "",rLow,rUp);
          bindstat<TProfile>(hpx,0.7797101,0.4298401,0.9797101,0.9360568);         
          fout->cd();
          hpx->Write();
          fin->cd();
          ff << "Slope=" << func->GetParameter("Slope") << std::endl;
          ff << "Constant=" << func->GetParameter("Constant") << std::endl;
        }
      }

   }
//Do quadratic fit
void doQfitRhoIso(TFile* fin,TFile* fout, int rLow,int rUp,std::ofstream& ff) {

      
      // open a file and get a histogram
      //TFile *fin = TFile::Open("test_HZZ4lmuIso.root");
      //TFile *fout = TFile::Open("test_HZZ4lmuRhoIsofit.root","RECREATE");
      TIter next(fin->GetListOfKeys());
      TKey *key;
      
      while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if( cl->InheritsFrom("TProfile") ) {
          if(std::string(key->GetName()).find("rho")==std::string::npos ) continue;
          ff<<key->GetName()<<std::endl;
          TProfile *hpx = dynamic_cast<TProfile*>(fin->Get(key->GetName()));
          TF1 *func = new TF1("fit",fitQf,0,100,3);
          //func->SetParameters(0.0768942,0.0215659);
          // give the parameters meaningful names
          func->SetParNames ("a","b","c");
          // call TH1::Fit with the name of the TF1 object
          hpx->Fit("fit","", "",rLow,rUp);
          bindstat<TProfile>(hpx,0.7797101,0.4298401,0.9797101,0.9360568);         
          fout->cd();
          hpx->Write();
          fin->cd();
          ff<< "a=" << func->GetParameter("a") <<"+/-" << func->GetParError(func->GetParNumber("a")) 
            << " b=" << func->GetParameter("b") <<"+/-" << func->GetParError(func->GetParNumber("b")) 
            << " c=" << func->GetParameter("c") <<"+/-" << func->GetParError(func->GetParNumber("c")) 
            << std::endl;
        }
      }

   }
//do linear fit
void doIsonVtxLfit(TFile* fin,TFile* fout, int rLow,int rUp,TString isotype,std::ofstream& ff) {
 ff<< isotype + "IsovsNVtx Linear Fit results"<<std::endl;
 for(int c=0;c<7;c++) {
  TProfile *hpx = dynamic_cast<TProfile*>(fin->Get(isotype+cone[c]+"vsnVtx"));
  TF1 *func = new TF1("fit",fitf,0,100,2);
  func->SetParNames ("c","m");
  hpx->Fit("fit","", "",rLow,rUp);
  ff<<cone[c]<<"\t";
  bindstat<TProfile>(hpx,0.7797101,0.4298401,0.9797101,0.9360568); 
  fout->cd();
  hpx->Write();
  fin->cd();
  ff<< "c"+isotype+cone[c]+"=" << func->GetParameter("c") <<"+/-" << func->GetParError(func->GetParNumber("c")) 
    << " m"+isotype+cone[c]+"=" << func->GetParameter("m") <<"+/-"<< func->GetParError(func->GetParNumber("m"))
    << std::endl;
 }
}
//do quadratic fit
void doIsonVtxQfit(TFile* fin,TFile* fout, int rLow,int rUp,TString isotype,std::ofstream& ff) {
 //ofstream ff("slopes.txt");
 ff<< isotype + "IsovsNVtx Quadratic Fit results"<<std::endl;
 for(int c=0;c<7;c++) {
  TProfile *hpx = dynamic_cast<TProfile*>(fin->Get(isotype+cone[c]+"vsnVtx"));
  TF1 *func = new TF1("fit",fitQf,0,100,3);
  func->SetParNames ("a","b","c");
  hpx->Fit("fit","", "",rLow,rUp);
  ff<<cone[c]<<"\t";
  bindstat<TProfile>(hpx,0.7797101,0.4298401,0.9797101,0.9360568); 
  fout->cd();
  hpx->Write();
  fin->cd();
  ff<< "a"+isotype+cone[c]+"=" << func->GetParameter("a") <<"+/-" << func->GetParError(func->GetParNumber("a")) 
    << " b"+isotype+cone[c]+"=" << func->GetParameter("b") <<"+/-"<< func->GetParError(func->GetParNumber("b"))
    << " c"+isotype+cone[c]+"=" << func->GetParameter("c") <<"+/-"<< func->GetParError(func->GetParNumber("c"))
    << std::endl;
 }
 //ff.close();
}

int main(int argc,char* argv[]) {
 TString inFile(argv[1]);
 TString outFile(argv[2]);
 int rlow = atoi(argv[3]);
 int rup = atoi(argv[4]);
 if( argc < 4 ) {
   std::cout<<"Usage: fitiso [inputFile] [outputFile] [LowerBound] [UpperBound]" << std::endl;
   return 0;
 } 
 std::cout<<"Infile>>"<<inFile
          <<"\nOutFile>>"<<outFile
          <<"\nLowerRange>>"<<rlow
          <<"\nUpperRange>>"<<rup<<std::endl;
 TFile *fin = TFile::Open(inFile);
 TFile *fout = TFile::Open(outFile,"RECREATE");
 ofstream ff("slopes.txt");
 //doQfitRhoIso(fin,fout,rlow,rup,ff);
 //doIsonVtxfit(fin,fout,8,28,"puIso",ff);
 //doIsonVtxfit(fin,fout,8,28,"NIso",ff);
 doIsonVtxQfit(fin,fout,rlow,rup,"puIso",ff);
 doIsonVtxQfit(fin,fout,rlow,rup,"NIso",ff);

 fin->Close();
 fout->Close();
 ff.close();
 return 0;
}
