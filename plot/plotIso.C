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
#include<TProfile.h>
#include<TPaveStats.h>
#include<TCanvas.h>
#include<TColor.h>
//plot different IsovsRho on same canvas
void doNisovRhoplots(TFile* fin,TFile* fout){
 TString cones[]={"15","20","25","30","35","40","45"};
 for(int i=0; i<7; i++ ) {
   TString chname="chIso0"+cones[i]+"vsrho";
   TString Nname="NIso0"+cones[i]+"vsrho";
   TString puname="puIso0"+cones[i]+"vsrho";
   TProfile *chP = dynamic_cast<TProfile*>(fin->Get(chname));
   chP->SetLineColor(kBlue); 
   chP->SetLineWidth(2);
   TProfile *NP = dynamic_cast<TProfile*>(fin->Get(Nname));
   NP->SetLineColor(kGreen);
   NP->SetLineWidth(2);
   TProfile *puP = dynamic_cast<TProfile*>(fin->Get(puname));
   puP->SetLineColor(kRed);
   puP->SetLineWidth(2);
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",323,75,1039,567);
   puP->Draw();
   NP->Draw("sames");
   chP->Draw("sames");
   fout->cd();
   Canvas_1->Write();
   fin->cd();
   Canvas_1->Delete();
 }
}
//plot Ch/N/Pu Relative Iso for different cones
void doIsoforconeplots(TString p,TFile* fin,TFile* fout){
 TString cones[]={"15","20","25","30","35","40","45"};
 for(int i=0; i<7; i++ ) {
   TString chname=p+"ChIso_0"+cones[i];
   TString Nname=p+"NIso_0"+cones[i];
   TString puname=p+"PuIso_0"+cones[i];
   TH1F *chP = dynamic_cast<TH1F*>(fin->Get(chname));
   chP->SetLineColor(kBlue); 
   chP->SetLineWidth(2);
   TH1F *NP = dynamic_cast<TH1F*>(fin->Get(Nname));
   NP->SetLineColor(TColor::GetColor("#009900"));
   NP->SetLineWidth(2);
   TH1F *puP = dynamic_cast<TH1F*>(fin->Get(puname));
   puP->SetLineColor(kRed);
   puP->SetLineWidth(2);
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",323,75,1039,567);
   chP->Draw();
   NP->Draw("sames");
   puP->Draw("sames");
   fout->cd();
   Canvas_1->Write();
   fin->cd();
   Canvas_1->Delete();
 }
}
//plot Same Relative Iso for different cones
void doIsoforconeplots(TString p,TFile* fin,TFile* fout){
 TString cones[]={"15","20","25","30","35","40","45"};
 for(int i=0; i<7; i++ ) {
   TString hname=p+"ChIso_0"+cones[i];
   TH1F *hP = dynamic_cast<TH1F*>(fin->Get(chname));
   hP->SetLineColor(i); 
   hP->SetLineWidth(2);
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",323,75,1039,567);
   chP->Draw();
   NP->Draw("sames");
   puP->Draw("sames");
   fout->cd();
   Canvas_1->Write();
   fin->cd();
   Canvas_1->Delete();
 }
}
void mainPlot() {
 TFile *fin = TFile::Open("test_HZZ4lmuIso.root");
 TFile *fout = TFile::Open("test_HZZ4lmuIsovsRhoplot.root","RECREATE");
 doNisovRhoplots(fin,fout);
 doIsoforconeplots("muon",fin,fout);
 fin->Close();
 fout->Close();
}
