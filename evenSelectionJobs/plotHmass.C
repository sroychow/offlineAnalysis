void plotHmass(){
 THStack *hs= new THStack("hm4l","4l mass");
 TFile* fin1=TFile::Open("test_evetSelectn_PHYS14DR_PU20bx25_sig_iso03.root");
 TFile* fin2=TFile::Open("test_evetSelectn_PHYS14DR_PU20bx25_bkg_iso03.root");
  
 TH1F* hbkg = dynamic_cast<TH1F*>(fin2->Get("newmassH")); 
 hbkg->Scale(13.78/1958600);
 hbkg->SetFillColor(kBlue);
 hs->Add(hbkg);


 TH1F* hsig = dynamic_cast<TH1F*>(fin1->Get("newmassH"));
 hsig->Scale(1./204864);
 hsig->SetFillColor(kRed);
 hs->Add(hsig);


 

 //hs->Draw();

 //fin1->Close();
 //fin2->Close();
}
