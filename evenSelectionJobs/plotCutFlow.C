void plotCutFlow(){
 ofstream fout("cutFlow.txt");
 TFile* fin=TFile::Open("test_evetSelectn_PHYS14DR_PU20bx25_bkg_iso03.root");
/*
 TH1F* htemp = dynamic_cast<TH1F*>(fin->Get("evntCutFlow_4mu"));
 std::fout<<"4mu cutflow"<<std::endl;
 for( int i = 1;i<15;i++)
   std::fout << "Bin:"<<i<<"="<<htemp->GetBinContent(i)<<"("<<(htemp->GetBinContent(i)*100)/htemp->GetBinContent(2)<<")"<<std::endl;
 //htemp->Draw();
 TH1F* htemp1 = dynamic_cast<TH1F*>(fin->Get("evntCutFlow_4e"));
 std::fout<<"4e cutflow"<<std::endl;
 for( int i = 1;i<15;i++)
   std::fout << "Bin:"<<i<<"="<<htemp1->GetBinContent(i)<<"("<<(htemp1->GetBinContent(i)*100)/htemp1->GetBinContent(2)<<")"<<std::endl;

 TH1F* htemp2 = dynamic_cast<TH1F*>(fin->Get("evntCutFlow_2e2mu"));
 std::fout<<"2e2mu cutflow"<<std::endl;
 for( int i = 1;i<15;i++)
   std::fout << "Bin:"<<i<<"="<<htemp2->GetBinContent(i)<<"("<<(htemp2->GetBinContent(i)*100)/htemp2->GetBinContent(2)<<")"<<std::endl;
 fout.close();
*/
 TH1F* htemp3 = dynamic_cast<TH1F*>(fin->Get("evntCutFlow_bkg"));
 std::fout<<"bkg cutflow"<<std::endl;
 for( int i = 1;i<15;i++)
   std::fout << "Bin:"<<i<<"="<<htemp3->GetBinContent(i)<<"("<<(htemp3->GetBinContent(i)*100)/htemp3->GetBinContent(1)<<")"<<std::endl;
 fout.close();
}
