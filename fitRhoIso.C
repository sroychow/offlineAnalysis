// define a function 
   Double_t fitf(Double_t *x,Double_t *par) {
      return fitval = par[0] + par[1]*x[0];
   }

void fitRhoIso() {

      
      // open a file and get a histogram
      TFile *fin = TFile::Open("test_HZZ4lmuIso.root");
      TFile *fout = TFile::Open("test_HZZ4lmuRhoIsofit.root","RECREATE");
      TIter next(fin->GetListOfKeys());
      TKey *key;
      
      while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if( cl->InheritsFrom("TProfile") ) {
          std::cout<<key->GetName()<<std::endl;
          TProfile *hpx = dynamic_cast<TProfile*>(fin->Get(key->GetName()));
          TF1 *func = new TF1("fit",fitf,0,100,2);
          //func->SetParameters(0.0768942,0.0215659);
          // give the parameters meaningful names
          func->SetParNames ("Constant","Slope");
          // call TH1::Fit with the name of the TF1 object
          hpx->Fit("fit","", "",5,40);
          fout->cd();
          hpx->Write();
          fin->cd();
          std::cout << "Slope=" << func->GetParameter("Slope") << std::endl;
          std::cout << "Constant=" << func->GetParameter("Constant") << std::endl;
        }
      }
      fin->Close();
      fout->Close();
   }
