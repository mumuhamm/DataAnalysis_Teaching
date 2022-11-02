void DataMCComparison(string fi1, string fi2, string fi3, float ratio){
  //f1 = full data file, f2= data SB file, f3= signal MC
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000);

  TFile* f1= new TFile(fi1.c_str());
  TFile* f2= new TFile(fi2.c_str());
  TFile* f3= new TFile(fi3.c_str());

  string name=string(f1->GetName());
  name="DataMC_Comparison_" + name;

  TFile* fout= new TFile(name.c_str(),"RECREATE");

  TIter nextkey(f1->GetListOfKeys());
  while (TKey* key = (TKey*)nextkey()){

    TObject* obj=key->ReadObj();
    if (!obj->InheritsFrom("TH1F")) continue;

    TH1F* hdata=(TH1F*)f1->Get(obj->GetName());
    TH1F* hsb=(TH1F*)f2->Get(obj->GetName());
    TH1F* hmc=(TH1F*)f3->Get(obj->GetName());
 
    hdata->SetLineColor(1);
    hsb->SetLineColor(2);
    hmc->SetLineColor(4);

    TH1F* hsb_clone=(TH1F*)hsb->Clone();
    hsb_clone->Scale(ratio);

    hdata->Add(hsb_clone,-1);

    double scale=hdata->Integral();
    hsb->Scale(scale/hsb->Integral());
    hmc->Scale(scale/hmc->Integral());

    double max=hdata->GetMaximum();
    if (hsb->GetMaximum()>max) max=hsb->GetMaximum();
    if (hmc->GetMaximum()>max) max=hmc->GetMaximum();

    fout->cd();
    string hname;

    TCanvas c;
    hdata->Draw();
    hsb->Draw("samehisto");
    hmc->Draw("samehisto");

    hname=string(hdata->GetName())+"Data_SB_Subtracted";
    hdata->SetName(hname.c_str());
    hdata->Write();
    hname="";

    hname=string(hsb->GetName())+"Data_SB";
    hsb->SetName(hname.c_str());
    hsb->Write();
    hname="";

    hname=string(hmc->GetName())+"Signal_MC";
    hmc->SetName(hname.c_str());
    hmc->Write();
    hname="";

    hdata->GetYaxis()->SetRangeUser(0, max*1.2);

    TLegend leg(0.65,0.65,0.85,0.85);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(hdata, "Data SB sub");
    leg.AddEntry(hsb, "Side-Bands");
    leg.AddEntry(hmc, "MC Signal");
    leg.Draw();

    string out="DataMc/"+string(obj->GetName())+".png";
    c.SaveAs(out.c_str());
  }

  fout->Close();

}
