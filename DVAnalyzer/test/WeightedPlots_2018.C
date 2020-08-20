using namespace std;

double lumiScale(int N_evt, double crossSection, double L_int)
{
  return L_int*crossSection/N_evt;
}

void make1DPlot(vector<TFile*> f, vector<double> lumiScale,TString name)
{
  TCanvas *c = new TCanvas("c_"+name, "c_"+name, 600, 600);
  c->cd();
  TH1F *h;
  for(int i=0; i<f.size(); ++i){
    TH1F *hi = (TH1F*)f[i]->Get(name);
    if(i==0){
      h = (TH1F*)hi->Clone();
      h->Scale(lumiScale[i]);
    }
    else{
      h->Add(hi, lumiScale[i]);
    }
  }
  h->Draw("hist");
  c->SetLogy();
}

void make1DStackPlot(vector<TFile*> f, vector<double> lumiScale, TFile* fsig, double sigLumiScale, TString name)
{
  const TString l[4] = {"QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};
  TCanvas *c = new TCanvas("cs_"+name, "cs_"+name, 600, 600);
  c->cd();
  THStack *hs = new THStack("hs"+name,"hs"+name);
  auto legend = new TLegend(0.65,0.7,0.9,0.9);
  for(int i=0; i<f.size(); ++i){
    TH1F *hi = (TH1F*)f[i]->Get(name);
    hi->Scale(lumiScale[i]);
    hi->SetFillColor(i+6);
    hi->SetLineColor(i+6);
    hs->Add(hi,"hist");
    legend->AddEntry(hi,l[i]);
  }
  if(name == "vtx_tkSize")
    hs->SetMinimum(0.1);
  else
    hs->SetMinimum(1);
  hs->SetMaximum(1e+06);
  hs->Draw();
  TH1F * h = (TH1F*)fsig->Get(name);
  h->Scale(sigLumiScale);
  h->SetLineColor(kRed);
  h->SetLineWidth(3);
  h->Draw("same hist");
  legend->AddEntry(h,"signal");
  c->SetLogy();
  legend->Draw();
}

void AddPlots2D(vector<TFile*> fBkg, vector<double> lumiScale, TFile* fsig, double sigLumiScale, TString name)
{
  TCanvas *c = new TCanvas("c_"+name, "c_"+name, 800,800);
  c->cd();
  TH2F *ts = new TH2F("ts","vertex y vs. x", 80, -4, 4, 80, -4, 4);
  TH2F *h = (TH2F*)fsig->Get(name);
  ts->Add(h, sigLumiScale);
  for(int i=0; i<fBkg.size(); ++i){
    TH2F *hi = (TH2F*)fBkg[i]->Get(name);
    ts->Add(hi, lumiScale[i]);
  }
  ts->Draw("colz");
}

void WeightedPlots()
{
  gStyle->SetOptStat(0); 
  const double lumi = 59.7;
  const TString signame = "ggToNN_800M_1mm.root";
  const vector<TString> Bkg {"QCD_HT700to1000.root",
                               "QCD_HT1000to1500.root",
                               "QCD_HT1500to2000.root",
                               "QCD_HT2000toInf.root",
                               "TTJets_HT600To800.root",
                               "TTJets_HT800To1200.root",
                               "TTJets_HT1200To2500.root",
                               "TTJets_HT2500ToInf.root"
  };
  const vector<int> Events_bkg {48158738, 15466225, 10955087, 5475677, 14149394, 10372802, 2779427, 1451104};
  const vector<double> crossSection {6.4e+06, 1.1e+06, 9.9e+04, 2.0e+04, 1.8e+03, 0.75e+03, 0.13e+03, 1.41};
  //const TString QCD[4] = {"QCD_HT700to1000.root",
  //                  "QCD_HT1000to1500.root",
  //                  "QCD_HT1500to2000.root",
  //                  "QCD_HT2000toInf.root"
  //};
  //const TString TTJet[4] = {"TTJets_HT600To800.root",
  //                          "TTJets_HT800To1200.root",
  //                          "TTJets_HT1200To2500.root",
  //                          "TTJets_HT2500ToInf.root"
  //};
  //const int Events_QCD[4] = {48158738, 15466225, 10955087, 5475677};
  //const double crossSection[4] = {6.4e+06, 1.1e+06, 9.9e+04, 2.0e+04};

  //const int Events_TTJets[4] = {14149394, 10372802, 2779427, 1451104};
  //const double crossSection_TTJets[4] = {1.8e+03, 0.75e+03, 0.13e+03, 1.41};
  const double sigLumiScale = lumi/100000;

  TFile* fsig = new TFile(signame);
  const int Nbkg = Bkg.size();
  vector<TFile*> f(Nbkg);
  vector<double> scaleFactor(Nbkg);
  for (int i =0;i<Nbkg;++i){
    f[i] = new TFile(Bkg[i]);
    scaleFactor[i] = lumiScale(Events_bkg[i], crossSection[i], lumi);
  }

  //make1DPlot(f, scaleFactor, "vtx_tkSize");
  //make1DPlot(f, scaleFactor, "vtx_dBV");
  //make1DPlot(f, scaleFactor, "vtx_sigma_dBV");

  make1DStackPlot(f, scaleFactor, fsig, sigLumiScale, "vtx_tkSize");
  make1DStackPlot(f, scaleFactor, fsig, sigLumiScale, "vtx_dBV");
  make1DStackPlot(f, scaleFactor, fsig, sigLumiScale, "vtx_sigma_dBV");
  AddPlots2D(f, scaleFactor, fsig, sigLumiScale, "vtx_xy");

}
