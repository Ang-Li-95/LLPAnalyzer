using namespace std;

vector<vector<TH2F*>> scalePlots2D(vector<TFile*> fsig, vector<double> sigLumiScale, vector<vector<TFile*>> fbkg, vector<vector<double>> lumiScale, TString name)
{
  vector<TH2F*> hist_sig;
  for (size_t i=0; i<fsig.size(); ++i){
    TH2F* h = (TH2F*)fsig[i]->Get(name);
    h->Scale(sigLumiScale[i]);
    hist_sig.push_back(h);
  }
  vector<TH2F*> hist_bkg;
  for (size_t i=0; i<fbkg.size(); ++i){
    for (size_t j=0; j<fbkg[i].size(); ++j){
      TH2F* h = (TH2F*)fbkg[i][j]->Get(name);
      h->Scale(lumiScale[i][j]);
      hist_bkg.push_back(h);
    }
  }
  vector<vector<TH2F*>> hist {hist_sig, hist_bkg};
  return hist;
}

void plotHists2D(vector<vector<TH2F*>> hist, TString name){
  TCanvas *c = new TCanvas("c_"+name, "c_"+name, 600, 600);
  c->cd();
  TH2F* hsig;
  for (size_t i=0; i<hist[0].size(); ++i){
    if (i==0){
      hsig = (TH2F*)hist[0][i]->Clone();
    }
    else{
      hsig->Add(hist[0][i]);
    }
  }
  TH2F* hbkg;
  for (size_t i=0; i<hist[1].size(); ++i){
    if(i==0){
      hbkg = hist[1][i];
    }
    else{
      hbkg->Add(hist[1][i]);
    }
  }
  hbkg->Draw("colz");
}

vector<vector<TH1F*>> scalePlots(vector<TFile*> fsig, vector<double> sigLumiScale, vector<vector<TFile*>> fbkg, vector<vector<double>> lumiScale, TString name)
{
  vector<TH1F*> hist_sig;
  for (size_t i=0; i<fsig.size(); ++i){
    TH1F* h = (TH1F*)fsig[i]->Get(name);
    h->Scale(sigLumiScale[i]);
    hist_sig.push_back(h);
  }
  vector<TH1F*> hist_bkg;
  for (size_t i=0; i<fbkg.size(); ++i){
    for (size_t j=0; j<fbkg[i].size(); ++j){
      TH1F* h = (TH1F*)fbkg[i][j]->Get(name);
      h->Scale(lumiScale[i][j]);
      hist_bkg.push_back(h);
    }
  }
  vector<vector<TH1F*>> hist {hist_sig, hist_bkg};
  return hist;
}

void plotHists(vector<vector<TH1F*>> hist, TString name){
  TCanvas *c = new TCanvas("c_"+name, "c_"+name, 600, 600);
  c->cd();
  TH1F* hsig;
  for (size_t i=0; i<hist[0].size(); ++i){
    if (i==0){
      hsig = (TH1F*)hist[0][i]->Clone();
    }
    else{
      hsig->Add(hist[0][i]);
    }
  }
  hsig->SetLineColor(kGreen);
  hsig->SetLineWidth(3);
  TH1F* hbkg;
  for (size_t i=0; i<hist[1].size(); ++i){
    if(i==0){
      hbkg = hist[1][i];
    }
    else{
      hbkg->Add(hist[1][i]);
    }
  }
  hbkg->SetLineColor(38);
  hbkg->SetFillColor(38);
  if(name == "vtx_tkSize")
    hbkg->SetMinimum(0.05);
  else
    hbkg->SetMinimum(0.5);
  hbkg->SetMaximum(1e+06);
  hbkg->Draw("hist");
  hsig->Draw("same hist");
  c->SetLogy();
}

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

void make2DPlot(vector<TFile*> f, vector<double> lumiScale,TString name)
{
  TCanvas *c = new TCanvas("c_"+name, "c_"+name, 600, 600);
  c->cd();
  TH2F *h;
  for(int i=0; i<f.size(); ++i){
    TH2F *hi = (TH2F*)f[i]->Get(name);
    if(i==0){
      h = (TH2F*)hi->Clone();
      h->Scale(lumiScale[i]);
    }
    else{
      h->Add(hi, lumiScale[i]);
    }
  }
  h->Draw("colz");
  //c->SetLogy();
}
void make1DStackPlot(vector<TFile*> f, vector<double> lumiScale, TFile* fsig, double sigLumiScale, TString name)
{
  const TString l[8] = {"QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf", "", "", "", ""};
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

void WeightedPlotsYears()
{
  gStyle->SetOptStat(0); 
  // dir[0] is 2017 directory
  // dir[1] is 2018 directory
  const vector<TString> dir {"/uscms/home/ali/nobackup/LLP/DVAnalyzer/CMSSW_9_4_15/src/LLPAnalyzer/DVAnalyzer/test/", "/uscms/home/ali/nobackup/LLP/DVAnalyzer/CMSSW_10_2_19/src/LLPAnalyzer/DVAnalyzer/test/"};
  const vector<double> lumi {41.5, 59.7};
  const TString signame = "ggToNN_800M_1mm.root";
  const vector<TString> Bkg {"QCD_HT700to1000.root",
                               "QCD_HT1000to1500.root",
                               "QCD_HT1500to2000.root",
                               "QCD_HT2000toInf.root",
                               "TTJets600To800.root",
                               "TTJets800To1200.root",
                               "TTJets1200To2500.root",
                               "TTJets2500ToInf.root"
  };
  //const vector<int> Events_bkg_2017 {47724800, 16595628, 11634434, 5941306, 81507662, 40187347, 13214871, 5155687};
  //const vector<int> Events_bkg_2018 {48158738, 15466225, 10955087, 5475677, 14149394, 10372802, 2779427, 1451104};
  const vector<int> Events_bkg_2017 {48042655, 16882838, 11634434, 5941306, 81565576, 40248127, 13214871, 5155687};
  const vector<int> Events_bkg_2018 {43523821, 15174716, 11082955, 5557453, 14363689, 10462756, 2897601, 1451104};
  const vector<vector<int>> Events_bkg {Events_bkg_2017, Events_bkg_2018};
  const vector<double> sigLumiScale {lumi[0]/100000,lumi[1]/100000};
  const vector<double> crossSection {6.4e+06, 1.1e+06, 9.9e+04, 2.0e+04, 1.8e+03, 0.75e+03, 0.13e+03, 1.41};

  const int Nbkg = Bkg.size();
  const int Nyear = dir.size();
  vector<TFile*> fsig(Nyear);
  vector<vector<TFile*>> fBkg;
  vector<vector<double>> scaleFactors;
  for(size_t iy=0; iy<Nyear; ++iy){
    fsig[iy] = new TFile(dir[iy] + signame);
    vector<TFile*> f(Nbkg);
    vector<double> scaleFactor(Nbkg);
    for (int i =0;i<Nbkg;++i){
      f[i] = new TFile(Bkg[i]);
      scaleFactor[i] = lumiScale(Events_bkg[iy][i], crossSection[i], lumi[iy]);
    }
    fBkg.push_back(f);
    scaleFactors.push_back(scaleFactor);

  }
  vector<vector<TH1F*>> tkSize = scalePlots(fsig, sigLumiScale, fBkg, scaleFactors, "vtx_tkSize");
  plotHists(tkSize, "vtx_tkSize");

  vector<vector<TH1F*>> dBV = scalePlots(fsig, sigLumiScale, fBkg, scaleFactors, "vtx_dBV");
  plotHists(dBV, "vtx_dBV");

  vector<vector<TH1F*>> sigma_dBV = scalePlots(fsig, sigLumiScale, fBkg, scaleFactors, "vtx_sigma_dBV");
  plotHists(sigma_dBV, "vtx_sigma_dBV");

  vector<vector<TH2F*>> vtx_xy = scalePlots2D(fsig, sigLumiScale, fBkg, scaleFactors, "vtx_xy");
  plotHists2D(vtx_xy, "vtx_xy");
  //make1DPlot(f, scaleFactor, "vtx_tkSize");
  //make1DPlot(f, scaleFactor, "vtx_dBV");
  //make1DPlot(f, scaleFactor, "vtx_sigma_dBV");

  //make1DStackPlot(f, scaleFactor, fsig, sigLumiScale, "vtx_tkSize");
  //make1DStackPlot(f, scaleFactor, fsig, sigLumiScale, "vtx_dBV");
  //make1DStackPlot(f, scaleFactor, fsig, sigLumiScale, "vtx_sigma_dBV");
  //AddPlots2D(f, scaleFactor, fsig, sigLumiScale, "vtx_xy");

}
