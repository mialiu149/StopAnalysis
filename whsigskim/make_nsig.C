#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

void make_nsig(TString dir) {

  TChain* t = new TChain("Events");
  t->Add(Form("%s/*.root",dir.Data()));

  TFile* fout = new TFile("nsig_TChiWH.root","RECREATE");

  // default: 25 GeV binning, m1 from 0-700, m2 from 0-500
  int x_nbins = 29;
  float x_min = -12.5;
  float x_max = 712.5;
  int y_nbins = 21;
  float y_min = -12.5;
  float y_max = 512.5;

  TH2D* h_nsig = new TH2D("h_nsig",";mass1 [GeV];mass2 [GeV]", x_nbins, x_min, x_max, y_nbins, y_min, y_max);

  t->Draw("sparm_values[1]:sparm_values[0]>>h_nsig");

  fout->Write();
  fout->Close();

}
