#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;
using namespace TMath;

TCanvas* GetSeparation(string);
int main(int argc, char** argv)
{
  string outName  = "result/separation.root";
  if(argc>=2) outName = argv[1];
  TFile* outFile = new TFile(outName.c_str(),"recreate");

  cout<<"Input  file: "<<endl;
  string input  = "23.00_45.00_2.00";
  if(argc>=3)
  {
    for(int i=1;i<argc-1;i++)
    {
      input  = argv[i+1];
      TCanvas* Canvas = GetSeparation(input);
      outFile->cd();
      Canvas->Write();
    }
  }
  outFile->Close();
  cout<<"Output file:\n             "<<outName<<endl;

  return 0;
}

TCanvas* GetSeparation(string input)
{
  double track[3];
  stringstream ss(input);
  string tmp;
  for(int i=0;i<3;i++)
  {
    getline(ss,tmp,'_');
    track[i]=atof(tmp.c_str());
  }

  string pionName = "result/"+input+"_pion.root";
  string kaonName = "result/"+input+"_kaon.root";
  cout<<"             "<<pionName<<endl;
  cout<<"             "<<kaonName<<endl;

  TCanvas *Canvas = new TCanvas(input.c_str(), input.c_str());
   Canvas->SetLeftMargin(0.18);
   Canvas->SetRightMargin(0.16);
   Canvas->SetBottomMargin(0.16);

  TFile* pionFile = new TFile(pionName.c_str());
  TH1F* pionHist = new TH1F(*((TH1F*) pionFile->Get("ProbL")));
  pionHist->Scale(1./pionHist->Integral("width"));
  pionHist->SetTitle("#pi/K separation");
  pionHist->SetName("pion");
  pionHist->SetStats(0);
  pionHist->SetFillColor(6);
  pionHist->SetFillStyle(3004);
  pionHist->SetLineColor(6);
  pionHist->SetLineWidth(2);
  pionHist->GetXaxis()->SetTitle("Log L_{#pi}- Log L_{K}");
  pionHist->GetXaxis()->CenterTitle(true);
  pionHist->GetXaxis()->SetLabelSize(0.05);
  pionHist->GetXaxis()->SetTitleSize(0.05);
  pionHist->GetYaxis()->SetTitle("probability density");
  pionHist->GetYaxis()->CenterTitle(true);
  pionHist->GetYaxis()->SetLabelSize(0.05);
  pionHist->GetYaxis()->SetTitleSize(0.05);

  TFile* kaonFile = new TFile(kaonName.c_str());
  TH1F* kaonHist = new TH1F(*((TH1F*) kaonFile->Get("ProbL")));
  kaonHist->Scale(1./kaonHist->Integral("width"));
  kaonHist->SetTitle("#pi/K separation");
  kaonHist->SetName("kaon");
  kaonHist->SetStats(0);
  kaonHist->SetFillColor(4);
  kaonHist->SetFillStyle(3004);
  kaonHist->SetLineColor(4);
  kaonHist->SetLineWidth(2);
  kaonHist->GetXaxis()->SetTitle("Log L_{#pi}- Log L_{K}");
  kaonHist->GetXaxis()->CenterTitle(true);
  kaonHist->GetXaxis()->SetLabelSize(0.05);
  kaonHist->GetXaxis()->SetTitleSize(0.05);
  kaonHist->GetYaxis()->SetTitle("probability density");
  kaonHist->GetYaxis()->CenterTitle(true);
  kaonHist->GetYaxis()->SetLabelSize(0.05);
  kaonHist->GetYaxis()->SetTitleSize(0.05);

  double counts = 0;
  for(int i=1;i<=pionHist->GetNbinsX();i++)
  {
    if(pionHist->GetBinContent(i)<kaonHist->GetBinContent(i))
      counts += pionHist->GetBinContent(i)*pionHist->GetBinWidth(i);
    else
      counts += kaonHist->GetBinContent(i)*kaonHist->GetBinWidth(i);
  }

  double separation = Sqrt(2)*ErfInverse(1.-counts)*2.; 
  cout<<"             π/K separation = "<<setprecision(3)<<separation<<"σ"<<endl;

  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.9);
  leg1->SetTextColor(2);
  leg1->SetLineWidth(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(pionHist,"#pi sample","lpf");
  leg1->AddEntry(kaonHist,"K sample","lpf");
  TLegend* leg2 = new TLegend(0.2,0.6,0.5,0.9);
  leg2->SetTextColor(2);
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry("NULL",Form("p = %.2fGeV",track[2]),"h");
  leg2->AddEntry("NULL",Form("#theta = %.2f^{o}, #varphi = %.2f^{o}",track[0],track[1]),"h");
  leg2->AddEntry("NULL","   ","h");
  leg2->AddEntry("NULL","#sigma_{T0}#oplus#sigma_{track}=40ps","h");
  leg2->AddEntry("NULL","#sigma_{TTS}#oplus#sigma_{Elec}=70ps","h");
  leg2->AddEntry("NULL","   ","h");
  leg2->AddEntry("NULL",Form("#pi/K separation = %.2f#sigma",separation),"h");

  Canvas->cd();
  pionHist->Draw("hist");
  kaonHist->Draw("histsame");
  leg1->Draw("same");
  leg2->Draw("same");
  return Canvas;
}
