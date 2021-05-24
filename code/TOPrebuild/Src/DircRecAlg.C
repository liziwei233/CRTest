#include "DircRecAlg.h"

using namespace TMath;
using namespace std;

DircRecAlg::DircRecAlg(string in, string out, TTree* tree):DircRecSvc(in,out,tree)
{
  mass["pion"]   =139.5701;
  mass["kaon"]   =493.677;
  mass["proton"] = 938.27231;
  PDF["pion"]   = new TF1("PDF_pion",  "gaus(0)+0.05",0,10);
  PDF["kaon"]   = new TF1("PDF_kaon",  "gaus(0)+0.05",0,10);
  double resolution = 0.1; // ns
  PDF["pion"] -> SetParameter(0,Sqrt(0.5/Pi())/resolution);
  PDF["pion"] -> SetParameter(1,5.224);
  PDF["pion"] -> SetParameter(2,resolution);
  PDF["kaon"] -> SetParameter(0,Sqrt(0.5/Pi())/resolution);
  PDF["kaon"] -> SetParameter(1,5.368);
  PDF["kaon"] -> SetParameter(2,resolution);

  fNtuple = new TTree("DircRec","Rec data of DIRC");
  book();

  hist1["ThetaC"]               = new TH1F("ThetaC",              "ThetaC",               400, 0.6, 1);
  hist1["TrackLength"]          = new TH1F("TrackLength",         "TrackLength",          400, 0,   1300);
  hist1["RecTrackLengthHyPi"]   = new TH1F("RecTrackLengthHyPi",  "RecTrackLengthHyPi",   400, 0,   1300);
  hist1["DeltaTrackLengthHyPi"] = new TH1F("DeltaTrackLengthHyPi","DeltaTrackLengthHyPi", 400, -50, 50);
  hist2["TrackLengthHyPi"]      = new TH2F("TrackLength2DHyPi",   "TrackLength2DHyPi",    400, 0,   1300, 400, 0,   1300);
  hist1["RecFlightTimeHyPi"]    = new TH1F("RecFlightTimeHyPi",   "RecFlightTimeHyPi",    400, 4,   7.5);
  hist1["MeanFlightTimeHyPi"]   = new TH1F("MeanFlightTimeHyPi",  "MeanFlightTimeHyPi",   400, 4,   7.5);
  hist1["RecTrackLengthHyK"]    = new TH1F("RecTrackLengthHyK",   "RecTrackLengthHyK",    400, 0,   1300);
  hist1["DeltaTrackLengthHyK"]  = new TH1F("DeltaTrackLengthHyK", "DeltaTrackLengthHyK",  400, -50, 50);
  hist2["TrackLengthHyK"]       = new TH2F("TrackLength2DHyK",    "TrackLength2DHyK",     400, 0,   1300, 400, 0,   1300);
  hist1["RecFlightTimeHyK"]     = new TH1F("RecFlightTimeHyK",    "RecFlightTimeHyK",     400, 4,   7.5);
  hist1["MeanFlightTimeHyK"]    = new TH1F("MeanFlightTimeHyK",   "MeanFlightTimeHyK",    400, 4,   7.5);
  hist1["ProbL"]                = new TH1F("ProbL",               "ProbL",                400, -100, 100);
  hist1["PhotoNu"]              = new TH1F("PhotoNu",             "PhotoNu",              1000, 0,  1000);
  hist2["TvsX"]                 = new TH2F("TvsX",                "TvsX",                 72, 0,   72,  400, -40,   60);
}

DircRecAlg::~DircRecAlg()
{
  f->cd();
  fNtuple->Write();
  for(std::map<string,TH1F*>::iterator it=hist1.begin();it!=hist1.end();it++)
    it->second->Write();
  for(std::map<string,TH2F*>::iterator it=hist2.begin();it!=hist2.end();it++)
    it->second->Write();
  for(std::map<string,TF1*>::iterator it=PDF.begin();it!=PDF.end();it++)
    it->second->Write();
}
void DircRecAlg::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  for (Long64_t tentry=0; tentry<nentries;tentry++) 
  {
    fChain->GetEntry(tentry);
    Sort();

    if(FlightLength->size()==0) continue;
    PDF["pion"]->SetParameter(1,FlightLength->at(0)/TMath::C()*Sqrt(1+Power(mass["pion"]/PrimaryMomentum,2))*1e6);
    PDF["kaon"]->SetParameter(1,FlightLength->at(0)/TMath::C()*Sqrt(1+Power(mass["kaon"]/PrimaryMomentum,2))*1e6);
    SetTrackHit();
    hist1["PhotoNu"]->Fill(PhotoNu);

    TOFHyPi.clear();
    TOFHyK.clear();
    long ventries = Type -> size();
    ProbL = 0;
    for(long i=0; i<ventries; i++)
    {
      /// some cut
      //if(Type->at(i)!=1) continue;
      if(GlobalTime->at(i)==0) continue;

      /// set hits
      SetPhotonHit(i);
      TruThetaC=ThetaC->at(i);
      TruPropLength=TrackLength->at(i);


      /// Fill fata
      hist1["ThetaC"]->Fill(TruThetaC);
      hist1["TrackLength"]->Fill(TruPropLength);

      /// reconstruction
      Reconstruction(mass["pion"]);
      ProbL += Log10(PDF["pion"]->Eval(BestFlightTime));
      RecPropLengthHyPi = BestPropLength;
      RecFlightTimeHyPi = BestFlightTime;
      TOFHyPi.push_back(RecFlightTimeHyPi);
      hist1["RecTrackLengthHyPi"]->Fill(BestPropLength);
      hist1["DeltaTrackLengthHyPi"]->Fill(BestPropLength-TruPropLength);
      hist2["TrackLengthHyPi"]->Fill(TruPropLength,BestPropLength);
      hist1["RecFlightTimeHyPi"]->Fill(BestFlightTime);

      Reconstruction(mass["kaon"]);
      ProbL -= Log10(PDF["kaon"]->Eval(BestFlightTime));
      RecPropLengthHyK = BestPropLength;
      RecFlightTimeHyK = BestFlightTime;
      TOFHyK.push_back(RecFlightTimeHyK);
      hist1["RecTrackLengthHyK"]->Fill(BestPropLength);
      hist1["DeltaTrackLengthHyK"]->Fill(BestPropLength-TruPropLength);
      hist2["TrackLengthHyK"]->Fill(TruPropLength,BestPropLength);
      hist1["RecFlightTimeHyK"]->Fill(BestFlightTime);

      hist2["TvsX"]->Fill(ChannelX->at(i),GlobalTime->at(i));

      fNtuple->Fill();
    }
    hist1["ProbL"]->Fill(ProbL);
    int cut = 2;
    if(TOFHyPi.size()>cut)
    {
      sort(TOFHyPi.begin(),TOFHyPi.end());
      double sum=0, count=0;
      for(int i=cut/2;i<TOFHyPi.size()-cut/2;i++)
      {
        if(TOFHyPi.at(i)==0) continue;
        sum +=TOFHyPi.at(i);
        count += 1;
      }
      if(count!=0) hist1["MeanFlightTimeHyPi"]->Fill(sum/count);
    }
    if(TOFHyK.size()>cut)
    {
      sort(TOFHyK.begin(),TOFHyK.end());
      double sum=0,count=0;
      for(int i=cut/2;i<TOFHyK.size()-cut/2;i++)
      {
        if(TOFHyK.at(i)==0) continue;
        sum +=TOFHyK.at(i);
        count += 1;
      }
      if(count!=0) hist1["MeanFlightTimeHyK"]->Fill(sum/count);
    }
  }
}

void DircRecAlg::book()
{
  fNtuple->Branch("ChX",        &ChX);
  fNtuple->Branch("ChY",        &ChY);
  fNtuple->Branch("ThetaC",     &TruThetaC);
  fNtuple->Branch("PropLength", &TruPropLength);
  fNtuple->Branch("RecPropLengthHyPi", &RecPropLengthHyPi);
  fNtuple->Branch("RecFlightTimeHyPi", &RecFlightTimeHyPi);
  fNtuple->Branch("RecPropLengthHyK",  &RecPropLengthHyK);
  fNtuple->Branch("RecFlightTimeHyK",  &RecFlightTimeHyK);
  fNtuple->Branch("ProbL",             &ProbL);
}

