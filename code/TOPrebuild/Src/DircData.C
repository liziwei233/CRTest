#include "DircData.h"

using namespace TMath;
using namespace std;


DircData::DircData(string in,string out,TTree *tree) : fChain(0),input(in),output(out) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile* fIn = new TFile(input.c_str());
    fIn->GetObject("DIRC",tree);
  }
  Init(tree);
}

DircData::~DircData()
{
  f->Close();
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t DircData::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t DircData::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void DircData::Init(TTree *tree)
{
  // Set object pointer
  TrackDx = 0;
  TrackDy = 0;
  TrackDz = 0;
  TrackPx = 0;
  TrackPy = 0;
  TrackPz = 0;
  FlightLength = 0;
  Type = 0;
  SectorID = 0;
  ChannelX = 0;
  ChannelY = 0;
  PhotonKinetic = 0;
  WaveLength = 0;
  GlobalTime = 0;
  LocalTime = 0;
  TrackLength = 0;
  ThetaX = 0;
  ThetaZ = 0;
  ThetaC = 0;
  PhotonDx = 0;
  PhotonDy = 0;
  PhotonDz = 0;
  PhotonPx = 0;
  PhotonPy = 0;
  PhotonPz = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("PrimaryMomentum", &PrimaryMomentum, &b_PrimaryMomentum);
  fChain->SetBranchAddress("TrackDx", &TrackDx, &b_TrackDx);
  fChain->SetBranchAddress("TrackDy", &TrackDy, &b_TrackDy);
  fChain->SetBranchAddress("TrackDz", &TrackDz, &b_TrackDz);
  fChain->SetBranchAddress("TrackPx", &TrackPx, &b_TrackPx);
  fChain->SetBranchAddress("TrackPy", &TrackPy, &b_TrackPy);
  fChain->SetBranchAddress("TrackPz", &TrackPz, &b_TrackPz);
  fChain->SetBranchAddress("FlightLength", &FlightLength, &b_FlightLength);
  fChain->SetBranchAddress("PhotoNu", &PhotoNu, &b_PhotoNu);
  fChain->SetBranchAddress("Type", &Type, &b_Type);
  fChain->SetBranchAddress("SectorID", &SectorID, &b_SectorID);
  fChain->SetBranchAddress("ChannelX", &ChannelX, &b_ChannelX);
  fChain->SetBranchAddress("ChannelY", &ChannelY, &b_ChannelY);
  fChain->SetBranchAddress("PhotonKinetic", &PhotonKinetic, &b_PhotonKinetic);
  fChain->SetBranchAddress("WaveLength", &WaveLength, &b_WaveLength);
  fChain->SetBranchAddress("GlobalTime", &GlobalTime, &b_GlobalTime);
  fChain->SetBranchAddress("LocalTime", &LocalTime, &b_LocalTime);
  fChain->SetBranchAddress("TrackLength", &TrackLength, &b_TrackLength);
  fChain->SetBranchAddress("ThetaX", &ThetaX, &b_ThetaX);
  fChain->SetBranchAddress("ThetaZ", &ThetaZ, &b_ThetaZ);
  fChain->SetBranchAddress("ThetaC", &ThetaC, &b_ThetaC);
  fChain->SetBranchAddress("PhotonDx", &PhotonDx, &b_PhotonDx);
  fChain->SetBranchAddress("PhotonDy", &PhotonDy, &b_PhotonDy);
  fChain->SetBranchAddress("PhotonDz", &PhotonDz, &b_PhotonDz);
  fChain->SetBranchAddress("PhotonPx", &PhotonPx, &b_PhotonPx);
  fChain->SetBranchAddress("PhotonPy", &PhotonPy, &b_PhotonPy);
  fChain->SetBranchAddress("PhotonPz", &PhotonPz, &b_PhotonPz);
  Notify();

  f = new TFile(output.c_str(),"recreate");
}

Bool_t DircData::Notify()
{
  return kTRUE;
}

void DircData::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t DircData::Cut(Long64_t entry)
{
  return 1;
}

void DircData::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
  }
}

void DircData::Sort()
{
  int tmp[216][4]={0};
  for(int i=0;i<GlobalTime->size();i++)
  {
    if(tmp[ChannelX->at(i)][ChannelY->at(i)]==0) 
      tmp[ChannelX->at(i)][ChannelY->at(i)]=i;
    else
    {
      PhotoNu--;
      if(GlobalTime->at(i)<GlobalTime->at(tmp[ChannelX->at(i)][ChannelY->at(i)]))
      {
        (*GlobalTime)[tmp[ChannelX->at(i)][ChannelY->at(i)]]=0;
        tmp[ChannelX->at(i)][ChannelY->at(i)]=i;
      }
      else (*GlobalTime)[i]=0; 
    }
  }
}
