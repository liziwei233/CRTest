#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "vector"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

class DIRC {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Double_t        PrimaryEnergy;
    vector<double>  *TrackDx;
    vector<double>  *TrackDy;
    vector<double>  *TrackDz;
    vector<double>  *TrackPx;
    vector<double>  *TrackPy;
    vector<double>  *TrackPz;
    Int_t           PhotoNu;
    vector<int>     *Type;
    vector<int>     *SectorID;
    vector<int>     *ChannelX;
    vector<int>     *ChannelY;
    vector<double>  *PhotonKinetic;
    vector<double>  *WaveLength;
    vector<double>  *GlobalTime;
    vector<double>  *LocalTime;
    vector<double>  *TrackLength;
    vector<double>  *ThetaX;
    vector<double>  *ThetaZ;
    vector<double>  *ThetaC;
    vector<double>  *PhotonDx;
    vector<double>  *PhotonDy;
    vector<double>  *PhotonDz;
    vector<double>  *PhotonPx;
    vector<double>  *PhotonPy;
    vector<double>  *PhotonPz;

    // List of branches
    TBranch        *b_PrimaryEnergy;   //!
    TBranch        *b_TrackDx;   //!
    TBranch        *b_TrackDy;   //!
    TBranch        *b_TrackDz;   //!
    TBranch        *b_TrackPx;   //!
    TBranch        *b_TrackPy;   //!
    TBranch        *b_TrackPz;   //!
    TBranch        *b_PhotoNu;   //!
    TBranch        *b_Type;   //!
    TBranch        *b_SectorID;   //!
    TBranch        *b_ChannelX;   //!
    TBranch        *b_ChannelY;   //!
    TBranch        *b_PhotonKinetic;   //!
    TBranch        *b_WaveLength;   //!
    TBranch        *b_GlobalTime;   //!
    TBranch        *b_LocalTime;   //!
    TBranch        *b_TrackLength;   //!
    TBranch        *b_ThetaX;   //!
    TBranch        *b_ThetaZ;   //!
    TBranch        *b_ThetaC;   //!
    TBranch        *b_PhotonDx;   //!
    TBranch        *b_PhotonDy;   //!
    TBranch        *b_PhotonDz;   //!
    TBranch        *b_PhotonPx;   //!
    TBranch        *b_PhotonPy;   //!
    TBranch        *b_PhotonPz;   //!

    DIRC(TTree *tree=0);
    virtual ~DIRC();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    //my function
    TFile* f;
};

DIRC::DIRC(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile* fIn = new TFile("../data/opticalphoton.root");
    fIn->GetObject("DIRC",tree);

  }
  Init(tree);
}

DIRC::~DIRC()
{
  f->Close();
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t DIRC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t DIRC::LoadTree(Long64_t entry)
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

void DIRC::Init(TTree *tree)
{
  // Set object pointer
  TrackDx = 0;
  TrackDy = 0;
  TrackDz = 0;
  TrackPx = 0;
  TrackPy = 0;
  TrackPz = 0;
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

  fChain->SetBranchAddress("PrimaryEnergy", &PrimaryEnergy, &b_PrimaryEnergy);
  fChain->SetBranchAddress("TrackDx", &TrackDx, &b_TrackDx);
  fChain->SetBranchAddress("TrackDy", &TrackDy, &b_TrackDy);
  fChain->SetBranchAddress("TrackDz", &TrackDz, &b_TrackDz);
  fChain->SetBranchAddress("TrackPx", &TrackPx, &b_TrackPx);
  fChain->SetBranchAddress("TrackPy", &TrackPy, &b_TrackPy);
  fChain->SetBranchAddress("TrackPz", &TrackPz, &b_TrackPz);
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

  f = new TFile("LUT.root","recreate");
}

Bool_t DIRC::Notify()
{
  return kTRUE;
}

void DIRC::Show(Long64_t entry)
{
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t DIRC::Cut(Long64_t entry)
{
  return 1;
}

void DIRC::Loop()
{
  if (fChain == 0) return;

  TProfile* pHvsY[7];
  for(int i=0;i<7;i++)
  {
    if(i==0) pHvsY[i]=new TProfile("Height_vs_ChannelY", "Height vs ChannelY",128,0,128);
    else if(i<4) pHvsY[i]=new TProfile(Form("positive_%d",i),Form("Z=(%d,%d),positive",15-5*i,20-5*i),128,0,128);
    else  pHvsY[i]=new TProfile(Form("negtive_%d",i),Form("Z=(%d,%d),negtive",5*i-20,5*i-15),128,0,128);
    pHvsY[i]->SetMinimum(125);
    pHvsY[i]->SetMaximum(215);
    pHvsY[i]->SetLineColor(i+1);
    pHvsY[i]->SetLineWidth(2);
    pHvsY[i]->SetMarkerColor(i+1);
    pHvsY[i]->SetMarkerStyle(20);
    pHvsY[i]->GetXaxis()->SetTitle("ChannelY");
    pHvsY[i]->GetYaxis()->SetTitle("Height");
  }
  TCanvas* cHvsY = new TCanvas("Height_vs_ChannelY","Height vs ChannelY");
  TProfile* pThZvsY = new TProfile("ThetaZ_vs_ChannelY","ThetaZ vs ChannelY",128,0,128);
  pThZvsY->SetMinimum(0.15);
  pThZvsY->SetMaximum(0.9);
  pThZvsY->SetLineColor(4);
  pThZvsY->SetLineWidth(2);
  pThZvsY->SetMarkerColor(4);
  pThZvsY->SetMarkerStyle(20);
  pThZvsY->GetXaxis()->SetTitle("ChannelY");
  pThZvsY->GetYaxis()->SetTitle("#theta_{z}(rad)");

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    long ventries = Type -> size();
    for(long i=0; i<ventries; i++)
    {
      double Z=1350.94015-PhotonPz->at(i);
      int direction = PhotonDz->at(i)>0? 1:-1;
      int j=Z/5.*direction+4.;
      double TrackHeight = TrackLength->at(i)/TMath::Sqrt(1./TMath::Power(TMath::Cos(ThetaX->at(i)),2)+1./TMath::Power(TMath::Cos(ThetaZ->at(i)),2)-1.);
      pHvsY[j]->Fill(ChannelY->at(i),TrackHeight);
      pHvsY[0]->Fill(ChannelY->at(i),TrackHeight);
      pThZvsY->Fill(ChannelY->at(i),ThetaZ->at(i));
    }
  }

  ofstream ofH("LUT/Height.txt");
  ofstream ofThZ("LUT/ThetaZ.txt");
  for(int i=1;i<=128;i++)
  {
    ofH<<pHvsY[0]->GetBinContent(i)<<endl;
    ofThZ<<pThZvsY->GetBinContent(i)<<endl;
  }
  ofH.close();
  ofThZ.close();

  cHvsY->cd();
  for(int i=1;i<7;i++)
    pHvsY[i]->Draw("same");

  f->cd();
  cHvsY->Write();
  pHvsY[0]->Write();
  pThZvsY->Write();
}

void makeLUT()
{
  DIRC d;
  d.Loop();
}
