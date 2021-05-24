#ifndef DircData_h
#define DircData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom3.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace TMath;
using namespace std;
class DircData {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Double_t        PrimaryMomentum;
    vector<double>  *TrackDx;
    vector<double>  *TrackDy;
    vector<double>  *TrackDz;
    vector<double>  *TrackPx;
    vector<double>  *TrackPy;
    vector<double>  *TrackPz;
    vector<double>  *FlightLength;
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
    TBranch        *b_PrimaryMomentum;   //!
    TBranch        *b_TrackDx;   //!
    TBranch        *b_TrackDy;   //!
    TBranch        *b_TrackDz;   //!
    TBranch        *b_TrackPx;   //!
    TBranch        *b_TrackPy;   //!
    TBranch        *b_TrackPz;   //!
    TBranch        *b_FlightLength;   //!
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

    DircData(string in="data/DIRC.root",
        string out="result/DircRec.root",
        TTree *tree=0);
    virtual ~DircData();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual void     Sort();

    //my function

    string input;
    string output;
    TFile* f;
};

#endif
