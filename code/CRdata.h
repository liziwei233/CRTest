//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 27 10:35:18 2021 by ROOT version 6.20/00
// from TTree data/restore analysed data  from G4
// found on file: 101RBdata.root
//////////////////////////////////////////////////////////

#ifndef CRdata_h
#define CRdata_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace TMath;
using namespace std;

// Header file for the classes stored in the TTree if any.

class CRdata {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //CRsysRBData     *data;
   vector<int>     T0photonid;
   vector<double>  T0photonE;
   vector<double>  T0photonTOP;
   vector<double>  T0photonx;
   vector<double>  T0photony;
   vector<double>  T0photonz;
   vector<double>  T0photont;
   Int_t           T0detid;
   Double_t        T0detx;
   Double_t        T0dety;
   Double_t        T0detz;
   Double_t        T0dett;
   vector<int>     FTOFphotonid;
   vector<double>  FTOFphotonE;
   vector<double>  FTOFphotonTOP;
   vector<double>  FTOFphotonx;
   vector<double>  FTOFphotony;
   vector<double>  FTOFphotonz;
   vector<double>  FTOFphotont;
   Int_t           FTOFdetid;
   Double_t        FTOFdetx;
   Double_t        FTOFdety;
   Double_t        FTOFdetz;
   Double_t        FTOFdett;
   Double_t        Trackerdetx[4];
   Double_t        Trackerdety[4];
   Double_t        Trackerdetz[4];
   Double_t        Trackerdett[4];
   Double_t        CRE;
   Double_t        CRpx;
   Double_t        CRpy;
   Double_t        CRpz;
   Double_t        CRtheta;
   Double_t        T0detRBx;
   Double_t        T0detRBy;
   Double_t        T0detRBz;
   Int_t           T0eleid[4];
   Double_t        T0elethrd[4];
   Double_t        T0eleU[4];
   Double_t        T0eletot[4];
   Double_t        T0elethtime1[4];
   Double_t        T0elethtime2[4];
   Double_t        T0elefittot[4];
   Double_t        T0elefittime1[4];
   Double_t        T0elefittime2[4];
   Double_t        FTOFdetRBx;
   Double_t        FTOFdetRBy;
   Double_t        FTOFdetRBz;
   Int_t           FTOFeleid[128];
   Double_t        FTOFelethrd[128];
   Double_t        FTOFeleU[128];
   Double_t        FTOFeletot[128];
   Double_t        FTOFelethtime1[128];
   Double_t        FTOFelethtime2[128];
   Double_t        FTOFelefittot[128];
   Double_t        FTOFelefittime1[128];
   Double_t        FTOFelefittime2[128];
   Double_t        CRRBpx;
   Double_t        CRRBpy;
   Double_t        CRRBpz;
   Double_t        CRRBtheta;

   // List of branches
   TBranch        *b_data_T0photonid;   //!
   TBranch        *b_data_T0photonE;   //!
   TBranch        *b_data_T0photonTOP;   //!
   TBranch        *b_data_T0photonx;   //!
   TBranch        *b_data_T0photony;   //!
   TBranch        *b_data_T0photonz;   //!
   TBranch        *b_data_T0photont;   //!
   TBranch        *b_data_T0detid;   //!
   TBranch        *b_data_T0detx;   //!
   TBranch        *b_data_T0dety;   //!
   TBranch        *b_data_T0detz;   //!
   TBranch        *b_data_T0dett;   //!
   TBranch        *b_data_FTOFphotonid;   //!
   TBranch        *b_data_FTOFphotonE;   //!
   TBranch        *b_data_FTOFphotonTOP;   //!
   TBranch        *b_data_FTOFphotonx;   //!
   TBranch        *b_data_FTOFphotony;   //!
   TBranch        *b_data_FTOFphotonz;   //!
   TBranch        *b_data_FTOFphotont;   //!
   TBranch        *b_data_FTOFdetid;   //!
   TBranch        *b_data_FTOFdetx;   //!
   TBranch        *b_data_FTOFdety;   //!
   TBranch        *b_data_FTOFdetz;   //!
   TBranch        *b_data_FTOFdett;   //!
   TBranch        *b_data_Trackerdetx;   //!
   TBranch        *b_data_Trackerdety;   //!
   TBranch        *b_data_Trackerdetz;   //!
   TBranch        *b_data_Trackerdett;   //!
   TBranch        *b_data_CRE;   //!
   TBranch        *b_data_CRpx;   //!
   TBranch        *b_data_CRpy;   //!
   TBranch        *b_data_CRpz;   //!
   TBranch        *b_data_CRtheta;   //!
   TBranch        *b_data_T0detRBx;   //!
   TBranch        *b_data_T0detRBy;   //!
   TBranch        *b_data_T0detRBz;   //!
   TBranch        *b_data_T0eleid;   //!
   TBranch        *b_data_T0elethrd;   //!
   TBranch        *b_data_T0eleU;   //!
   TBranch        *b_data_T0eletot;   //!
   TBranch        *b_data_T0elethtime1;   //!
   TBranch        *b_data_T0elethtime2;   //!
   TBranch        *b_data_T0elefittot;   //!
   TBranch        *b_data_T0elefittime1;   //!
   TBranch        *b_data_T0elefittime2;   //!
   TBranch        *b_data_FTOFdetRBx;   //!
   TBranch        *b_data_FTOFdetRBy;   //!
   TBranch        *b_data_FTOFdetRBz;   //!
   TBranch        *b_data_FTOFeleid;   //!
   TBranch        *b_data_FTOFelethrd;   //!
   TBranch        *b_data_FTOFeleU;   //!
   TBranch        *b_data_FTOFeletot;   //!
   TBranch        *b_data_FTOFelethtime1;   //!
   TBranch        *b_data_FTOFelethtime2;   //!
   TBranch        *b_data_FTOFelefittot;   //!
   TBranch        *b_data_FTOFelefittime1;   //!
   TBranch        *b_data_FTOFelefittime2;   //!
   TBranch        *b_data_CRRBpx;   //!
   TBranch        *b_data_CRRBpy;   //!
   TBranch        *b_data_CRRBpz;   //!
   TBranch        *b_data_CRRBtheta;   //!

   CRdata(string in="data/DIRC.root",
        string out="result/DircRec.root",
        TTree *tree=0);
   virtual ~CRdata();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   string input;
   string output;
   TFile* fout;
};

#endif

#ifdef CRdata_cxx
CRdata::CRdata(string in,string out,TTree *tree) : fChain(0),input(in),output(out) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
    TFile* fIn = new TFile(input.c_str());
    fIn->GetObject("DIRC",tree);
  }
      
   Init(tree);
}

CRdata::~CRdata()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CRdata::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CRdata::LoadTree(Long64_t entry)
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

void CRdata::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("T0photonid", &T0photonid, &b_data_T0photonid);
   fChain->SetBranchAddress("T0photonE", &T0photonE, &b_data_T0photonE);
   fChain->SetBranchAddress("T0photonTOP", &T0photonTOP, &b_data_T0photonTOP);
   fChain->SetBranchAddress("T0photonx", &T0photonx, &b_data_T0photonx);
   fChain->SetBranchAddress("T0photony", &T0photony, &b_data_T0photony);
   fChain->SetBranchAddress("T0photonz", &T0photonz, &b_data_T0photonz);
   fChain->SetBranchAddress("T0photont", &T0photont, &b_data_T0photont);
   fChain->SetBranchAddress("T0detid", &T0detid, &b_data_T0detid);
   fChain->SetBranchAddress("T0detx", &T0detx, &b_data_T0detx);
   fChain->SetBranchAddress("T0dety", &T0dety, &b_data_T0dety);
   fChain->SetBranchAddress("T0detz", &T0detz, &b_data_T0detz);
   fChain->SetBranchAddress("T0dett", &T0dett, &b_data_T0dett);
   fChain->SetBranchAddress("FTOFphotonid", &FTOFphotonid, &b_data_FTOFphotonid);
   fChain->SetBranchAddress("FTOFphotonE", &FTOFphotonE, &b_data_FTOFphotonE);
   fChain->SetBranchAddress("FTOFphotonTOP", &FTOFphotonTOP, &b_data_FTOFphotonTOP);
   fChain->SetBranchAddress("FTOFphotonx", &FTOFphotonx, &b_data_FTOFphotonx);
   fChain->SetBranchAddress("FTOFphotony", &FTOFphotony, &b_data_FTOFphotony);
   fChain->SetBranchAddress("FTOFphotonz", &FTOFphotonz, &b_data_FTOFphotonz);
   fChain->SetBranchAddress("FTOFphotont", &FTOFphotont, &b_data_FTOFphotont);
   fChain->SetBranchAddress("FTOFdetid", &FTOFdetid, &b_data_FTOFdetid);
   fChain->SetBranchAddress("FTOFdetx", &FTOFdetx, &b_data_FTOFdetx);
   fChain->SetBranchAddress("FTOFdety", &FTOFdety, &b_data_FTOFdety);
   fChain->SetBranchAddress("FTOFdetz", &FTOFdetz, &b_data_FTOFdetz);
   fChain->SetBranchAddress("FTOFdett", &FTOFdett, &b_data_FTOFdett);
   fChain->SetBranchAddress("Trackerdetx[4]", Trackerdetx, &b_data_Trackerdetx);
   fChain->SetBranchAddress("Trackerdety[4]", Trackerdety, &b_data_Trackerdety);
   fChain->SetBranchAddress("Trackerdetz[4]", Trackerdetz, &b_data_Trackerdetz);
   fChain->SetBranchAddress("Trackerdett[4]", Trackerdett, &b_data_Trackerdett);
   fChain->SetBranchAddress("CRE", &CRE, &b_data_CRE);
   fChain->SetBranchAddress("CRpx", &CRpx, &b_data_CRpx);
   fChain->SetBranchAddress("CRpy", &CRpy, &b_data_CRpy);
   fChain->SetBranchAddress("CRpz", &CRpz, &b_data_CRpz);
   fChain->SetBranchAddress("CRtheta", &CRtheta, &b_data_CRtheta);
   fChain->SetBranchAddress("T0detRBx", &T0detRBx, &b_data_T0detRBx);
   fChain->SetBranchAddress("T0detRBy", &T0detRBy, &b_data_T0detRBy);
   fChain->SetBranchAddress("T0detRBz", &T0detRBz, &b_data_T0detRBz);
   fChain->SetBranchAddress("T0eleid[4]", T0eleid, &b_data_T0eleid);
   fChain->SetBranchAddress("T0elethrd[4]", T0elethrd, &b_data_T0elethrd);
   fChain->SetBranchAddress("T0eleU[4]", T0eleU, &b_data_T0eleU);
   fChain->SetBranchAddress("T0eletot[4]", T0eletot, &b_data_T0eletot);
   fChain->SetBranchAddress("T0elethtime1[4]", T0elethtime1, &b_data_T0elethtime1);
   fChain->SetBranchAddress("T0elethtime2[4]", T0elethtime2, &b_data_T0elethtime2);
   fChain->SetBranchAddress("T0elefittot[4]", T0elefittot, &b_data_T0elefittot);
   fChain->SetBranchAddress("T0elefittime1[4]", T0elefittime1, &b_data_T0elefittime1);
   fChain->SetBranchAddress("T0elefittime2[4]", T0elefittime2, &b_data_T0elefittime2);
   fChain->SetBranchAddress("FTOFdetRBx", &FTOFdetRBx, &b_data_FTOFdetRBx);
   fChain->SetBranchAddress("FTOFdetRBy", &FTOFdetRBy, &b_data_FTOFdetRBy);
   fChain->SetBranchAddress("FTOFdetRBz", &FTOFdetRBz, &b_data_FTOFdetRBz);
   fChain->SetBranchAddress("FTOFeleid[128]", FTOFeleid, &b_data_FTOFeleid);
   fChain->SetBranchAddress("FTOFelethrd[128]", FTOFelethrd, &b_data_FTOFelethrd);
   fChain->SetBranchAddress("FTOFeleU[128]", FTOFeleU, &b_data_FTOFeleU);
   fChain->SetBranchAddress("FTOFeletot[128]", FTOFeletot, &b_data_FTOFeletot);
   fChain->SetBranchAddress("FTOFelethtime1[128]", FTOFelethtime1, &b_data_FTOFelethtime1);
   fChain->SetBranchAddress("FTOFelethtime2[128]", FTOFelethtime2, &b_data_FTOFelethtime2);
   fChain->SetBranchAddress("FTOFelefittot[128]", FTOFelefittot, &b_data_FTOFelefittot);
   fChain->SetBranchAddress("FTOFelefittime1[128]", FTOFelefittime1, &b_data_FTOFelefittime1);
   fChain->SetBranchAddress("FTOFelefittime2[128]", FTOFelefittime2, &b_data_FTOFelefittime2);
   fChain->SetBranchAddress("CRRBpx", &CRRBpx, &b_data_CRRBpx);
   fChain->SetBranchAddress("CRRBpy", &CRRBpy, &b_data_CRRBpy);
   fChain->SetBranchAddress("CRRBpz", &CRRBpz, &b_data_CRRBpz);
   fChain->SetBranchAddress("CRRBtheta", &CRRBtheta, &b_data_CRRBtheta);
   Notify();

  fout = new TFile(output.c_str(),"recreate");
   
}

Bool_t CRdata::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CRdata::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CRdata::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CRdata_cxx
