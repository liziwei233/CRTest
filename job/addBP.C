#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

TRandom3 r;

using namespace std;

void addBP(const char* rootname=""){

        char name[1024];        ;
        char buff[1024];
        sprintf(name,"%s",rootname);
        sprintf(buff,"%s.root",name);

        int N=0,temp=0,idN=0;
        vector<int>* IDR;
        vector<int>* IDL;
        vector<int>* phID;
        vector<double>* phX;
        vector<double>* phY;
        vector<double>* phZ;

        vector<double>* bpx_R;
        vector<double>* bpy_R;
        vector<double>* bpz_R;

        vector<double>* bpx_L;
        vector<double>* bpy_L;
        vector<double>* bpz_L;

        IDR = new vector<int>;
        IDL = new vector<int>;
        phID = new vector<int>;
        phX = new vector<double>;
        phY = new vector<double>;
        phZ = new vector<double>;

        bpx_R = new vector<double>;
        bpy_R = new vector<double>;
        bpz_R = new vector<double>;
        
        bpx_L = new vector<double>;
        bpy_L = new vector<double>;
        bpz_L = new vector<double>;

        TFile *f1 = new TFile(buff,"update");
        TTree *t1 = (TTree*)f1->Get("Run");
    
    // add branch to root file
        TBranch *bbpx_R=t1->Branch("PmtR.bpx",&bpx_R);
        TBranch *bbpy_R=t1->Branch("PmtR.bpy",&bpy_R);
        TBranch *bbpz_R=t1->Branch("PmtR.bpz",&bpz_R);
        TBranch *bbpx_L=t1->Branch("PmtL.bpx",&bpx_L);
        TBranch *bbpy_L=t1->Branch("PmtL.bpy",&bpy_L);
        TBranch *bbpz_L=t1->Branch("PmtL.bpz",&bpz_L);

        t1->SetBranchAddress("PmtR.trackID",&IDR);
        t1->SetBranchAddress("PmtL.trackID",&IDL);
        t1->SetBranchAddress("ph.ID", &phID);
        t1->SetBranchAddress("ph.x", &phX);
        t1->SetBranchAddress("ph.y", &phY);
        t1->SetBranchAddress("ph.z", &phZ);

        N = t1->GetEntries();
        cout<<"Entries = "<<N<<endl;

        for(int i = 0; i < N; i++){

            //-----------initial----------------------//
            IDR->clear();
            IDL->clear();
            phX->clear();
            phY->clear();
            phZ->clear();
            phID->clear();

            bpx_R->clear();
            bpy_R->clear();
            bpz_R->clear();
        
            bpx_L->clear();
            bpy_L->clear();
            bpz_L->clear();

            t1->GetEntry(i);
            
            temp = IDR->size();
            //cout<<"counterR = "<< temp <<endl;
            for(int k=0;k<temp;k++){
               
                //find the birthplace of  photons those hit right PMT.
                idN = (*IDR)[k];
                for(int p=idN-20;p<idN;p++)
                {
                    if((*phID)[p]==idN) {
                        //cout<<"RID= "<<idN<<endl;
                        //cout<<"phID[p]= "<<(*phID)[p]<<",p= "<<p<<endl;
                        bpx_R->push_back((*phX)[p]);
                        bpy_R->push_back((*phY)[p]);
                        bpz_R->push_back((*phZ)[p]);
                    }
                }

            }
            
            temp = IDL->size();
            for(int k=0;k<temp;k++){
               
                //find the birthplace of  photons those hit right PMT.
                idN = (*IDL)[k];
                for(int p=idN-20;p<idN;p++)
                {
                    if((*phID)[p]==idN) {
                        //cout<<"LID= "<<idN<<endl;
                        //cout<<"phID[p]= "<<(*phID)[p]<<",p= "<<p<<endl;
                        bpx_L->push_back((*phX)[p]);
                        bpy_L->push_back((*phY)[p]);
                        bpz_L->push_back((*phZ)[p]);
                    }
                }

            }
            bbpx_R->Fill();
            bbpy_R->Fill();
            bbpz_R->Fill();
            bbpx_L->Fill();
            bbpy_L->Fill();
            bbpz_L->Fill();
        }
        t1->Print();
        t1->Write();
        delete f1;
}