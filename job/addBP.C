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
        double theta=0,phi=0;
        int initial=0;
        int limit=0;
        //int counter=0;
        //int counter2=0;
        TVector3 MD_L;
        TVector3 MD_R;
        
        vector<int>* IDR;
        vector<int>* IDL;
        vector<int>* phID;
        vector<double>* phX;
        vector<double>* phY;
        vector<double>* phZ;

        vector<double>* phPX;
        vector<double>* phPY;
        vector<double>* phPZ;

        vector<double>* bpx_R;
        vector<double>* bpy_R;
        vector<double>* bpz_R;

        vector<double>* bpx_L;
        vector<double>* bpy_L;
        vector<double>* bpz_L;
        
        vector<double>* bpPx_R;
        vector<double>* bpPy_R;
        vector<double>* bpPz_R;

        vector<double>* bpPx_L;
        vector<double>* bpPy_L;
        vector<double>* bpPz_L;
        
        vector<double>* bptheta_L;
        vector<double>* bpphi_L;
        
        vector<double>* bptheta_R;
        vector<double>* bpphi_R;


        IDR = new vector<int>;
        IDL = new vector<int>;
        phID = new vector<int>;
        phX = new vector<double>;
        phY = new vector<double>;
        phZ = new vector<double>;
        phPX = new vector<double>;
        phPY = new vector<double>;
        phPZ = new vector<double>;

        bpx_R = new vector<double>;
        bpy_R = new vector<double>;
        bpz_R = new vector<double>;
        bpx_L = new vector<double>;
        bpy_L = new vector<double>;
        bpz_L = new vector<double>;
        
        bpPx_L = new vector<double>;
        bpPy_L = new vector<double>;
        bpPz_L = new vector<double>;

        bpPx_R = new vector<double>;
        bpPy_R = new vector<double>;
        bpPz_R = new vector<double>;
        
        bptheta_R = new vector<double>;
        bpphi_R = new vector<double>;
        bptheta_L = new vector<double>;
        bpphi_L = new vector<double>;
        
        TFile *f1 = new TFile(buff,"update");
        TTree *t1 = (TTree*)f1->Get("Run");
    
    // add branch to root file
        TBranch *bbpx_R=t1->Branch("PmtR.bpx",&bpx_R);
        TBranch *bbpy_R=t1->Branch("PmtR.bpy",&bpy_R);
        TBranch *bbpz_R=t1->Branch("PmtR.bpz",&bpz_R);
        TBranch *bbpx_L=t1->Branch("PmtL.bpx",&bpx_L);
        TBranch *bbpy_L=t1->Branch("PmtL.bpy",&bpy_L);
        TBranch *bbpz_L=t1->Branch("PmtL.bpz",&bpz_L);
        
        /*
        TBranch *bbpPx_R=t1->Branch("PmtR.bpPx",&bpPx_R);
        TBranch *bbpPy_R=t1->Branch("PmtR.bpPy",&bpPy_R);
        TBranch *bbpPz_R=t1->Branch("PmtR.bpPz",&bpPz_R);
        TBranch *bbpPx_L=t1->Branch("PmtL.bpPx",&bpPx_L);
        TBranch *bbpPy_L=t1->Branch("PmtL.bpPy",&bpPy_L);
        TBranch *bbpPz_L=t1->Branch("PmtL.bpPz",&bpPz_L);
        
        TBranch *bbptheta_L=t1->Branch("PmtL.theta",&bptheta_L);
        TBranch *bbpphi_L=t1->Branch("PmtL.phi",&bpphi_L);
        TBranch *bbptheta_R=t1->Branch("PmtR.theta",&bptheta_R);
        TBranch *bbpphi_R=t1->Branch("PmtR.phi",&bpphi_R);
*/
        t1->SetBranchAddress("PmtR.trackID",&IDR);
        t1->SetBranchAddress("PmtL.trackID",&IDL);
        t1->SetBranchAddress("ph.ID", &phID);
        t1->SetBranchAddress("ph.x", &phX);
        t1->SetBranchAddress("ph.y", &phY);
        t1->SetBranchAddress("ph.z", &phZ);
        t1->SetBranchAddress("ph.px", &phPX);
        t1->SetBranchAddress("ph.py", &phPY);
        t1->SetBranchAddress("ph.pz", &phPZ);

        N = t1->GetEntries();
        cout<<"Entries = "<<N<<endl;

        for(int i = 0; i < N; i++){

            //-----------initial----------------------//
            IDR->clear();
            IDL->clear();
            phX->clear();
            phY->clear();
            phZ->clear();
            phPX->clear();
            phPY->clear();
            phPZ->clear();
            phID->clear();

            bpx_R->clear();
            bpy_R->clear();
            bpz_R->clear();
        
            bpx_L->clear();
            bpy_L->clear();
            bpz_L->clear();
            
            bpPx_R->clear();
            bpPy_R->clear();
            bpPz_R->clear();
        
            bpPx_L->clear();
            bpPy_L->clear();
            bpPz_L->clear();
        
            bptheta_R->clear();
            bpphi_R->clear();
            bptheta_L->clear();
            bpphi_L->clear();

            t1->GetEntry(i);
            
            temp = IDR->size();
            //cout<<"counterR = "<< temp <<endl;
            for(int k=0;k<temp;k++){
               
                //find the birthplace of  photons those hit right PMT.
                idN = (*IDR)[k];
                if(idN-20<0) initial=0;
                else initial=idN-20;
                if(phID->size()<idN) limit=phID->size();
                else limit=idN;
                for(int p=initial;p<limit;p++)
                {
                    if((*phID)[p]==idN) {
                        //cout<<"RID= "<<idN<<endl;
                        //cout<<"phID[p]= "<<(*phID)[p]<<",p= "<<p<<endl;
                        bpx_R->push_back((*phX)[p]);
                        bpy_R->push_back((*phY)[p]);
                        bpz_R->push_back((*phZ)[p]);
                        bpPx_R->push_back((*phPX)[p]);
                        bpPy_R->push_back((*phPY)[p]);
                        bpPz_R->push_back((*phPZ)[p]);
                        MD_R.SetXYZ((*phPX)[p],(*phPY)[p],(*phPZ)[p]);
                        theta=MD_R.Theta()*TMath::RadToDeg();
                        phi=MD_R.Phi()*TMath::RadToDeg();
                        bptheta_R->push_back(theta);
                        bpphi_R->push_back(phi);

                    }
                }

            }
            if(temp!=bpx_R->size()) cout<<"entry$="<<i<<", R counts= "<<temp<<", bp counts= "<<bpx_R->size()<<endl;
            temp = IDL->size();
            //counter2+=temp;
            for(int k=0;k<temp;k++){
               
                //find the birthplace of  photons those hit right PMT.
                idN = (*IDL)[k];
                for(int p=idN-20;p<phID->size();p++)
                {
                    if((*phID)[p]==idN) {
                        //cout<<"p ="<<p<<", phid= "<<(*phID)[p]<<", IDL="<<(*IDL)[k]<<", idN="<<idN<<endl;
                       // cout<<"Entry$= "<<i<<", temp="<<temp<<", k="<<k<<endl;
                        //cout<<"counter="<<counter<<endl;
                        //cout<<"counter2="<<counter2<<endl;
                        //cout<<"LID= "<<idN<<endl;
                        //cout<<"phID[p]= "<<(*phID)[p]<<",p= "<<p<<endl;
                        bpx_L->push_back((*phX)[p]);
                        bpy_L->push_back((*phY)[p]);
                        bpz_L->push_back((*phZ)[p]);
                        bpPx_L->push_back((*phPX)[p]);
                        bpPy_L->push_back((*phPY)[p]);
                        bpPz_L->push_back((*phPZ)[p]);
                        MD_L.SetXYZ((*phPX)[p],(*phPY)[p],(*phPZ)[p]);
                        theta=MD_L.Theta()*TMath::RadToDeg();
                        phi=MD_L.Phi()*TMath::RadToDeg();
                        bptheta_L->push_back(theta);
                        bpphi_L->push_back(phi);
                        //counter++;
                    }
                }

            }
            if(temp!=bpx_L->size()) cout<<"R counts= "<<temp<<", bp counts= "<<bpx_L->size()<<endl;
            bbpx_R->Fill();
            bbpy_R->Fill();
            bbpz_R->Fill();
            bbpx_L->Fill();
            bbpy_L->Fill();
            bbpz_L->Fill();
            
            /*
            bbpPx_R->Fill();
            bbpPy_R->Fill();
            bbpPz_R->Fill();
            bbpPx_L->Fill();
            bbpPy_L->Fill();
            bbpPz_L->Fill();
            bbptheta_L->Fill();
            bbpphi_L->Fill();
            bbptheta_R->Fill();
            bbpphi_R->Fill();
            */
        }
        //TBranch *b1 = t1->GetBranch("ph.E");
        //t1->GetListOfBranches()->Remove(b1);

        //t1->Print();
        t1->Write("", TObject::kOverwrite);
        delete f1;
}
