#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>

void MultiTiers_TACor_sep(const char *rootname=""){

    TF1 *profilefit(TH2* Rt,double rbU,double rbt,double tRL,double tRR,double URL,double URR,char* name);
    TF1 *gausfit(TH1* h,int rbU,double* U_RL, double* U_RR);
    void twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1);
    //void MygStyle();
    //MygStyle();
    gStyle->SetOptFit(111);

    char name[100];
    char buff[1024];


    //***************************************************//
    //--------------Configuration-----------------------//
    //***************************************************//
    double tL = -1.e-9;
    double tR = 2.e-9;

    double uL = -15e3;
    double uR = 0; 

    int rbt = 2;
    int rbu = 80;

    // the range set After Correct
    //double L2 = -100e-12;
    //double R2 = 100e-12;
    int bint = (tR-tL)/1e-12;
    int binu = (uR-uL)/1;
    const int T = 1;
    const int iter = 5; //the number of iteration
    //---------------------------------------------------//	
    //***************************************************//


    //double RL,RR;
    //double U_RL,U_RR;
    //int binT = (init_tR-init_tL)/0.5e-12;
    //int binU = (init_UR-init_UL)/1;

    cout<<"Start excuate TA correction procedure =====>>>>"<<endl;
    //return;
    double thrd;
    for (int j = 14;j<15;j++)
    {
    //cout<<"<<---- Succeed excuating ---->>"<<endl; 
    //return;
        //thrd = (s+1)*0.2;
        sprintf(name,"%s",rootname);
        sprintf(buff,"%s.root",name);
        TFile *f = new TFile(buff,"READ");
        TTree *t = (TTree*)f->Get("data");
        double UL[T],UR[T],T0L[T],T0R[T],T0[T];
        double Utemp=0,Ttemp=0;
        double T0_cor=0;
        float TT[1000000]={0.};
        
        //t->SetBranchAddress("UL",&UL);
        t->SetBranchAddress("UR",UR);
        //t->SetBranchAddress("T0L",&T0L);
        t->SetBranchAddress("T0R",T0R);
        //t->SetBranchAddress("T0",T0);
        /*

*/
        TCanvas *c1 = new TCanvas("c1","",800,600);
        TCanvas *c2 = new TCanvas("c2","",800,600);
        TCanvas *c3 = new TCanvas("c3","",800,600);

        TCanvas *cU = new TCanvas("cU","",800,600);





        TH1D *h = new TH1D("h",";Time (s);Counts",bint,tL,tR);
        TH1D *hU = new TH1D("hU",";Amplitude (mV);Counts",binu,uL,uR);
        TH1D *hpfx;        
        TH2D *ht = new TH2D("ht",";Amp (mV);Time (s)",binu,uL,uR,bint,tL,tR);
        TF1 *fitAT = new TF1("fitAT","0",uL,uR);
        //t->Draw("T0>>h");

        for(int s=0;s<iter;s++){
            int N = 0;
            N = t->GetEntries();
            cout<<"the Entries is :"<<N<<endl;
            for(int i=0;i<N;i++){
                t->GetEntry(i);


                Utemp=0;
                Ttemp=0;
                for (int iT=0;iT<T;iT++){
                    Utemp+=UR[iT];
                    Ttemp+=T0R[iT];

                }
                Utemp=Utemp/T;
                Ttemp=Ttemp/T;
                T0_cor = fitAT->Eval(Utemp);
               
               if(Utemp<0&&Ttemp!=0){
                    if (s==0)
                    {
                        TT[i]=Ttemp-T0_cor;
                    }
                    else{
                        TT[i] = TT[i]-T0_cor;
                    }
                    //cout<<"T0_cor = "<<T0_cor<<endl;
                    //cout<<"TT = "<<TT[i]<<endl;
                    //cout<<"UU = "<<Utemp<<endl;

                    h->Fill(TT[i]);
                    hU->Fill(Utemp);
                    ht->Fill(Utemp,TT[i]);
                }
            }

            c1->cd();
            c1->Clear();
            //RL = h->GetMean();
            //cout<<"RL = "<<RL<<endl;
            //return;
            h->Draw();
            tL = h->GetMean()-5*h->GetRMS();
            tR = h->GetMean()+5*h->GetRMS();
            TF1 *fit1 = gausfit(h,rbt,&tL,&tR);
            //h->Draw();
            //h->Fit(fit1);
            sprintf(buff,"%s_T0_iter%d.png",name,s);
            c1->SaveAs(buff);
            //return;
            sprintf(buff,"%s_Tres.dat",name);
            ofstream output(buff,ios::app);
            output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;

            cU->cd();
            cU->Clear();
            hU->Draw();
            if(uL<hU->GetMean()-3*hU->GetRMS())
                uL = hU->GetMean()-3*hU->GetRMS();
            if(uR>hU->GetMean()+3*hU->GetRMS())
                uR = hU->GetMean()+3*hU->GetRMS();
            cout<<"uL="<<uL<<",uR="<<uR<<endl;
            //uR=-30;
            TF1 *fitU = gausfit(hU,rbu,&uL,&uR);
            //if(s==2) return;
            //hU->Draw();
            sprintf(buff,"%s_U_iter%d.png",name,s);
            cU->SaveAs(buff);	
            //return;
            cU->Close();

            c2->cd();
            c2->Clear();
            //t->Draw("T0:UR>>ht","","colz");
            sprintf(buff,"%s_T0U_iter%d",name,s);
            fitAT = profilefit(ht,rbu*2,rbt*2,tL,tR,uL,uR,buff);
            ht->Reset();
            h->Reset();
            hU->Reset();
        }
    }
}

TF1* gausfit(TH1* h,int rbU,double* U_RL, double* U_RR){



    TH1* hU = (TH1*)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    TF1 *fitU = new TF1("fitU","gaus",*U_RL,*U_RR);
    fitU->SetParameter(1,hU->GetMean());
    hU->Fit(fitU,"R");
	if(*U_RL<fitU->GetParameter(1)-5*fitU->GetParameter(2))
	*U_RL = fitU->GetParameter(1)-5*fitU->GetParameter(2);
	if(*U_RR>fitU->GetParameter(1)+5*fitU->GetParameter(2))
	*U_RR = fitU->GetParameter(1)+5*fitU->GetParameter(2);

    hU->GetXaxis()->SetRangeUser(*U_RL,*U_RR);
    hU->Fit(fitU);
    return fitU;

}


void twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1){
    //First fit for ensuring the rangement of histgram;
    TH1* h = (TH1*)ht->Clone();
    h->Rebin(rbt);
    TF1 *fit = new TF1("fit","gaus",*tRL,*tRR);        
    h->GetXaxis()->SetRangeUser(*tRL,*tRR);
    h->Fit(fit);
    double mean = fit->GetParameter(1);
    double sigma = TMath::Abs(fit->GetParameter(2));
    if(*tRL<mean-5*sigma||sigma>1)
    {*tRL = mean-5*sigma;
        *tRR = mean+5*sigma;}
    cout<<h->GetName()<<"\t"<<*tRL<<"\t"<<*tRR<<endl;

    h->GetXaxis()->SetRangeUser(*tRL,*tRR);

    TF1 *fit2 = new TF1("fit2","gaus(0)+gaus(3)",*tRL,*tRR);
    fit2->SetParNames("C_{TR}","#mu_{TR}","#sigma_{TR}","C_{nz}","#mu_{nz}","#sigma_{nz}");
    fit2->SetParameter(1,mean);
    fit2->SetParameter(2,sigma);
    fit2->SetParLimits(3,0,fit->GetParameter(0)*fac);
    fit2->SetParameter(4,mean);
    fit2->SetParameter(5,3*sigma);

    h->Fit(fit2,"","",mean-5*sigma,mean+5*sigma);
    TF1 *fit_tr = new TF1("fit_tr","gaus",*tRL,*tRR);
    fit_tr->SetParameter(0,fit2->GetParameter(0));
    fit_tr->SetParameter(1,fit2->GetParameter(1));
    fit_tr->SetParameter(2,fit2->GetParameter(2));
    TF1 *fit_nz = new TF1("fit_nz","gaus",*tRL,*tRR);
    fit_nz->SetParameter(0,fit2->GetParameter(3));
    fit_nz->SetParameter(1,fit2->GetParameter(4));
    fit_nz->SetParameter(2,fit2->GetParameter(5));
    fit_tr->SetLineColor(3);
    fit_tr->SetLineStyle(7);
    fit_nz->SetLineColor(3);
    fit_nz->SetLineStyle(7);
    fit_tr->Draw("same");
    fit_nz->Draw("same");

    if(fit2->GetParameter(3)>=fit2->GetParameter(0)&&TMath::Abs(fit2->GetParameter(5))>=0.05)
    {
        mean = fit2->GetParameter(4);
        sigma = TMath::Abs(fit2->GetParameter(5));
    }
    else if(TMath::Abs(fit2->GetParameter(2))>=0.05){
        mean = fit2->GetParameter(1);
        sigma = TMath::Abs(fit2->GetParameter(2));
    }

    //if(*tRL<mean-5*sigma)
    //{
    *tRL = mean-5*sigma;
    *tRR = mean+5*sigma;
    //}

    cout<<h->GetName()<<"\t"<<*tRL<<"\t"<<*tRR<<endl;
    h->GetXaxis()->SetRangeUser(*tRL,*tRR);
}

TF1* profilefit(TH2* Rt,double rbU,double rbt,double tRL,double tRR,double URL,double URR,char* name){
    char buff[1024];


    TCanvas *c5 = new TCanvas("c5","c5",1600,600);
    c5->Clear();
    c5->Divide(2,1);
    c5->cd(1);
    TH2* Qt = (TH2*)Rt->Clone();
    //TH2* Qt = (TH2*) Rt->Clone("tmp");

    Qt->Draw("colz");

    Qt->RebinX(rbU);
    Qt->RebinY(rbt);
    Qt->GetYaxis()->SetRangeUser(tRL,tRR);
    Qt->GetXaxis()->SetRangeUser(URL,URR);
    Qt->ProfileX();
    c5->cd(2);


    TH1* Qpfx = Qt->ProfileX();
    //Qpfx->Reset();


    //Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
    Qpfx->Draw();
    Qpfx->GetYaxis()->SetRangeUser(tRL,tRR);
    Qpfx->GetXaxis()->SetRangeUser(URL,URR);

    TF1* fitQt = new TF1("fitQt","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)",URL,URR);

    Qpfx->Fit(fitQt,"R");

    sprintf(buff,"%s.png",name);
    c5->SaveAs(buff);
    Qpfx->Reset();
    //delete Qt;
    //delete Qpfx;
    //delete c5;
    return fitQt;


}



/*
   void MygStyle(){
//gStyle->SetOptStat(10);
gStyle->SetOptFit(1111);
//gStyle->SetOptTitle(0);
gStyle->SetStatX(0.96);
gStyle->SetStatY(0.97);
gStyle->SetStatH(0.17);
gStyle->SetStatW(0.20);
gStyle->SetStatStyle(0);
gStyle->SetPadGridX(0);
gStyle->SetPadGridY(0);
*/
/*gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleFont(42,"T");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetTitleSize(0.05,"T");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetTitleOffset(1,"z");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetLabelColor(kRed,"XYZ");
  gStyle->SetTitleColor(kRed,"XYZ");
  gStyle->SetTitleColor(kRed,"T");
  gStyle->SetAxisColor(kRed,"XYZ");
  gStyle->SetFrameLineColor(kRed);
//gStyle->SetFrameFillStyle(3003);
//gStyle->SetFillColor(kGreen);
//gStyle->SetFrameFillColor(kGreen);*/
//gStyle->SetTitleX(0.55);

/*
   gStyle->SetFuncWidth(1);
   gStyle->SetHistLineWidth(1);
   gStyle->SetFuncColor(kRed);
   gStyle->SetEndErrorSize(0);
//TGaxis::SetMaxDigits(3);

}
*/	

