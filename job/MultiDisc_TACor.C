#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>

struct EVENT{
	double time[20],A;
};

void MultiDisc_TACor(const char *rootname=""){

    TF1 *profilefit(TH2* Rt,double rbU,double rbt,double tRL,double tRR,double URL,double URR,char* name);
    TF1 *gausfit(TH1* h,int rbU,double* U_RL, double* U_RR);
    void twoguasfit(TH1 *ht,double* tRL,double* tRR,double fac=0.4, int rbt=1);
    TF1* CH2Correction(TTree* t1,EVENT* A,EVENT* B,double* p,int rbU,int rbt,char* name,double fac,int iter);
    TF1* CH1Correction(TTree* t1,EVENT* A,EVENT* B,double* p,int rbU,int rbt,char* name,double fac,int iter);
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
    int rbu = 16;

    // the range set After Correct
    //double L2 = -100e-12;
    //double R2 = 100e-12;
    int bint = (tR-tL)/1e-12;
    int binu = (uR-uL)/1;
    const int T = 100;
   // double fac;
    const int iter = 5; //the number of iteration
    //---------------------------------------------------//	
    //***************************************************//


    //double RL,RR;
    //double U_RL,U_RR;
    //int binT = (init_tR-init_tL)/0.5e-12;
    //int binU = (init_UR-init_UL)/1;

    cout<<"Start excuate TA correction procedure =====>>>>"<<endl;
    //return;
    
    //cout<<"<<---- Succeed excuating ---->>"<<endl; 
    //return;
        //thrd = (s+1)*0.2;
        sprintf(name,"%s",rootname);
        sprintf(buff,"%s.root",name);
        TFile *f = new TFile(buff,"READ");
        TTree *t = (TTree*)f->Get("data");

        EVENT L,R;
     
        //double ULtemp=0,TLtemp=0;
        //double URtemp=0,TRtemp=0;
       // double T0L_cor=0;
        //double T0R_cor=0;
        double p[20];
       // float TT[1000000]={0.};
        
        t->SetBranchAddress("UL",&L.A);
        t->SetBranchAddress("UR",$R.A);
        t->SetBranchAddress("T0L",L.time);
        t->SetBranchAddress("T0R",R.time);
        t->SetBranchAddress("percent",p);

        sprintf(buff,"L_%s",name);
        CH1Correction(t,&L,&R,p,2*rbu,4*rbt,buff,fac,iter);

        sprintf(buff,"R_%s",name);
        CH1Correction(t,&R,&L,p,2*rbu,4*rbt,buff,fac,iter);

        sprintf(buff,"L-R_%s",name);
	    CH2Correction(t,&L,&R,p,2*rbU,4*rbt,buff,fac,iter);

        
       
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
TF1* CH3Correction(TTree* t1,EVENT* A,EVENT* B,EVENT* MCP,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 




	TCanvas* c4 = new TCanvas("c4","c4",800,600); 	
	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-10, 10);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",200,0,1.6e3,2e3,-10,10);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=2;
	double tR=6;
	double QL=0;
	double QR=1.6e3;

	double TAcor=0;
	double Q1=0,Q2=0;
	double T[1000000]={0};
	const int chN=2;
	EVENT** ch[chN];

	ch[0]=&A;
	ch[1]=&B;
	//ch[2]=&MCP;
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;

		//correction of A1
		for (int h= 0;h<chN;h++){
		tL=-2;
		tR=2;
		
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			//Q1=(*A).Q;
			//Q2=(*B).Q;
			if (h-1<0)
				TAcor=fitAT->Eval((**ch[chN-1]).Q);
			else
				TAcor=fitAT->Eval((**ch[h-1]).Q);
			
			//TAcor=fitAT->Eval(Q2);
			if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
				//if(1){
				if((*A).Q<600&&(*B).Q<600){
				if (s==0) 
					T[i]=((*A).p20+(*B).p20)/2-(*MCP).p20-TAcor;
				else 
					T[i]=T[i]-TAcor;


				ht->Fill(T[i]);		
				hAT->Fill((**ch[h]).Q,T[i]);
			}
			}

		}


		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		twoguasfit(ht,&tL,&tR,fac,rbt);
		sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
		c4->SaveAs(buff);

		sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
		}

	}


	return fitAT;

}

TF1* CH2Correction(TTree* t1,EVENT* A,EVENT* B,double* p,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 




	TCanvas* c4 = new TCanvas("c4","c4",800,600); 	
	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-2, 2);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",600,-3e3,0,2e3,-2,2);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=-2;
	double tR=2;
	double QL=0;
	double QR=1.6e3;

	double TAcor=0;
	double Q1=0,Q2=0;
	double T[1000000]={0};
	const int chN=2;
	EVENT** ch[chN];

	ch[0]=&A;
	ch[1]=&B;
	//ch[2]=&MCP;
    double L=sizeof(p) / sizeof(p[0]);
    for(int j=0;j<L;j++){
    
	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;

		//correction of A1
		for (int h= 0;h<chN;h++){
		tL=-2;
		tR=2;
		
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			//Q1=(*A).Q;
			//Q2=(*B).Q;
			if (h-1<0)
				TAcor=fitAT->Eval((**ch[chN-1]).A[j]);
			else
				TAcor=fitAT->Eval((**ch[h-1]).A[j]);
			
			//TAcor=fitAT->Eval(Q2);
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4){
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*A).rise<4&&(*B).rise<4&&(*MCP).Q>0.3){
				//if(1){
				//if(1){
				//if((*A).A>0.3&&(*B).A>0.3){
                if(p[j]==0) return fitAT;
                else if((*A).A[j]<0&&(*A).time!=0&&(*B).A[j]<0&&(*B).time!=0){
    
				if (s==0&&h==0) 
					//T[i]=(*A).p20-(*B).p20-TAcor;
					T[i]=(*A).time-(*B).time;
				else 
					T[i]=T[i]-TAcor;


				ht->Fill(T[i]);		
				hAT->Fill((**ch[h]).A[j],T[i]);
			
			}

		}


		c4->cd();
		c4->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;

		//twoguasfit(ht,&tL,&tR,fac,rbt);
        gausfit(ht,rbt,&tL,&tR);
		//sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
        sprintf(buff,"%s_th%g_TR_cor%d.png",name,p[j],s);
		c4->SaveAs(buff);
        sprintf(buff,"%s_th%gAt_pfx_cor%d",name,p[j],s);
		//sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
		}
        
	}
    fitAT=new TF1("fitAT","0",-1e4,1e4);
    }

	return fitAT;

}

TF1* CH1Correction(TTree* t1,EVENT* A,EVENT* B,double* p,int rbU,int rbt,char* name,double fac,int iter){
	//Correction QT
	char buff[1024]; 


	TCanvas* c1 = new TCanvas("c1","c1",800,600); 	


	sprintf(buff,"%s_ht",name);
	TH1D *ht = new TH1D(buff,"Time Resolution;T (ps);Counts",2e3,-2, 2);




	sprintf(buff,"%s_hAT",name);
	TH2D *hAT = new TH2D(buff,"",600,-3e3,0,2e3,-2,2);


	TF1 *fitAT=new TF1("fitAT","0",-1e4,1e4);





	double tL=0;
	double tR=2;
	double QL=0;
	double QR=2e3;

	double TAcor=0;
	double Q=0;
	double T[1000000]={0};
    double L=sizeof(p) / sizeof(p[0]);
    for(int j=0;j<L;j++){

	for(int s =0;s<iter;s++){



		int N = t1->GetEntries();
		cout<< "Total entries is :"<<N<<endl;
		for(int i=0; i<N; i++)
		{
			t1->GetEntry(i);

			Q=(*A).A[j];
			TAcor=fitAT->Eval(Q);
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
			//if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
            if(p[j]==0) return fitAT;
            else if((*A).A[p]<0&&(*A).time!=0){

				//if((*A).Q<600){
				//if(1){
					if (s==0) 
						//T[i]=(*A).p20-(*MCP).p20-TAcor;
						T[i]=(*A).A;
					else 
						T[i]=T[i]-TAcor;


					ht->Fill(T[i]);		
//				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
					hAT->Fill(Q,T[i]);
				}
			
		}

//cout<<"<====process check====>"<<endl;

		c1->cd();
		c1->Clear();
		ht->Draw();
		//return;
		//tRL=-5;
		//tRR=5;
		tL=ht->GetMean()-3*ht->GetRMS();
		tR=ht->GetMean()+3*ht->GetRMS();
		//cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
		//cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
		//if (s==1) return;
		//twoguasfit(ht,&tL,&tR,fac,rbt);
        gausfit(ht,rbt,&tL,&tR);

		sprintf(buff,"%s_th%g_TR_cor%d.png",name,p[j],s);
		c1->SaveAs(buff);

		sprintf(buff,"%s_th%gAt_pfx_cor%d",name,p[j],s);

		fitAT=profilefit(hAT,rbU,rbt*2,tL,tR,QL,QR,buff);

		ht->Reset();
		hAT->Reset();
	}
    fitAT=new TF1("fitAT","0",-1e4,1e4);
}

	return fitAT;

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

