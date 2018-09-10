#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>

void MultiTiers_TACor_sep(const char *rootname=""){
	
	
	TF1* gausfit(TH1* hU,int rbU,float* U_RL,float* U_RR);
	TF1* prfxfit(TH2* ht,int rbU,int rbt,float RL, float RR, float U_RL, float U_RR, char* name);
	//void MygStyle();
	//MygStyle();
	gStyle->SetOptFit(111);
	
	char name[100];
	char buff[1024];


//***************************************************//
//--------------Configuration-----------------------//
//***************************************************//

	
	int rbt = 5;
	int rbU = 12;

	float RL = -2.e-9;
	float RR = 16.e-9;
	float U_RL = -6e3;
	float U_RR = 0;
	// the range set After Correct
	float L2 = -500e-12;
	float R2 = 500e-12;
	int binT2 = 2e3;
//---------------------------------------------------//	
//***************************************************//



	int binT = (RR-RL)/10e-12;
	int binU = (U_RR-U_RL)/1;
	
	
	const int T = 4;
	float thrd;
	for (int s = 14;s<15;s++)
	{
	//thrd = (s+1)*0.2;
	sprintf(name,"%s",rootname);
	sprintf(buff,"%s.root",name);
	TFile *f = new TFile(buff,"READ");
	TTree *t = (TTree*)f->Get("data");
	double UL[T],UR[T],T0L[T],T0R[T],T0[T];
	double T0_cor;
	double UU=0.,TT=0.;

	//t->SetBranchAddress("UL",&UL);
	t->SetBranchAddress("UR",UR);
	t->SetBranchAddress("T0L",T0L);
	t->SetBranchAddress("T0R",T0R);
	t->SetBranchAddress("T0",T0);
	/*
	
	*/
	TCanvas *c1 = new TCanvas("c1","",800,600);
	TCanvas *c2 = new TCanvas("c2","",800,600);
	TCanvas *c3 = new TCanvas("c3","",800,600);
	
	TCanvas *cU = new TCanvas("cU","",800,600);
	

	
	
	
	TH1D *h = new TH1D("h","",binT, RL, RR);
	TH1D *hU = new TH1D("hU","",binU, U_RL, U_RR);
	TH2D *ht = new TH2D("ht","",binU, U_RL, U_RR, binT, RL, RR);
	//t->Draw("T0>>h");

	int N = 0;
	N = t->GetEntries();
	for(int i=0;i<N;i++){
		t->GetEntry(i);
		TT=0;
		UU=0;
		for (int iT=0;iT<T;iT++){
			TT+=T0[iT];
			UU+=UR[iT];
		}
		TT=TT/T;
		UU=UU/T;

		h->Fill(TT);
		hU->Fill(UU);
		ht->Fill(UU,TT);
	}
	
	c1->cd();
	TF1* fit1 = gausfit(h, rbt, &RL, &RR);
	sprintf(buff,"%s_T0.png",name);
	c1->SaveAs(buff);
	//return;
	sprintf(buff,"%s_Tres.dat",name);
	ofstream output(buff,ios::app);
	output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;
	
	cU->cd();
	gausfit(hU, rbU, &U_RL, &U_RR);
	sprintf(buff,"%s_U.png",name);
	cU->SaveAs(buff);	
	//return;
	cU->Close();
	
	
	sprintf(buff,"%s_beforeCor.png",name);
	TF1* fit2 = prfxfit(ht, rbU, rbt, RL, RR, U_RL, U_RR, buff);
	//return;
	//h->Reset();

//*****************************************************************//
//--------------------After Correction-----------------------------//
//*****************************************************************//	
	
	ht->Reset();
	ht->GetYaxis()->SetRangeUser(L2,R2);
	//ht->Draw();
	//return;
	h->Reset();
	h->GetXaxis()->SetRangeUser(L2,R2);
	//TH2D *ht_cor= new TH2D("ht_cor","",binU,U_RL,U_RR,binT2,L2,R2);
	//hpfx->GetYaxis()->SetRangeUser(-0.3e-9,0.3e-9);
	
	//TH1D *h_cor = new TH1D("h_cor","",3600,-1e-9,2.6e-9);
	//nt n=t->GetEntries();
	for(int i=0; i<N; i++)
	{
		t->GetEntry(i);
		TT = 0;
		UU = 0;
		for (int iT=0;iT<T;iT++){
			TT+=T0[iT];
			UU+=UR[iT];
		}
		TT=TT/T;
		UU=UU/T;
		
		T0_cor = fit2->Eval(UU);
		//cout<<"T0-T0_cor="<<T0<<"-"<<T0_cor<<"="<<T0-T0_cor<<endl;
		
		h->Fill(TT-T0_cor);
		ht->Fill(UU,TT-T0_cor);
		//cout<<"U:T0-T0_cor="<<UR<<"\t"<<T0-T0_cor<<endl;
	}
	c1->cd();
	c1->Clear();
	fit1 = gausfit(h, 1, &L2, &R2);
	//return;
	//h->GetXaxis()->SetRangeUser(-0.1e-9,0.1e-9);
	sprintf(buff,"%s_T0_AfterCorrect.png",name);
	c1->SaveAs(buff);
	
	//output("timeres.dat",ios::app);
	output<<fit1->GetParameter(2)<<"\t"<<fit1->GetParError(2)<<endl;

	sprintf(buff,"%s_AfterCor",name);	
	prfxfit(ht, 1, 1, L2, R2, U_RL, U_RR, buff);

	}
	
}

TF1* gausfit(TH1* hU,int rbU,float* U_RL, float* U_RR){
	
	
	
	hU->Draw();
	hU->Rebin(rbU);
	TF1 *fitU = new TF1("fitU","gaus",*U_RL,*U_RR);
	fitU->SetParameter(1,hU->GetMean());
	hU->Fit(fitU,"R");
	*U_RL = fitU->GetParameter(1)-5*fitU->GetParameter(2);
	*U_RR = fitU->GetParameter(1)+5*fitU->GetParameter(2);

	hU->GetXaxis()->SetRangeUser(*U_RL,*U_RR);
	hU->Fit(fitU);
	return fitU;

}

TF1* prfxfit(TH2* ht,int rbU,int rbt,float RL, float RR, float U_RL, float U_RR, char* name){
	char buff[1024];
	TCanvas *c2 = new TCanvas("c2","c2",1600,600);	
	c2->Clear();
	c2->Divide(2,1);
	c2->cd(1);
	//t->Draw("T0:UR>>ht","","colz");
	ht->Draw("colz");
	ht->RebinX(rbU);
	ht->RebinY(rbt);
	ht->GetYaxis()->SetRangeUser(RL,RR);
	ht->GetXaxis()->SetRangeUser(U_RL,U_RR);
	//return;

	c2->cd(2);
	TH1* hpfx=ht->ProfileX();
	hpfx->GetYaxis()->SetRangeUser(RL,RR);
	hpfx->GetXaxis()->SetRangeUser(U_RL,U_RR);
	hpfx->Draw();
	
	TF1 *fit2 = new TF1("fit2","[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)",U_RL,U_RR);
	hpfx->Fit(fit2,"R");

	sprintf(buff,"%s.png",name);
	c2->SaveAs(buff);
	hpfx->Reset();
	return fit2;
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
	
