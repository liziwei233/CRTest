#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

TRandom3 r;

using namespace std;

Double_t outputfunc(Double_t x, vector<double> par);
Double_t response(Double_t x, Double_t par[7]);


Double_t response(Double_t x, Double_t par[7]){
	Double_t G = par[0];
	Double_t Q = par[1];
	Double_t Ca = par[2];
	Double_t C = par[3];
	Double_t R = par[4];
	Double_t Z = par[5];
	Double_t rise = par[6];

	Double_t val = 0;
	Double_t tt = 550e12;
	Double_t t0;
	Double_t a,b,beta,gamma;
	
	// how to accumulate these single photon signal?
	// which distribution does t0 sample? 
	//!! get t0 by result of simulation!!
	
	beta = ((R+Z)*C+2.0*Ca*R)/(2.0*C*Ca*R*Z);
	gamma = 1.0/(2.0*C*Ca*R*Z);
	a = -(beta+TMath::Sqrt(beta*beta-4.0*gamma))/2.0;
	b = -(beta-TMath::Sqrt(beta*beta-4.0*gamma))/2.0;
	
	val = -0.5*G*Q*(b*TMath::Exp(b*x+0.5*b*b*rise*rise)*TMath::Erfc((-b*rise-x/rise)/TMath::Sqrt(2.0))-a*TMath::Exp(a*x+0.5*a*a*rise*rise)*TMath::Erfc((-a*rise-x/rise)/TMath::Sqrt(2.0)))/(Ca*(b-a))*1e3;
	
	
	
	return val;


	
}

Double_t outputfunc(Double_t x, vector<double> par){


	//************SiPM S13360-30-50***************
	//***Rise Time 1.5ns
	//***Fall Time 40ns
	//***Gain 1.7e6
	//***SPE Umax =0.71mV
	//***TTS =170ps
	//*************************************************

	
	Double_t val = 0;
	Double_t tts = 0;
    double SPEpar[7];
	
	r.SetSeed(0);
	tts = r.Gaus(0,170e-12);
	//tts = 0;

    SPEpar[0]=1.7e6;  //Gain
	SPEpar[1]=1.6e-19; //e
	//SPEpar[2]=150e-12;  //Ca  ??
	//SPEpar[3]=3e-12; //C
	SPEpar[2]=35e-11;  //Ca  ??
	SPEpar[3]=80e-11; //C
	SPEpar[4]=10e3;    //R   ??
	SPEpar[5]=50;      //Z
	SPEpar[6]=0.65e-9;  //rise time 1.5ns
	//int N;
	//N=sizeof(par)/sizeof(par[0]);
	int n=0;
	for (int n=0;n<par.size();n++){
	//while(par[n]>5e-9){
		if (x-par.at(n)<-30e-9){
		val+=0;}
		else val+=response(x-tts-par.at(n),SPEpar);
		//cout<<"    [-] x : par : val "<<x<<"\t"<<par[n]<<"\t"<<val<<endl;
		
	}
	//cout<<"n = "<<n<<endl;
	
	return val;
	


}


void MultiTiersOutputfun_SiPM(const char *rootname=""){
	
	
	TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);

	TLatex* DrawMyLatex(char* text, Double_t x=0.65, Double_t y=0.5, Int_t textFont=62, Size_t textSize=0.05, Color_t colorIndex=2);

	//void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);

	
	/*--------------SPE----------------
	TF1 *myFun = new TF1("myFun",response,RL,RR,7);
	myFun->SetParNames("Gain","Q","#C_{a}","C","R","Z", "#t_{Rise}");
	myFun->SetParameter(0,2e5);  //Gain
	myFun->SetParameter(1,1.6e-19); //e
	myFun->SetParameter(2,10e-12);  //Ca  ??
	myFun->SetParameter(3,3e-12); //C
	myFun->SetParameter(4,10e3);    //R   ??
	myFun->SetParameter(5,50);      //Z
	myFun->SetParameter(6,80e-12);  //rise time
	*/
	

	gStyle->SetOptFit(1111);

	//TRandom3 r;
	//r.SetSeed(0);
	char name[1024];
	char buff[1024];

	const int T = 4;//amounts of Tiers 

	//Double_t parR[500]={};
	//Double_t parL[500]={};
	vector<vector<double> > parR(T);


	//vector<double> parR2;
	vector<vector<double> > parL(T);
	//vector<double> parL2;

	Double_t RL = 0e-9;
	Double_t RR = 15e-9;
	int binNum=0;
	binNum = (RR-RL)/1e-12;

    Int_t range = 5e3; //3ps/Sample
	Double_t thrd = -2;//Umax = -2
	double Rate=0;
	
	bool flagR=0,flagL=0;
	double xT0_L[T]={0},xT0_R[T]={0},xT0[T]={0};
	double UL[T]={0},UR[T]={0};
	const int certain=5;	
	
	vector<double>* TR;
	vector<double>* TL;
	vector<int>* IDL;
	vector<int>* IDR;
	TL = new vector<double>;
	TR = new vector<double>;
	IDL = new vector<int>;
	IDR = new vector<int>;
	//count = new vector<int>;
	int N=0,temp=0;
	double Temp = 0.;
	Double_t x[5000]={};
	Double_t yR[T][5000]={};
	Double_t yL[T][5000]={};
	

	//TString a = TString(rootname);
	//TString b;
	//b.Append(a,a.Length()-5);
	sprintf(name,"%s",rootname);
	sprintf(buff,"%s.root",name);

	TFile *f1 = new TFile(buff,"READ");
	TTree *t1 = (TTree*)f1->Get("Run");

	t1->SetBranchAddress("PmtR.t",&TR);
	t1->SetBranchAddress("PmtL.t",&TL);
	t1->SetBranchAddress("PmtL.id",&IDL);
	t1->SetBranchAddress("PmtR.id",&IDR);

	//sprintf(name,"Thrd_%g",abs(thrd));	
	
	sprintf(buff,"%sdata.root",name);

	TFile *f2 = new TFile(buff,"RECREATE");
	TTree *t2 = new TTree("data","restore analysed data  from G4");
	t2->Branch("UL",UL,"UL[4]/D");
	t2->Branch("UR",UR,"UR[4]/D");
	t2->Branch("T0L",xT0_L,"T0L[4]/D");
	t2->Branch("T0R",xT0_R,"T0R[4]/D");
	t2->Branch("T0",xT0,"T0[4]/D");	
	//for(int s = 0; s<4;s++){

	
	
	double t_L = RL;
	double t_R = RR;
	
	
	//f1->cd();
	
	
	//TF1 *myFun;
	TH1D *h[2];
	h[0] = new TH1D("hR","",binNum,RL,RR);
	h[1] = new TH1D("hL","",binNum,RL,RR);

	TH1D *hSIG = new TH1D("hSIG","PMT's total signal",binNum,RL,RR);
	TH1D *hSig[2];
	hSig[0] = new TH1D("hRSig","",binNum,RL,RR);
	hSig[1] = new TH1D("hLSig","",binNum,RL,RR);

	
	TF1 *fSIG = new TF1("fSIG","gaus",RL,RR);
	TF1 *fSig[2];
	fSig[0] = new TF1("fRSig","gaus",RL,RR);
	fSig[1] = new TF1("fLSig","gaus",RL,RR);

	N = t1->GetEntries();
	cout<<"Entries = "<<N<<endl;


	//count->clear();
	//for(int i = certain; i < certain+1; i++){
	for(int i = 0; i < 25; i++){
		
		//-----------initial----------------------//
		TL->clear();
		TR->clear();
		IDL->clear();
		IDR->clear();

		h[0]->Reset();
		h[1]->Reset();

		//parR.clear();
		//parL.clear();
		vector<vector<double> >(T).swap(parR);
		vector<vector<double> >(T).swap(parL);
	

		//memset(parL,0,sizeof(parL));
		//memset(parR,0,sizeof(parR));
		
	//for(int i = certain; i < certain+1; i++){
		//par[i]=r.Gaus(2.4e-9,0.5e-9);
		//par[i]=4e-9;
		t1->GetEntry(i);
		temp = TR->size();
		
		//cout<<"counterR = "<< temp <<endl;
		//myFun = new TF1("myFun",outputfunc,RL,RR,temp);

		
		
		for(int k=0;k<temp;k++){
			//cout<< T[][k] <<endl;
			Temp = (*TR)[k]*1e-9;
			
			for(int iT=0;iT<T;iT++){
				//cout<< "enter the LOOP"<<endl;
				
				//parR[iT].push_back(0);
				//cout<<parR[iT][k]<<endl;
				//

				if((*IDR)[k]==iT) {
					parR[iT].push_back(Temp);
				//parR[iT].push_back(0);
			//parR[k]=8.3e-9;
			//cout<<" [+] parR1 "<<k<<"\t = "<<Temp<<endl;
			//cout<<"par"<<k<<" = "<<parR[iT][k]<<endl;
					}
			//return;
				if(iT==0) h[0]->Fill(Temp);
				}
				//return;
			
			
			//myFun->SetParameter(k,par[k]);
			}
			//return;
			
		temp = TL->size();	
		//cout<<"counterL = "<< temp <<endl;
		for(int k=0;k<temp;k++){
			//cout<< T[][k] <<endl;
			Temp = (*TL)[k]*1e-9;
			for(int iT=0; iT<T; iT++){
				if((*IDL)[k]==iT) parL[iT].push_back(Temp);
			//parR[k]=8.3e-9;
			//cout<<" [+] parR1 "<<k<<"\t = "<<Temp<<endl;
			//cout<<"par"<<k<<" = "<<par[k]<<endl;
				if(iT==0) h[1]->Fill(Temp);
				}			
			//myFun->SetParameter(k,par[k]);
			}
			
			
		for(int q=0;q<T;q++){
		flagR = 1;
		flagL = 1;
		
		xT0_R[q] = 0;
		//xT0_R[1] = 0;
		xT0_L[q] = 0;
		//xT0_L[1] = 0;
		xT0[q] = 0;
		//xT0[1] = 0;

		//cout<<"hello"<<endl;
		//cout<<"parL.size() = "<<parL.size()<<endl;
		//cout<<"parR.size() = "<<parR.size()<<endl;
		for(int j=0;j<range;j++){
			x[j]=(RR-RL)/range*j+RL;

			yL[q][j]=outputfunc(x[j],parL[q]);
			yR[q][j]=outputfunc(x[j],parR[q]);

			if(yR[q][j]<thrd&&flagR) {xT0_R[q]=x[j];flagR = false;}
			if(yL[q][j]<thrd&&flagL) {xT0_L[q]=x[j];flagL = false;}
			//cout<<"[+] x"<<j<<":y"<<j<<"=\t"<<x[j]<<"\t"<<yR[j]<<"\t"<<yL[j]<<endl;
		}

		
		UR[q] = TMath::MinElement(range,yR[q]);
		UL[q] = TMath::MinElement(range,yL[q]);
		/*
		cout<<"[+] PMT_Right Messages :"<<endl;
		//Float_t U0 = TMath::MinElement(range,yR);
		//Float_t t0 = g->GetY(U0);
		cout<<"		[-]U0 = "<<UR<<"mV"<<endl;
	
		cout<<"[+] PMT_Left Messages :"<<endl;
		//U0 = TMath::MinElement(range,yL);
		//Float_t t0 = g->GetY(U0);
		cout<<"		[-]U0 = "<<UL<<"mV"<<endl;	
		*/
		
		/*
		for( int q = 0 ;q < range; q++){
			if(yR[q]<thrd&&flagR)
			//if(yR[q]<Rate*UR && flagR && yR[q]<thrd) 
			{
				xT0_R=x[q];
				flagR = 0;
				//cout<<"		[+] selected xR = "<<xT0_R<<"\t"<<yR[q]<<endl;
				}
			if(yL[q]<thrd&&flagL) 
			//if(yL[q]<Rate*UL && flagL && yL[q]<thrd) 
			{
				xT0_L=x[q];
				flagL = 0;
				//cout<<"		[+] selected xL = "<<xT0_L<<"\t"<<yL[q]<<endl;
				}
				//cout<<" q value"<<q<<endl;
		}
		*/
		
		
		//if(xT0_L[q]&&xT0_R[q]) 
		//cout<<"xT0_L[q]&&xT0_R[q] = "<<(xT0_L[q] && xT0_R[q])<<endl;
		xT0[q] = (xT0_L[q]+xT0_R[q])*0.5;
		/*cout<<"q="<<q<<"\t, xT0_L= "<<xT0_L[q]<<endl;
		cout<<"q="<<q<<"\t, xT0_R= "<<xT0_R[q]<<endl;
		cout<<"xT0_L[q]&&xT0_R[q] = "<<(xT0_L[q] && xT0_R[q])<<endl;
		cout<<"q="<<q<<"\t, xT0= "<<xT0[q]<<endl;*/
		
		}
		hSig[1]->Fill(xT0_L[0]);
		hSig[0]->Fill(xT0_R[0]);
		hSIG->Fill(xT0[0]);
		//cout<<"[-] Filled  xR:xL:x0 = "<<xT0_R<<"\t"<<xT0_L<<"\t"<<xT0<<endl;
		
		t2->Fill();
		

		
		//cout<<"loop k = "<<k<<endl;
		

	}
	/*
	cout<<xT0_L[0]<<"\t"<<xT0_L[1]<<endl;
	cout<<xT0_R[0]<<"\t"<<xT0_R[1]<<endl;
	cout<<xT0[0]<<"\t"<<xT0[1]<<endl;
	return;
	*/
	f1->Close();
	TGraph *gR = new TGraph(range,x,yR[0]);
	TGraph *gL = new TGraph(range,x,yL[0]);
	
	
	TCanvas *c = new TCanvas("c","",1600,600);
	
	
	//c->cd();
	gPad->Clear();
	c->Divide(2,1);
	c->cd(1);
	gR->Draw();
	TAxis *xA = gR->GetXaxis();
	TAxis *yA = gR->GetYaxis();
	gR->SetLineWidth(2);
	gR->SetLineColor(kRed);
	xA->SetTitle("Time (s)");
	xA->SetRangeUser(RL,RR);
	xA->SetRangeUser(5,10);
	yA->SetTitle("Amplitude (mV)");
	//myFun->SetNpx(5e3);
	//myFun->Draw();
	//gR->Draw("AC");
	
	gL->SetLineWidth(2);
	gL->SetLineColor(kBlue);
	gL->Draw("same");
	//gL->GetXaxis()->SetTitle("Time (s)");
	//gL->Draw("same");
	
	TLegend* leg; 
	leg = DrawMyLeg(0.6,0.2,0.75,0.32);
	leg->AddEntry(gR,"PMT_Right","lp");
	leg->AddEntry(gL,"PMT_Left","lp");
	leg->Draw();
	//gR->SetHistogram(hRSig);
	//gL->SetHistogram(hLSig);
	//hRSig->Draw();
	//hLSig->Draw("same");

	
	c->cd(2);
	h[0]->SetTitle("t0");
	h[0]->GetXaxis()->SetTitle("Time (s)");
	h[0]->GetYaxis()->SetTitle("Counts");
	h[0]->SetLineColor(kRed);
	h[0]->Draw();
	h[1]->SetLineColor(kBlue);
	h[1]->Draw("SAMES");
	
	gPad->Update();
	
	
	TPaveStats *stL = (TPaveStats*)h[1]->FindObject("stats");
	//stR->SetName("PMT_Right");
	stL->SetY1NDC(0.6);
	stL->SetY2NDC(0.78);
	gPad->Modified();
	gPad->Update();
	

	//t1->Draw("PmtR.t[0]>>ht(200,0,20)");
	
	//cout<<myFun->GetParameter(0)<<endl;

	//Float_t tRise = g->GetX(0.9*U0,RL,t0)-g->GetX(0.1*U0,RL,t0);
	//Float_t tFall = g->GetX(0.1*U0,t0,RR)-g>GetX(0.9*U0,t0,RR);
	//cout<<"The Rise time is "<< tRise*1e12 <<"ps"<<endl;
	//cout<<"The fall time is "<< tFall*1e12 <<"ps"<<endl;

	
	
	
	sprintf(buff,"%s_Signal.png",name);
	c->SaveAs(buff);
	//return ;

	//hSig[0]->FillRandom("gaus",1000);
	//hSig[0]->Draw();
	
	TCanvas *c1 = new TCanvas("c1","",1600,600);
	
	c1->Divide(2,1);
	c1->cd(1);
	hSig[0]->Draw();
	hSig[0]->Rebin(6);
	hSig[0]->GetXaxis()->SetRangeUser(t_L,t_R);
	
	hSig[0]->Fit(fSig[0],"R");
	c1->cd(2);
	hSig[1]->Draw();
	hSig[1]->Rebin(6);
	hSig[1]->GetXaxis()->SetRangeUser(t_L,t_R);
	
	hSig[1]->Fit(fSig[1],"R");
	sprintf(buff,"%s_Twosides_timeresolution.png",name);
	c1->SaveAs(buff);
	
	TCanvas *c2 = new TCanvas("c2","",800,600);
	c2->cd();
	hSIG->Draw();

	
	hSIG->Rebin(6);
	hSIG->GetXaxis()->SetRangeUser(t_L,t_R);
	hSIG->Fit(fSIG,"R");
	sprintf(buff,"%s_timeresolution.png",name);
	c2->SaveAs(buff);

	f2->cd();
	t2->Write();
	
	//f2->Close();
	//sprintf(buff,"%s_TimeRes.dat",name);
	ofstream outputdata("TimeRes.dat",ios::app);
	outputdata<<fSIG->GetParameter(2)<<"\t"<<fSIG->GetParError(2)<<endl;
	


    
	//}
	cout<<"The process is over,THANK YOU!"<<endl;

	//c->Delete();
	vector<double>().swap(*TR);
	vector<double>().swap(*TL);
	vector<int>().swap(*IDL);
	vector<int>().swap(*IDR);
	delete TR;
	delete TL;
	delete IDL;
	delete IDR;
	//delete count;
	
}

TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05){
  TLegend *leg = new TLegend(xlow,ylow,xup,yup);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(10);
  leg->SetTextFont(textFont);
  leg->SetTextSize(textSize);
  //leg->Draw("same");
  return leg;
}

TLatex* DrawMyLatex(char* text, Double_t x=0.65, Double_t y=0.5, Int_t textFont=62, Size_t textSize=0.05, Color_t colorIndex=2){
  TLatex *latex = new TLatex(x,y,text);
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(textSize);
  latex->SetTextColor(colorIndex);
  latex->Draw("same");
  return latex;
}

/*
void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1){
     datahist->SetLineColor( LColor );
     datahist->SetLineWidth( LWidth );
 
     datahist->GetXaxis()->SetTitle( xtitle);
     datahist->GetYaxis()->SetTitle( ytitle);
     datahist->GetXaxis()->SetAxisColor(1);
     datahist->GetYaxis()->SetAxisColor(1);
     datahist->GetXaxis()->SetLabelColor(1);
     datahist->GetYaxis()->SetLabelColor(1);
     datahist->GetXaxis()->SetLabelFont( 42 );
     datahist->GetYaxis()->SetLabelFont( 42 );
     datahist->GetXaxis()->SetLabelSize( 0.06 );
     datahist->GetYaxis()->SetLabelSize( 0.06 );
     datahist->GetXaxis()->SetLabelOffset( 0.01 );
     datahist->GetYaxis()->SetLabelOffset( 0.01 );
     datahist->GetXaxis()->SetTitleFont( 42 );
     datahist->GetYaxis()->SetTitleFont( 42 );
     datahist->GetXaxis()->SetTitleColor( TitleColor);
     datahist->GetYaxis()->SetTitleColor( TitleColor );
     datahist->GetXaxis()->SetTitleSize(0.06);
     datahist->GetYaxis()->SetTitleSize(0.06);
     datahist->GetXaxis()->SetTitleOffset(1.0);
     datahist->GetYaxis()->SetTitleOffset(1.0);
     //datahist->GetXaxis()->SetBorderSize(5);
     datahist->GetXaxis()->SetNdivisions(510);
     datahist->GetYaxis()->SetNdivisions(510);
     datahist->GetXaxis()->CenterTitle();
     datahist->GetYaxis()->CenterTitle();
 
}
*/








