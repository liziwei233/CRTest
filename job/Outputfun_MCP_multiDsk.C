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
TF1* pol3fit(TGraph* g,float U_RL, float U_RR);

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
    if(par.empty()) return 0;

    Double_t val = 0;
    Double_t tts = 0;
    double SPEpar[7];
    double Tmark=0;
    bool flag;
    double Trecept=0.7e-9; //waiting the photons hit
    double Treject=1.0e-9; //recover time,during this time any photons are rejected.
    //tts = r.Gaus(0,10.6e-12);
    //
    //----MCP R10754------
    //------------------------------
    //
    SPEpar[0]=1.1e6;  //Gain
    SPEpar[1]=1.6e-19; //e
    SPEpar[2]=2.65e-12;  //Ca  ??
    SPEpar[3]=12e-12; //C
    SPEpar[4]=10e3;    //R   ??
    SPEpar[5]=50;      //Z
    SPEpar[6]=90e-12;  //rise time
    //int N;
    //N=sizeof(par)/sizeof(par[0]);
    sort(par.begin(),par.end());
    int n=0;
    Tmark=par.at(0);
    for (int n=0;n<par.size();n++){
        //while(par[n]>5e-9){

        if (x-par.at(n)<-30.e-9){
            val+=0;
        }
        else{

            if(par.at(n)-Tmark<Trecept)
            {
                r.SetSeed(par.at(n));
                tts = r.Gaus(0,10.6e-12); //TTS of MCP-R3805U
                //cout<<"tts= "<<tts<<endl;
                val+=response(x-tts-par.at(n),SPEpar);
            }
            else if(par.at(n)-Tmark<(Trecept+Treject))
            {
                val+=0;
            }
            else 
            {
                Tmark=par.at(n);
                r.SetSeed(par.at(n));
                tts = r.Gaus(0,10.6e-12); //TTS of MCP-R3805U
                //cout<<"tts= "<<tts<<endl;
                val+=response(x-tts-par.at(n),SPEpar);
            }


        }

    }
    //cout<<"n = "<<n<<endl;

    return val;



    }
    TF1* pol3fit(TGraph* g,float U_RL, float U_RR){



        g->Draw();

        TF1 *fitU = new TF1("fitU","pol3",U_RL,U_RR);
        //fitU->SetParameter(1,g->GetMean());
        TFitResultPtr p3 = g->Fit(fitU,"R");
        //cout<<"p3=\t"<<p3<<endl;
        if(p3)	{
            TF1 *fit2 = new TF1("fit2","pol2",U_RL,U_RR);
            TFitResultPtr p2 = g->Fit(fit2,"R");
            //cout<<"p2=\t"<<p2<<endl;

            if(p2){
                TF1 *fit1 = new TF1("fit1","pol1",U_RL,U_RR);
                TFitResultPtr p1 = g->Fit(fit1,"R");
                //	cout<<"p1=\t"<<p1<<endl;
                return fit1;

            }
            return fit2;

        }
        //U_RL = fitU->GetParameter(1)-5*fitU->GetParameter(2);
        //U_RR = fitU->GetParameter(1)+5*fitU->GetParameter(2);

        //g->GetXaxis()->SetRangeUser(U_RL,U_RR);
        //g->Fit(fitU);
        return fitU;

    }



    void Outputfun_MCP_multiDsk(const char *rootname="",double fac = -30, const char* ParType="FIX"){

        cout<<"fac="<<fac<<",dicriminate:"<<ParType<<endl;

        /*===========================================
         * ============Procedure timing start========
         * =========================================*/
        clock_t start,finish;
        double totaltime;
        start=clock();

        TLegend* DrawMyLeg(Double_t xlow=0.2, Double_t ylow=0.2, Double_t xup=0.5, Double_t yup=0.5, Int_t textFont=62, Size_t textSize=0.05);

        TLatex* DrawMyLatex(char* text, Double_t x=0.65, Double_t y=0.5, Int_t textFont=62, Size_t textSize=0.05, Color_t colorIndex=2);
        void DrawMyGraph(TGraph *datagraph, char *xtitle, char *ytitle, Size_t MSize=1, Style_t MStyle =28, Color_t MColor=1, Color_t LColor=1, Width_t LWidth=1, Style_t LStyle=1, Color_t FColor=16);
        void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom);
        void DrawMyPad(TVirtualPad *pad,const char* xname,const char* yname,float x1,float x2,float y1,float y2);

        //void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);




        gStyle->SetOptStat(1);
        gStyle->SetOptFit(1111);
        gStyle->SetOptTitle(0);

        //TRandom3 r;
        //r.SetSeed(0);
        char name[1024];
        char buff[1024];


        //Double_t parR[500]={};
        //Double_t parL[500]={};
        vector<double> parR;
        vector<double> parL;

        //Double_t RL = -5e-9;
        //Double_t RR = 20e-9;
        Double_t RL = -0e-9;
        Double_t RR = 10e-9;
        Double_t zoomRL = -0e-9;
        Double_t zoomRR = 10e-9;
        int binNum=0;
        binNum = (RR-RL)/25e-12;

        const int range =1e3;  // 25ps/sample
        double Rate=0;

        bool flagR=0,flagL=0;
        double xT0_L=0,xT0_R=0,xT0=0;
        int indexL=0,indexR=0;
        double keypointL=0,keypointR=0;
        double UL=0,UR=0;
        const int certain=2;	

        vector<double>* TR;
        vector<double>* TL;
        TL = new vector<double>;
        TR = new vector<double>;
        //count = new vector<int>;
        int N=0,temp=0;
        Double_t xR[range]={};
        Double_t xL[range]={};
        Double_t yR[range]={};
        Double_t yL[range]={};

        double Uspe = -30; //Umax = -28.94mV
        double ratio[5]={0.1,0.15,0.2,0.25,0.3};
        double thrd[5]={1,2,3,4,5};
        double xFIX[2][5]={};
        double xCFD[2][5]={};	
        double amp[2]={};
        int index_pkL=0,index_pkR=0;
        int idxCFD[2][5]={};
        int idxFIX[2][5]={};


        //sprintf(name,"%s_fac%g_type%s",rootname,fac,ParType);
        sprintf(name,"%s",rootname);
        sprintf(buff,"%s.root",rootname);

        TFile *f1 = new TFile(buff,"READ");
        TTree *t1 = (TTree*)f1->Get("Run");

        t1->SetBranchAddress("PmtR.t",&TR);
        t1->SetBranchAddress("PmtL.t",&TL);

        //sprintf(name,"Thrd_%g",abs(thrd));	

        sprintf(buff,"%sdata.root",name);

        TFile *f2 = new TFile(buff,"RECREATE");
        TTree *t2 = new TTree("dsk","restore analysed data  from G4");

        t2->Branch("amp",amp,"amp[2]/D");	
        t2->Branch("T0_CFD",xCFD,"xCFD[2][5]/D");	
        t2->Branch("CFD_fac",ratio,"ratio[5]/D");
        t2->Branch("T0_FIX",xFIX,"xFIX[2][5]/D");	
        t2->Branch("FIX_fac",thrd,"thrd[5]/D");


        //for(int s = 0; s<4;s++){



        double t_L = RL;
        double t_R = RR;


        //f1->cd();


        //TF1 *myFun;
        TH1D *h[2];
        h[0] = new TH1D("hR","",binNum,RL,RR);
        h[1] = new TH1D("hL","",binNum,RL,RR);

        N = t1->GetEntries();
        cout<<"Entries = "<<N<<endl;


        //count->clear();
        //for(int i = certain; i < certain+1; i++){
        for(int i = 0; i < N; i++){

            //-----------initial----------------------//
            TL->clear();
            TR->clear();

            h[0]->Reset();
            h[1]->Reset();

            //parR.clear();
            //parL.clear();
            vector<double>().swap(parR);
            vector<double>().swap(parL);

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
                parR.push_back((*TR)[k]*1e-9);
                //parR[k]=8.3e-9;
                //cout<<" [+] par "<<k<<"\t"<<parR.at(k)<<endl;
                //cout<<"par"<<k<<" = "<<par[k]<<endl;

                h[0]->Fill(parR.at(k));

                //myFun->SetParameter(k,par[k]);
            }

            temp = TL->size();	
            //cout<<"counterL = "<< temp <<endl;
            for(int k=0;k<temp;k++){
                //cout<< T[][k] <<endl;
                parL.push_back((*TL)[k]*1e-9);
                //cout<<" [+] par "<<k<<"\t"<<parL.at(k)<<endl;


                h[1]->Fill(parL.at(k));

                //myFun->SetParameter(k,par[k]);
            }
            flagR = 1;
            flagL = 1;

            //cout<<"hello"<<endl;
            //cout<<"parL.size() = "<<parL.size()<<endl;
            //cout<<"parR.size() = "<<parR.size()<<endl;

            for(int j=0;j<range;j++){
                xR[j]=(RR-RL)/range*j+RL;
                yR[j]=outputfunc(xR[j],parR);
                //cout<<"process check======>"<<endl;
                xL[j]=(RR-RL)/range*j+RL;
                yL[j]=outputfunc(xL[j],parL);



                //if(yR[j]<thrd&&flagR) {xT0_R=x[j];flagR = false;}
                //if(yL[j]<thrd&&flagL) {xT0_L=x[j];flagL = false;}
                //cout<<"[+] x"<<j<<":y"<<j<<"=\t"<<x[j]<<"\t"<<yR[j]<<"\t"<<yL[j]<<endl;
            }



            /*================================ 
             *=======ZOOM OUT the leading egde;
             *=================================*/
            zoomRR = xR[TMath::LocMin(range,yR)]+1e-9;
            zoomRL = xR[TMath::LocMin(range,yR)]-1e-9;
            //cout<<zoomRL<<"\t"<<zoomRR<<endl;
            //return;
            for(int j=0;j<range;j++){
                xR[j]=((zoomRR-zoomRL)/range*j+zoomRL)*1;
                yR[j]=outputfunc(xR[j],parR);
                xR[j]=((zoomRR-zoomRL)/range*j+zoomRL)*1e9;

            }

            zoomRR = xL[TMath::LocMin(range,yL)]+1e-9;
            zoomRL = xL[TMath::LocMin(range,yL)]-1e-9;
            for(int j=0;j<range;j++){
                xL[j]=((zoomRR-zoomRL)/range*j+zoomRL)*1;
                yL[j]=outputfunc(xL[j],parL);
                xL[j]=((zoomRR-zoomRL)/range*j+zoomRL)*1e9;

            }
            UR = TMath::MinElement(range,yR);
            UL = TMath::MinElement(range,yL);
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


            /*=================================
             *=======Discriminate the signal===
             *=================================*/
            xT0_R = 0;
            xT0_L = 0 ;
            xT0 = 0;
            keypointL = 0;
            keypointR = 0;
            index_pkL = TMath::LocMin(range,yL);
            index_pkR = TMath::LocMin(range,yR);
            amp[0]=UL;
            amp[1]=UR;
            for( int q = 1; q < range; q++)
            {
                for( int s = 0; s < 5; s++)
                {
                    //Left PMT discriminator

                    if(yL[q]<ratio[s]*UL&&yL[q-1]>ratio[s]*UL&&q<=index_pkL)
                    {
                        xCFD[0][s]=xL[q];
                        idxCFD[0][s]=q;
                    }

                    if(yL[q]<thrd[s]*Uspe&&yL[q-1]>thrd[s]*Uspe&&q<=index_pkL)
                    {

                        xFIX[0][s]=xL[q];
                        idxFIX[0][s]=q;
                    }

                    //Right PMT discrimanator

                    if(yR[q]<ratio[s]*UR&&yR[q-1]>ratio[s]*UR&&q<=index_pkR)
                    {

                        xCFD[1][s]=xR[q];
                        idxCFD[1][s]=q;
                    }
                    if(yR[q]<thrd[s]*Uspe&&yR[q-1]>thrd[s]*Uspe&&q<=index_pkR)
                    {

                        xFIX[1][s]=xR[q];
                        idxFIX[1][s]=q;
                    }



                }
            }
            /*
               for (int s = 0;s < 5; s++){
               for(
               int p = 0;p < 2;p++)
               {
               cout<<"xCFD_"<<p<<"="<<xCFD[p][s]<<endl;   
               cout<<"xFIX_"<<p<<"="<<xFIX[p][s]<<endl;   
               }
               }
               return;
               */
            /*
               if (strcmp(ParType,"FIX")==0) //if the discriminate way is fix threshold discrim
               {
               keypointL = fac;
               keypointR = fac;
               }
               else {
               keypointR = fac*UR;

               keypointL = fac*UL;
               }


               for( int q = 0 ;q < range; q++){
               if(yR[q]<keypointR&&flagR)
            //if(yR[q]<Rate*UR && flagR && yR[q]<thrd) 
            {

            indexR=q;
            //xT0_R=xR[q];
            flagR = 0;
            //cout<<"		[+] selected xR = "<<xT0_R<<"\t"<<yR[q]<<endl;
            }
            if(yL[q]<keypointL&&flagL) 
            //if(yL[q]<Rate*UL && flagL && yL[q]<thrd) 
            {

            indexL=q;
            //xT0_L=xL[q];
            flagL = 0;
            //cout<<"		[+] selected xL = "<<xT0_L<<"\t"<<yL[q]<<endl;
            }
            //cout<<" q value"<<q<<endl;
            }
            //cout<<"find the time stamp (ns) = "<<xR[indexR]<<",and the corrresponding amp = "<<yR[indexR]<<endl;
            //cout<<"index="<<indexR<<endl;
            */
            /*=====================================================
             * =======Fit the signal and find the timestamp========
             * ==================================================*/
            TGraph *gR = new TGraph(range,xR,yR);
            TGraph *gL = new TGraph(range,xL,yL);
            for(int s=0;s<5;s++)
            {

                // PMT Left
                TF1* fitCFDL = pol3fit(gL,xL[idxCFD[0][s]]-2e-2,xL[idxCFD[0][s]]+3e-2);
                xCFD[0][s]=fitCFDL->GetX(ratio[s]*UL)*1e-9;
                TF1* fitFIXL = pol3fit(gL,xL[idxFIX[0][s]]-2e-2,xL[idxFIX[0][s]]+3e-2);
                xFIX[0][s]=fitFIXL->GetX(thrd[s]*Uspe)*1e-9;


                //PMT right
                TF1* fitCFDR = pol3fit(gR,xR[idxCFD[1][s]]-2e-2,xR[idxCFD[1][s]]+3e-2);
                xCFD[1][s]=fitCFDR->GetX(ratio[s]*UR)*1e-9;


                TF1* fitFIXR = pol3fit(gR,xR[idxFIX[1][s]]-2e-2,xR[idxFIX[1][s]]+3e-2);
                xFIX[1][s]=fitFIXR->GetX(thrd[s]*UR)*1e-9;
            }

            //return;
            //xT0_L = Discriminate(xL,yL,indexL);

            //cout<<"[-] Event No. Filled  xR:xL:x0 = "<<i<<"\t"<<xR[indexR]<<"\t"<<xL[indexL]<<"\t"<<(xR[indexR]+xL[indexL])/2<<endl;
            //cout<<"[-] Event No. Filled  xR:xL:x0 = "<<i<<"\t"<<xT0_R<<"\t"<<xT0_L<<"\t"<<xT0<<endl;

            t2->Fill();



            //cout<<"loop k = "<<k<<endl;


        }

        f1->Close();
        TGraph *gR = new TGraph(range,xR,yR);
        TGraph *gL = new TGraph(range,xL,yL);


        TCanvas *c = new TCanvas("c","",1600,600);


        //c->cd();
        gPad->Clear();
        c->Divide(2,1);
        c->cd(1);
        SetMyPad(gPad,0.12,0.05,0.05,0.12);
        gR->Draw("AL");
        gPad->Update();

        float ymin = 0;
        float ymax = 0;
        ymax = gPad->GetUymax();
        if(UR<UL) ymin = UR*1.2;
        else ymin = UL*1.2;
        DrawMyPad(gPad,"Time (ns)","Amp (mV)",RL*1e9, RR*1e9,ymin, ymax);

        //return;
        //gR->Draw();
        DrawMyGraph(gR,"Time (ns)","Amp (mV)",1,28,1,kRed,2,1);
        gR->Draw("L");
        //return;

        //xA->SetRangeUser(5,10);
        //yA->SetTitle("Amplitude (mV)");
        //myFun->SetNpx(5e3);
        //myFun->Draw();
        //gR->Draw("AC");

        DrawMyGraph(gL,"Time (ns)","Amp (mV)",1,28,1,kBlue,2,1);
        gL->Draw("L");
        //gL->GetXaxis()->SetTitle("Time (s)");
        //gL->Draw("same");

        TLegend* leg; 
        leg = DrawMyLeg(0.6,0.2,0.75,0.32);
        leg->AddEntry(gR,"PMT_Right","l");
        leg->AddEntry(gL,"PMT_Left","l");
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
        f2->cd();
        t2->Write();

        //f2->Close();
        //sprintf(buff,"%s_TimeRes.dat",name);




        //}
        cout<<"The process is over,THANK YOU!"<<endl;

        //c->Delete();
        vector<double>().swap(*TR);
        vector<double>().swap(*TL);
        delete TR;
        delete TL;

        /*=======================================================*
         * ================Procedure timing end==================*
         * ======================================================*/
        finish=clock();
        totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
        cout<<"\nthe whole time through the procedure is "<<totaltime<<"s!!"<<endl;//delete count;

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

    void DrawMyPad(TVirtualPad *pad,const char* xname,const char* yname,float x1,float x2,float y1,float y2){


        TH1F* hpad = pad->DrawFrame(x1,y1,x2,y2);
        hpad->GetXaxis()->SetTitle( xname);
        hpad->GetYaxis()->SetTitle( yname);
        hpad->GetXaxis()->SetAxisColor(1);
        hpad->GetYaxis()->SetAxisColor(1);
        hpad->GetXaxis()->SetLabelColor(1);
        hpad->GetYaxis()->SetLabelColor(1);
        hpad->GetXaxis()->SetLabelFont( 42 );
        hpad->GetYaxis()->SetLabelFont( 42 );
        hpad->GetXaxis()->SetLabelSize( 0.05 );
        hpad->GetYaxis()->SetLabelSize( 0.05 );
        hpad->GetXaxis()->SetLabelOffset( 0.01 );
        hpad->GetYaxis()->SetLabelOffset( 0.01 );
        hpad->GetXaxis()->SetTitleFont( 42 );
        hpad->GetYaxis()->SetTitleFont( 42 );
        //hpad->GetXaxis()->SetTitleColor( TitleColor);
        //hpad->GetYaxis()->SetTitleColor( TitleColor );
        hpad->GetXaxis()->SetTitleSize(0.06);
        hpad->GetYaxis()->SetTitleSize(0.06);
        hpad->GetXaxis()->SetTitleOffset(1.0);
        hpad->GetYaxis()->SetTitleOffset(1.0);
        pad->Modified();
        pad->Update();

    }


    void SetMyPad(TVirtualPad *pad,float left, float right, float top, float bottom){
        pad->SetFillColor(10);
        pad->SetBorderMode(0);
        pad->SetBorderSize(0);
        pad->SetFrameFillColor(10);
        //pad->SetFrameFillStyle(3003);
        pad->SetFrameBorderMode(0);
        pad->SetFrameBorderSize(0);
        pad->SetLeftMargin(left);
        pad->SetRightMargin(right);
        pad->SetTopMargin(top);
        pad->SetBottomMargin(bottom);
    }



    void DrawMyGraph(TGraph *datagraph, char *xtitle, char *ytitle, Size_t MSize=1, Style_t MStyle =28, Color_t MColor=1, Color_t LColor=1, Width_t LWidth=1, Style_t LStyle=1, Color_t FColor=16){
        datagraph->SetLineColor( LColor );
        datagraph->SetLineWidth( LWidth );
        datagraph->SetLineStyle( LStyle );
        datagraph->SetMarkerSize( MSize );
        datagraph->SetMarkerStyle( MStyle );
        datagraph->SetMarkerColor( MColor );
        datagraph->SetFillColor( FColor );
        //datagraph->SetFillStyle( FStyle );
        datagraph->GetXaxis()->SetTitle( xtitle);
        datagraph->GetYaxis()->SetTitle( ytitle);
        datagraph->GetXaxis()->SetAxisColor(1);
        datagraph->GetYaxis()->SetAxisColor(1);
        datagraph->GetXaxis()->SetLabelColor(1);
        datagraph->GetYaxis()->SetLabelColor(1);
        datagraph->GetXaxis()->SetLabelFont( 42 );
        datagraph->GetYaxis()->SetLabelFont( 42 );
        datagraph->GetXaxis()->SetLabelSize( 0.05 );
        datagraph->GetYaxis()->SetLabelSize( 0.05 );
        datagraph->GetXaxis()->SetLabelOffset( 0.01 );
        datagraph->GetYaxis()->SetLabelOffset( 0.01 );
        datagraph->GetXaxis()->SetTitleFont( 42 );
        datagraph->GetYaxis()->SetTitleFont( 42 );
        //datagraph->GetXaxis()->SetTitleColor( TitleColor);
        //datagraph->GetYaxis()->SetTitleColor( TitleColor );
        datagraph->GetXaxis()->SetTitleSize(0.06);
        datagraph->GetYaxis()->SetTitleSize(0.06);
        datagraph->GetXaxis()->SetTitleOffset(1.0);
        datagraph->GetYaxis()->SetTitleOffset(1.0);
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








