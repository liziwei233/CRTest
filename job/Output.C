#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

TRandom3 r;
TH1D *hr =new TH1D("hr","hr",2000,-0.1e-9,0.1e-9);
using namespace std;

Double_t outputfunc(Double_t x, vector<double> par, vector<double> tts);
Double_t response(Double_t x, Double_t par[7]);
TF1 *pol3fit(TGraph *g, float U_RL, float U_RR);

Double_t response(Double_t x, Double_t par[7])
{
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
    Double_t a, b, beta, gamma;

    // how to accumulate these single photon signal?
    // which distribution does t0 sample?
    //!! get t0 by result of simulation!!

    beta = ((R + Z) * C + 2.0 * Ca * R) / (2.0 * C * Ca * R * Z);
    gamma = 1.0 / (2.0 * C * Ca * R * Z);
    a = -(beta + TMath::Sqrt(beta * beta - 4.0 * gamma)) / 2.0;
    b = -(beta - TMath::Sqrt(beta * beta - 4.0 * gamma)) / 2.0;

    val = -0.5 * G * Q * (b * TMath::Exp(b * x + 0.5 * b * b * rise * rise) * TMath::Erfc((-b * rise - x / rise) / TMath::Sqrt(2.0)) - a * TMath::Exp(a * x + 0.5 * a * a * rise * rise) * TMath::Erfc((-a * rise - x / rise) / TMath::Sqrt(2.0))) / (Ca * (b - a)) * 1e3;

    return val;
}

Double_t outputfunc(Double_t x, vector<double> par, vector<double> tts)
{
    if (par.empty())
        return 0;

    Double_t val = 0;
    //Double_t tts = 0;
    double SPEpar[7];
    double Tmark = 0;
    bool flag;
    double Trecept = 10e-12;  //waiting the photons hit
    double Treject = 500e-12; //recover time,during this time any photons are rejected.
    //double ttssigma = 20e-12;
    //tts = r.Gaus(0,10.6e-12);
    //
    //----MCP R10754------
    //------------------------------
    //
    SPEpar[0] = 1.1e6;    //Gain
    SPEpar[1] = 1.6e-19;  //e
    SPEpar[2] = 2.65e-12; //Ca  ??
    SPEpar[3] = 12e-12;   //C
    SPEpar[4] = 10e3;     //R   ??
    SPEpar[5] = 50;       //Z
    SPEpar[6] = 90e-12;   //rise time
    //int N;
    //N=sizeof(par)/sizeof(par[0]);
    sort(par.begin(), par.end()); // in order from small to large.
    int n = 0;
    Tmark = par.at(0);
    for (int n = 0; n < par.size(); n++)
    {
        //while(par[n]>5e-9){

        if (x - par.at(n) < -30.e-9)
        {
            val += 0;
        }
        else
        {

            if (par.at(n) - Tmark < Trecept)
            {
                //r.SetSeed(n);
                //r.SetSeed(0);
                //tts = r.Gaus(0, ttssigma); //TTS of MCP-R3805U
                //hr->Fill(tts);
                //cout<<"tts= "<<tts<<endl;
                val += response(x - tts.at(n) - par.at(n), SPEpar);
            }
            else if (par.at(n) - Tmark < (Trecept + Treject))
            {
                val += 0;
            }
            else
            {
                Tmark = par.at(n);
                //r.SetSeed(0);
                //tts = r.Gaus(0, ttssigma); //TTS of MCP-R3805U
                //hr->Fill(tts);
                //cout<<"tts= "<<tts<<endl;
                val += response(x - tts.at(n) - par.at(n), SPEpar);
            }
        }
    }
    //cout<<"n = "<<n<<endl;

    return val;
}
TF1 *pol3fit(TGraph *g, float U_RL, float U_RR)
{

    g->Draw();

    TF1 *fitU = new TF1("fitU", "pol3", U_RL, U_RR);
    //fitU->SetParameter(1,g->GetMean());
    TFitResultPtr p3 = g->Fit(fitU, "R");
    //cout<<"p3=\t"<<p3<<endl;
    if (p3)
    {
        TF1 *fit2 = new TF1("fit2", "pol2", U_RL, U_RR);
        TFitResultPtr p2 = g->Fit(fit2, "R");
        //cout<<"p2=\t"<<p2<<endl;

        if (p2)
        {
            TF1 *fit1 = new TF1("fit1", "pol1", U_RL, U_RR);
            TFitResultPtr p1 = g->Fit(fit1, "R");
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

void Output(const char *rootname = "", double fac = 0.2, const char *ParType = "CFD",unsigned long processN=1)
{
    //void Outputfun_MCP(const char *rootname="",double fac = -30, const char* ParType="FIX"){

    cout << "fac=" << fac << ",dicriminate:" << ParType << endl;

    /*===========================================
         * ============Procedure timing start========
         * =========================================*/
    clock_t start, finish;
    double totaltime;
    start = clock();

    TLegend *DrawMyLeg(Double_t xlow = 0.2, Double_t ylow = 0.2, Double_t xup = 0.5, Double_t yup = 0.5, Int_t textFont = 62, Size_t textSize = 0.05);

    TLatex *DrawMyLatex(char *text, Double_t x = 0.65, Double_t y = 0.5, Int_t textFont = 62, Size_t textSize = 0.05, Color_t colorIndex = 2);
    
    void DrawMyGraph(TGraph * datagraph, char *xtitle, char *ytitle, Size_t MSize = 1, Style_t MStyle = 28, Color_t MColor = 1, Color_t LColor = 1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t FColor = 16);
    
    void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
    
    void SetMyPad(TVirtualPad * pad, float left, float right, float top, float bottom);

    void DrawMyPad(TVirtualPad * pad, const char *xname, const char *yname, float x1, float x2, float y1, float y2);

    //void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    //TRandom3 r;
    //r.SetSeed(0);
    char name[1024];
    char buff[1024];

    const int T = 2; //the number of copyvolume

    //Double_t parR[500]={};
    //Double_t parL[500]={};

    //Double_t RL = -5e-9;
    //Double_t RR = 20e-9;
    Double_t RL = -2e-9;
    Double_t RR = 8e-9;
    Double_t zoomRL = -2e-9;
    Double_t zoomRR = 8e-9;
    int binNum = 0;
    binNum = (RR - RL) / 5e-12;

    double ttssigma=20e-12;
    const int range = 400; // 25ps/sample
    Double_t thrd = -30;   //Umax = -28.94mV
    double Rate = 0;

    bool flag = 0;
    int index = 0;
    double keypoint = 0;
    double xT0[T] = {0};
    double U[T] = {0};
    double InX[T]={0};
    double InY[T]={0};
    double T0_1stpe[T]={0};
    const int certain = 2;

	vector<vector<double> > par(T);
	vector<vector<double> > tts(T);

    vector<double> *hitT;
    vector<double> *IncidX; //the position of incident event
    vector<double> *IncidY;
    vector<int> *ID;

    hitT = new vector<double>;
    IncidX = new vector<double>;
    IncidY = new vector<double>;
    ID = new vector<int>;

    //count = new vector<int>;
    int N = 0, hitN=0;
    double temp = 0.;
    Double_t x[range] = {};
    Double_t y[T][range] = {};

    sprintf(name, "%s", rootname);
    sprintf(buff, "%s.root", name);

    TFile *f1 = new TFile(buff, "READ");
    TTree *t1 = (TTree *)f1->Get("Run");

    t1->SetBranchAddress("PmtS.t", &hitT);
    t1->SetBranchAddress("PmtS.id", &ID);
    t1->SetBranchAddress("ph.x", &IncidX);
    t1->SetBranchAddress("ph.y", &IncidY);

    //sprintf(name,"Thrd_%g",abs(thrd));

    sprintf(buff, "%sdata.root", name);

    TFile *f2 = new TFile(buff, "RECREATE");
    TTree *t2 = new TTree("data", "restore analysed data  from G4");
    sprintf(buff, "U[%d]/D", T);
    t2->Branch("U", &U, buff);       //-mV
    sprintf(buff, "xT0[%d]/D", T);  // ns
    t2->Branch("xT0", &xT0, buff);
    sprintf(buff, "T0_1stpe[%d]/D", T); //ns
    t2->Branch("T0_1stpe", &T0_1stpe, buff);
    sprintf(buff, "InX[%d]/D", T);
    t2->Branch("InX", &InX, buff);  //mm
    sprintf(buff, "InY[%d]/D", T);
    t2->Branch("InY", &InY, buff);  //mm
    //for(int s = 0; s<4;s++){

    double t_L = -2;
    double t_R = 2;

    //f1->cd();

    //TF1 *myFun;
    TH1D *h[T];
    TH1D *hSig[T];
    TF1 *fSig[T];
    TGraph *g[T];
    for (int iT = 0; iT < T; iT++)
    {
        sprintf(buff, "h%d", iT);
        h[iT] = new TH1D(buff, buff, binNum, RL*1e9, RR*1e9);
        sprintf(buff, "hSig%d", iT);
        hSig[iT] = new TH1D(buff, buff, binNum, t_L, t_R);
        sprintf(buff, "fSig%d", iT);
        fSig[iT] = new TF1(buff, "gaus", t_L, t_R);
    }

    N = t1->GetEntries();
    //N = 1;
    cout << "Entries = " << N << endl;

    //count->clear();
    //for(int i = certain; i < certain+1; i++){
    for (int i = 0; i < N; i++)
    {

        //-----------initial----------------------//
        hitT->clear();
        ID->clear();
        IncidX->clear();
        IncidY->clear();

        //parR.clear();
        //parL.clear();
        vector<vector<double> >(T).swap(par);
        vector<vector<double> >(T).swap(tts);

        //memset(parL,0,sizeof(parL));
        //memset(parR,0,sizeof(parR));

        //for(int i = certain; i < certain+1; i++){
        //par[i]=r.Gaus(2.4e-9,0.5e-9);
        //par[i]=4e-9;
        t1->GetEntry(i);
        hitN = hitT->size();
        //cout<<"counter = "<< hitN <<endl;
        //myFun = new TF1("myFun",outputfunc,RL,RR,temp);

        for (int iT = 0; iT < T; iT++)
        {
            InX[iT] = (*IncidX)[0];
            InY[iT] = (*IncidY)[0];
            for (int k = 0; k < hitN; k++)
            {
                //cout<< T[][k] <<endl;
                temp = (*hitT)[k] * 1e-9;
                if ((*ID)[k] == iT){
                    //r.SetSeed(time(NULL)*processN+k);
                    par[iT].push_back(temp);
                    tts[iT].push_back(r.Gaus(0,ttssigma));
                    if(i==N-1) h[iT]->Fill(temp*1e9);
                }
            }
            //parR[k]=8.3e-9;
            //cout<<" [+] par "<<k<<"\t"<<parR.at(k)<<endl;
            //cout<<"par"<<k<<" = "<<par[k]<<endl;
            sort(par[iT].begin(), par[iT].end());
        //cout<<">>> progress check <<<"<<endl;
           if(!par[iT].empty()) T0_1stpe[iT] = par[iT].at(0)*1e9;
            //myFun->SetParameter(k,par[k]);
        }

        //cout<<"hello"<<endl;
        //cout<<"parL.size() = "<<parL.size()<<endl;
        //cout<<"parR.size() = "<<parR.size()<<endl;
        for (int iT = 0; iT < T; iT++)
        {
            // Initial these variable
            memset(x, 0, sizeof(x));
            memset(y[iT], 0, sizeof(y[iT]));

            for (int j = 0; j < range; j++)
            {
                x[j] = (RR - RL) / range * j + RL;
                //cout<<"process check======>"<<endl;
                y[iT][j] = outputfunc(x[j], par[iT],tts[iT]);
                x[j] = ((RR - RL) / range * j + RL) * 1e9;

                //if(yR[j]<thrd&&flagR) {xT0_R=x[j];flagR = false;}
                //if(yL[j]<thrd&&flagL) {xT0_L=x[j];flagL = false;}
                //cout<<"[+] x"<<j<<":y"<<j<<"=\t"<<x[j]<<"\t"<<yR[j]<<"\t"<<yL[j]<<endl;
            }

            /*================================ 
             *=======ZOOM OUT the leading egde;
             *=================================
            zoomRR = x[TMath::LocMin(range, y[iT])] + 1e-9;
            zoomRL = x[TMath::LocMin(range, y[iT])] - 1e-9;
            //cout<<zoomRL<<"\t"<<zoomRR<<endl;
            //return;
            for (int j = 0; j < range; j++)
            {
                x[j] = ((zoomRR - zoomRL) / range * j + zoomRL) * 1;
                y[iT][j] = outputfunc(x[j], par[iT], tts[iT]);
                x[j] = ((zoomRR - zoomRL) / range * j + zoomRL) * 1e9;
            }
             *================================ 
             *=======ZOOM OUT the leading egde;
             *=================================*/

            U[iT] = TMath::MinElement(range, y[iT]);
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
            //Initial these variables
            flag = 1;
            xT0[iT] = 0;
            index = 0;
            keypoint = 0;
            if (strcmp(ParType, "FIX") == 0) //if the discriminate way is fix threshold discrim
            {
                keypoint = fac;
            }
            else
            {
                keypoint = fac * U[iT];
            }

            for (int q = 0; q < range; q++)
            {
                if (y[iT][q] < keypoint && flag)
                //if(yR[q]<Rate*UR && flagR && yR[q]<thrd)
                {
                    index = q;
                    //xT0_R=xR[q];
                    flag = 0;
                    //cout<<"		[+] selected xR = "<<xT0_R<<"\t"<<yR[q]<<endl;
                }
                //cout<<" q value"<<q<<endl;
            }
            //cout<<"find the time stamp (ns) = "<<xR[indexR]<<",and the corrresponding amp = "<<yR[indexR]<<endl;
            //cout<<"index="<<indexR<<endl;
            xT0[iT] = x[index];

            /*=====================================================
             * =======Fit the signal and find the timestamp========
             * ==================================================*/
            g[iT] = new TGraph(range, x, y[iT]);
            TF1 *fit = pol3fit(g[iT], x[index] - 60e-3, x[index] + 80e-3);
            xT0[iT] = fit->GetX(keypoint);
            /*=====================================================
             * =======Fit the signal and find the timestamp========
             * ==================================================*/


            //return;
            //xT0_L = Discriminate(xL,yL,indexL);

            hSig[iT]->Fill(xT0[iT]);
            //cout<<"[-] Event No. Filled  xR:xL:x0 = "<<i<<"\t"<<xR[indexR]<<"\t"<<xL[indexL]<<"\t"<<(xR[indexR]+xL[indexL])/2<<endl;
            //cout<<"[-] Event No. Filled  xR:xL:x0 = "<<i<<"\t"<<xT0_R<<"\t"<<xT0_L<<"\t"<<xT0<<endl;


            //cout<<"loop k = "<<k<<endl;
        }
            t2->Fill();
    }
    f1->Close();

    TCanvas *c = new TCanvas("c", "", 1600, 600);

    //c->cd();
    gPad->Clear();
    c->Divide(2, 1);
    c->cd(1);
    SetMyPad(gPad, 0.12, 0.05, 0.05, 0.12);
    gPad->Update();

    float ymin = 0;
    float ymax = 0;
    //ymax = gPad->GetUymax();
    ymin = TMath::MinElement(T,U) * 1.2;
    ymax = -0.2*ymin;
    DrawMyPad(gPad, "Time (ns)", "Amp (mV)", RL * 1e9, RR * 1e9, ymin, ymax);

    TLegend *leg;
    leg = DrawMyLeg(0.6, 0.2, 0.75, 0.32);
    //return;
    for(int iT;iT<T;iT++){
    DrawMyGraph(g[iT], "Time (ns)", "Amp (mV)", 1, 28, 1, iT+1, 2, 1);
    g[iT]->Draw("L");
    sprintf(buff,"PMT%d",iT);
    leg->AddEntry(g[iT], buff, "l");
    }
    leg->Draw();
    //return;

    c->cd(2);
    SetMyPad(gPad, 0.12, 0.05, 0.05, 0.12);
    TLatex* l[T];
    for(int iT=0;iT<T;iT++){

        DrawMyHist1(h[iT],"Time (ns)","Counts",iT+1,2);
        if (iT==0) h[iT]->Draw();
        else h[iT]->Draw("SAMES");
        sprintf(buff,"pmt%d hitN=%g",iT,h[iT]->GetSum());
        l[iT]=DrawMyLatex(buff,0.65,0.6-0.1*iT,62,0.05,iT+1);
    }

    gPad->Modified();
    gPad->Update();

    //t1->Draw("PmtR.t[0]>>ht(200,0,20)");

    //cout<<myFun->GetParameter(0)<<endl;

    //Float_t tRise = g->GetX(0.9*U0,RL,t0)-g->GetX(0.1*U0,RL,t0);
    //Float_t tFall = g->GetX(0.1*U0,t0,RR)-g>GetX(0.9*U0,t0,RR);
    //cout<<"The Rise time is "<< tRise*1e12 <<"ps"<<endl;
    //cout<<"The fall time is "<< tFall*1e12 <<"ps"<<endl;

    sprintf(buff, "%s_Signal.png", name);
    c->SaveAs(buff);
    //return ;

    //hSig[0]->FillRandom("gaus",1000);
    //hSig[0]->Draw();

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    SetMyPad(gPad, 0.12, 0.05, 0.05, 0.12);

    c1->cd();
    for(int iT=0;iT<T;iT++){

        DrawMyHist1(hSig[iT],"Time (ns)","Counts",iT+1,2);
        if (iT==0) hSig[iT]->Draw();
        else hSig[iT]->Draw("SAMES");
        hSig[iT]->Rebin(1);
        hSig[iT]->Fit(fSig[iT], "R");
        
        sprintf(buff,"pmt%d TR=%0.2f",iT,fSig[iT]->GetParameter(2));
        l[iT]=DrawMyLatex(buff,0.65,0.6-0.1*iT,62,0.05,iT+1);
    }


    sprintf(buff, "%s_Twosides_timeresolution.png", name);
    c1->SaveAs(buff);

/*
    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    c2->cd();
    hSIG->Draw();

    hSIG->Rebin(6);
    hSIG->GetXaxis()->SetRangeUser(t_L, t_R);
    hSIG->Fit(fSIG, "R");
    sprintf(buff, "%s_timeresolution.png", name);
    c2->SaveAs(buff);
*/
    f2->cd();
    t2->Write();

    //f2->Close();
    //sprintf(buff,"%s_TimeRes.dat",name);
    //ofstream outputdata("TimeRes.dat", ios::app);
    //outputdata << fSIG->GetParameter(2) << "\t" << fSIG->GetParError(2) << endl;

    //}
    //TCanvas *c3 = new TCanvas("c3","c3",800,600);
    //c3->cd();
    //hr->Draw();
    cout << "The process is over,THANK YOU!" << endl;

    //c->Delete();
    vector<double>().swap(*hitT);
    vector<double>().swap(*IncidX);
    vector<double>().swap(*IncidY);
    vector<int>().swap(*ID);
    delete hitT;
    delete IncidX;
    delete IncidY;
    delete ID;

    /*=======================================================*
         * ================Procedure timing end==================*
         * ======================================================*/
    finish = clock();
    totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << "\nthe whole time through the procedure is " << totaltime << "s!!" << endl; //delete count;
}

TLegend *DrawMyLeg(Double_t xlow = 0.2, Double_t ylow = 0.2, Double_t xup = 0.5, Double_t yup = 0.5, Int_t textFont = 62, Size_t textSize = 0.05)
{
    TLegend *leg = new TLegend(xlow, ylow, xup, yup);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextFont(textFont);
    leg->SetTextSize(textSize);
    //leg->Draw("same");
    return leg;
}

TLatex *DrawMyLatex(char *text, Double_t x = 0.65, Double_t y = 0.5, Int_t textFont = 62, Size_t textSize = 0.05, Color_t colorIndex = 2)
{
    TLatex *latex = new TLatex(x, y, text);
    latex->SetNDC();
    latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);
    latex->SetTextColor(colorIndex);
    latex->Draw("same");
    return latex;
}

void DrawMyPad(TVirtualPad *pad, const char *xname, const char *yname, float x1, float x2, float y1, float y2)
{

    TH1F *hpad = pad->DrawFrame(x1, y1, x2, y2);
    hpad->GetXaxis()->SetTitle(xname);
    hpad->GetYaxis()->SetTitle(yname);
    hpad->GetXaxis()->SetAxisColor(1);
    hpad->GetYaxis()->SetAxisColor(1);
    hpad->GetXaxis()->SetLabelColor(1);
    hpad->GetYaxis()->SetLabelColor(1);
    hpad->GetXaxis()->SetLabelFont(42);
    hpad->GetYaxis()->SetLabelFont(42);
    hpad->GetXaxis()->SetLabelSize(0.05);
    hpad->GetYaxis()->SetLabelSize(0.05);
    hpad->GetXaxis()->SetLabelOffset(0.01);
    hpad->GetYaxis()->SetLabelOffset(0.01);
    hpad->GetXaxis()->SetTitleFont(42);
    hpad->GetYaxis()->SetTitleFont(42);
    //hpad->GetXaxis()->SetTitleColor( TitleColor);
    //hpad->GetYaxis()->SetTitleColor( TitleColor );
    hpad->GetXaxis()->SetTitleSize(0.06);
    hpad->GetYaxis()->SetTitleSize(0.06);
    hpad->GetXaxis()->SetTitleOffset(1.0);
    hpad->GetYaxis()->SetTitleOffset(1.0);
    pad->Modified();
    pad->Update();
}

void SetMyPad(TVirtualPad *pad, float left, float right, float top, float bottom)
{
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

void DrawMyGraph(TGraph *datagraph, char *xtitle, char *ytitle, Size_t MSize = 1, Style_t MStyle = 28, Color_t MColor = 1, Color_t LColor = 1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t FColor = 16)
{
    datagraph->SetLineColor(LColor);
    datagraph->SetLineWidth(LWidth);
    datagraph->SetLineStyle(LStyle);
    datagraph->SetMarkerSize(MSize);
    datagraph->SetMarkerStyle(MStyle);
    datagraph->SetMarkerColor(MColor);
    datagraph->SetFillColor(FColor);
    //datagraph->SetFillStyle( FStyle );
    datagraph->GetXaxis()->SetTitle(xtitle);
    datagraph->GetYaxis()->SetTitle(ytitle);
    datagraph->GetXaxis()->SetAxisColor(1);
    datagraph->GetYaxis()->SetAxisColor(1);
    datagraph->GetXaxis()->SetLabelColor(1);
    datagraph->GetYaxis()->SetLabelColor(1);
    datagraph->GetXaxis()->SetLabelFont(42);
    datagraph->GetYaxis()->SetLabelFont(42);
    datagraph->GetXaxis()->SetLabelSize(0.05);
    datagraph->GetYaxis()->SetLabelSize(0.05);
    datagraph->GetXaxis()->SetLabelOffset(0.01);
    datagraph->GetYaxis()->SetLabelOffset(0.01);
    datagraph->GetXaxis()->SetTitleFont(42);
    datagraph->GetYaxis()->SetTitleFont(42);
    //datagraph->GetXaxis()->SetTitleColor( TitleColor);
    //datagraph->GetYaxis()->SetTitleColor( TitleColor );
    datagraph->GetXaxis()->SetTitleSize(0.06);
    datagraph->GetYaxis()->SetTitleSize(0.06);
    datagraph->GetXaxis()->SetTitleOffset(1.0);
    datagraph->GetYaxis()->SetTitleOffset(1.0);
}

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
    
