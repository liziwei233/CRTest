#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>

TRandom3 r;

using namespace std;

Double_t outputfunc(Double_t x, vector<double> par, vector<double> tts, int *npe);
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

Double_t outputfunc(Double_t x, vector<double> par, vector<double> tts, int *npe)
{

    //************SiPM S13360-60-25***************
    //***Rise Time 3.6ns
    //***Fall Time 103ns
    //***Gain 2.4e6
    //***SPE Umax =0.35mV
    //***SPE Qmax =0.35mV
    //***TTS =280ps
    //*************************************************
    /*
       if(par.empty()) {
       cout<<"par matrix is empty!"<<endl;
       return 0;
       }
       */
    Double_t val = 0;
    double SPEpar[7];

    double Tmark = 0;
    bool flag;
    //tts = 0;

    //
    // *-----------------SiPM+BJT650 output----------//
    double Trecept = 500e-9; //waiting the photons hit
    double Treject = 500e-9; //recover time,during this time any
    SPEpar[0] = 2.4e6;       //Gain
    SPEpar[1] = 1.6e-19;     //e
    //SPEpar[2]=150e-12;  //Ca  ??
    //SPEpar[3]=3e-12; //C
    SPEpar[2] = 1e-9;   //Ca  ??
    SPEpar[3] = 1.8e-9; //C
    SPEpar[4] = 10e3;   //R   ??
    SPEpar[5] = 50;     //Z
    SPEpar[6] = 1.5e-9; //rise time 1.5ns
    // **
    //

    /*
    //
    // 
    // *-----------------MCP-PMT R3809 output----------//
    double Trecept = 500e-9;	//waiting the photons hit
    double Treject = 500e-9; //recover time,during this time any
    SPEpar[0] = 1.1e6;	 //Gain
    SPEpar[1] = 1.6e-19; //e
    //SPEpar[2]=150e-12;  //Ca  ??
    //SPEpar[3]=3e-12; //C
    SPEpar[2] = 2.65e-12; //Ca  ??
    SPEpar[3] = 12e-12; //C
    SPEpar[4] = 10e3;		//R   ??
    SPEpar[5] = 50;			//Z
    SPEpar[6] = 90e-12; //rise time 1.5ns
    // **
    //
     */

    //int N;
    //N=sizeof(par)/sizeof(par[0]);
    if (par.empty())
        return val;
    sort(par.begin(), par.end());
    int counter = 0;
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
                //r.SetSeed(par.at(n));
                //tts = r.Gaus(0,ttssigma); //TTS of MCP-R3805U
                //cout<<"tts= "<<tts<<endl;
                val += response(x - tts.at(n) - par.at(n), SPEpar);
                counter++;
            }
            else if (par.at(n) - Tmark < (Trecept + Treject))
            {
                val += 0;
            }
            else
            {
                Tmark = par.at(n);
                //r.SetSeed(par.at(n));
                //tts = r.Gaus(0,ttssigma); //TTS of MCP-R3805U
                //cout<<"tts= "<<tts<<endl;
                val += response(x - tts.at(n) - par.at(n), SPEpar);
                counter++;
            }
        }
    }
    //cout<<"n = "<<n<<endl;
    *npe = counter;
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

void MultiTiersOutputfun_SiPM(const char *rootname = "", double fac = 0.2, const char *ParType = "CFD")
{

    TLegend *DrawMyLeg(Double_t xlow = 0.2, Double_t ylow = 0.2, Double_t xup = 0.5, Double_t yup = 0.5, Int_t textFont = 62, Size_t textSize = 0.05);

    TLatex *DrawMyLatex(char *text, Double_t x = 0.65, Double_t y = 0.5, Int_t textFont = 62, Size_t textSize = 0.05, Color_t colorIndex = 2);

    gStyle->SetOptFit(1111);

    //TRandom3 r;
    //r.SetSeed(0);
    char name[1024];
    char buff[1024];

    const int T = 1; //amounts of Tiers

    //Double_t parR[500]={};
    //Double_t parL[500]={};
    vector<vector<double>> parR(T);
    vector<vector<double>> ttsR(T);

    //vector<double> parR2;
    vector<vector<double>> parL(T);
    vector<vector<double>> ttsL(T);
    //vector<double> parL2;

    Double_t RL = -5e-9;
    Double_t RR = 30e-9;
    int binNum = 0;
    binNum = (RR - RL) / 50e-12;

    double ttssigma = 0e-12;
    const int range = (RR - RL) / 25e-12; //25ps/Sample
    Double_t Umax = -1.7;                 //Umax = -1.71mV
    double thrd = 0;
    int TH[6] = {1, 3, 5, 10, 20, 30};
    double Rate = 0;

    bool flagR = 0, flagL = 0;
    int indexL = 0, indexR = 0;
    double keypointL = 0, keypointR = 0;

    double xT0_L[T] = {0}, xT0_R[T] = {0}, xT0[T] = {0};
    double UL[T] = {0}, UR[T] = {0};
    int npeL[T] = {0}, npeR[T] = {0};
    const int certain = 5;

    vector<double> *TR;
    vector<double> *TL;
    vector<int> *IDL;
    vector<int> *IDR;
    TL = new vector<double>;
    TR = new vector<double>;
    IDL = new vector<int>;
    IDR = new vector<int>;
    //count = new vector<int>;
    int N = 0, temp = 0;
    double Temp = 0.;

    double xL[T][range];
    double xR[T][range];
    double yR[T][range];
    double yL[T][range];

    sprintf(name, "%s", rootname);
    sprintf(buff, "%s.root", name);
    cout << "Open the rootfile:" << buff << endl;
    TFile *f1 = new TFile(buff, "READ");
    TTree *t1 = (TTree *)f1->Get("Run");

    t1->SetBranchAddress("PmtR.t", &TR);
    t1->SetBranchAddress("PmtL.t", &TL);
    t1->SetBranchAddress("PmtL.id", &IDL);
    t1->SetBranchAddress("PmtR.id", &IDR);
    sprintf(buff, "%sdata.root", name);
    cout << "output data file name: " << buff << endl;
    TFile *f2 = new TFile(buff, "RECREATE");
    TTree *t2 = new TTree("data", "restore analysed data  from G4");
    sprintf(buff, "UL[%d]/D", T);
    t2->Branch("UL", UL, buff);

    sprintf(buff, "UR[%d]/D", T);
    t2->Branch("UR", UR, buff);

    sprintf(buff, "T0L[%d]/D", T);
    t2->Branch("T0L", xT0_L, buff);

    sprintf(buff, "T0R[%d]/D", T);
    t2->Branch("T0R", xT0_R, buff);

    sprintf(buff, "T0[%d]/D", T);
    t2->Branch("T0", xT0, buff);

    sprintf(buff, "npeL[%d]/I", T);
    t2->Branch("npeL", npeL, buff);

    sprintf(buff, "npeR[%d]/I", T);
    t2->Branch("npeR", npeR, buff);
    //for(int s = 0; s<4;s++){

    double t_L = RL;
    double t_R = RR;

    //f1->cd();

    //TF1 *myFun;
    TH1D *h[2];
    h[0] = new TH1D("hR", "", binNum, RL, RR);
    h[1] = new TH1D("hL", "", binNum, RL, RR);

    TH1D *hSIG = new TH1D("hSIG", "PMT's total signal", binNum, RL, RR);
    TH1D *hSig[2];
    hSig[0] = new TH1D("hRSig", "", binNum, RL, RR);
    hSig[1] = new TH1D("hLSig", "", binNum, RL, RR);

    TF1 *fSIG = new TF1("fSIG", "gaus", RL, RR);
    TF1 *fSig[2];
    fSig[0] = new TF1("fRSig", "gaus", RL, RR);
    fSig[1] = new TF1("fLSig", "gaus", RL, RR);

    N = t1->GetEntries();
    cout << "Entries = " << N << endl;

    //count->
    //for(int i = certain; i < certain+1; i++){
    for (int i = 0; i < N; i++)
    {

        //-----------initial----------------------//
        TL->clear();
        TR->clear();
        IDL->clear();
        IDR->clear();

        h[0]->Reset();
        h[1]->Reset();


        vector<vector<double>>(T).swap(parR);
        vector<vector<double>>(T).swap(parL);
        vector<vector<double>>(T).swap(ttsR);
        vector<vector<double>>(T).swap(ttsL);

        t1->GetEntry(i);
        temp = TR->size();



        for (int k = 0; k < temp; k++)
        {
            //cout<< T[][k] <<endl;
            Temp = (*TR)[k] * 1e-9;

            for (int iT = 0; iT < T; iT++)
            {
             
                if ((*IDR)[k] == iT)
                {

                    parR[iT].push_back(Temp);
                    ttsR[iT].push_back(r.Gaus(0, ttssigma));
               
                }
            }

            if ((*IDR)[k] == 0)
                h[0]->Fill(Temp);

        }

        temp = TL->size();
        for (int k = 0; k < temp; k++)
        {
            Temp = (*TL)[k] * 1e-9;
            for (int iT = 0; iT < T; iT++)
            {
                if ((*IDL)[k] == iT)
                {
                    parL[iT].push_back(Temp);
                    ttsL[iT].push_back(r.Gaus(0, ttssigma));
                }
                
            }
            if ((*IDL)[k] == 0)
                h[1]->Fill(Temp);
        }

        memset(xL, 0, sizeof(xL));
        memset(xR, 0, sizeof(xR));
        memset(yL, 0, sizeof(yL));
        memset(yR, 0, sizeof(yR));
        for (int q = 0; q < T; q++)
        {
            flagR = 1;
            flagL = 1;

            xT0_R[q] = 0;
            xT0_L[q] = 0;
            xT0[q] = 0;

            
            for (int j = 0; j < range; j++)
            {
                xL[q][j] = (RR - RL) / range * j + RL;
                xR[q][j] = xL[q][j];
                yL[q][j] = outputfunc(xL[q][j], parL[q], ttsL[q], npeL + q);
                yR[q][j] = outputfunc(xR[q][j], parR[q], ttsR[q], npeR + q);

                //change the unit from "s" to "ns"
                xL[q][j] = xL[q][j] * 1e9;
                xR[q][j] = xL[q][j];
            }

            UR[q] = TMath::MinElement(range, yR[q]);
            UL[q] = TMath::MinElement(range, yL[q]);

            cout << "Discriminate the signal=========" << endl;
            //* ==============================================
            //* =======Discriminate the signal=========
            //* =================================================
            xT0_R[q] = 0;
            xT0_L[q] = 0;
            xT0[q] = 0;

            flagR = 1;
            flagL = 1;
            indexL = 0;
            indexR = 0;
            keypointL = 0;
            keypointR = 0;
            if (strcmp(ParType, "FIX") == 0) //if the discriminate way is fix threshold discrim
            {
                keypointL = fac;
                keypointR = fac;
            }
            else
            {
                keypointR = fac * UR[q];

                keypointL = fac * UL[q];
            }

            for (int j = 0; j < range; j++)
            {
                if (yR[q][j] < keypointR && flagR)
                {
                    indexR = j;
                    flagR = 0;
                }
                if (yL[q][j] < keypointL && flagL)
                {

                    indexL = j;
                    flagL = 0;
                }
            }

            //*=====================================================
            //* =======Fit the signal and find the timestamp========
            //* ====================================================
            TGraph *gR = new TGraph(range, xR[q], yR[q]);
            TF1 *fitR = pol3fit(gR, xR[q][indexR] - 60e-3, xR[q][indexR] + 80e-3);
            xT0_R[q] = fitR->GetX(keypointR) * 1e-9;

            //return;
            TGraph *gL = new TGraph(range, xL[q], yL[q]);
            TF1 *fitL = pol3fit(gL, xL[q][indexL] - 60e-3, xL[q][indexL] + 80e-3);
            xT0_L[q] = fitL->GetX(keypointL) * 1e-9; 

        }
        hSig[1]->Fill(xT0_L[0]);
        hSig[0]->Fill(xT0_R[0]);
        hSIG->Fill(xT0[0]);

        t2->Fill();

    }

    f1->Close();
    TGraph *gR = new TGraph(range, xR[0], yR[0]);
    TGraph *gL = new TGraph(range, xL[0], yL[0]);

    TCanvas *c = new TCanvas("c", "", 1600, 600);

    //c->cd();
    gPad->Clear();
    c->Divide(2, 1);
    c->cd(1);
    gR->Draw();
    TAxis *xA = gR->GetXaxis();
    TAxis *yA = gR->GetYaxis();
    gR->SetLineWidth(2);
    gR->SetLineColor(kRed);
    xA->SetTitle("Time (s)");
    xA->SetRangeUser(RL, RR);
    xA->SetRangeUser(5, 10);
    yA->SetTitle("Amplitude (mV)");
    //myFun->SetNpx(5e3);
    //myFun->Draw();
    //gR->Draw("AC");

    gL->SetLineWidth(2);
    gL->SetLineColor(kBlue);
    gL->Draw("same");
    //gL->GetXaxis()->SetTitle("Time (s)");
    //gL->Draw("same");

    TLegend *leg;
    leg = DrawMyLeg(0.6, 0.2, 0.75, 0.32);
    leg->AddEntry(gR, "PMT_Right", "lp");
    leg->AddEntry(gL, "PMT_Left", "lp");
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

    TPaveStats *stL = (TPaveStats *)h[1]->FindObject("stats");
    stL->SetY1NDC(0.6);
    stL->SetY2NDC(0.78);
    gPad->Modified();
    gPad->Update();

    

    sprintf(buff, "%s_Signal.png", name);
    c->SaveAs(buff);
    

    TCanvas *c1 = new TCanvas("c1", "", 1600, 600);

    c1->Divide(2, 1);
    c1->cd(1);
    hSig[0]->Draw();
    hSig[0]->Rebin(6);
    hSig[0]->GetXaxis()->SetRangeUser(t_L, t_R);

    hSig[0]->Fit(fSig[0], "R");
    c1->cd(2);
    hSig[1]->Draw();
    hSig[1]->Rebin(6);
    hSig[1]->GetXaxis()->SetRangeUser(t_L, t_R);

    hSig[1]->Fit(fSig[1], "R");
    sprintf(buff, "%s_Twosides_timeresolution.png", name);
    c1->SaveAs(buff);

    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    c2->cd();
    hSIG->Draw();

    hSIG->Rebin(6);
    hSIG->GetXaxis()->SetRangeUser(t_L, t_R);
    hSIG->Fit(fSIG, "R");
    sprintf(buff, "%s_timeresolution.png", name);
    c2->SaveAs(buff);

    f2->cd();
    t2->Write();

    //f2->Close();
    //sprintf(buff,"%s_TimeRes.dat",name);
    cout << "The process is over,THANK YOU!" << endl;

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


