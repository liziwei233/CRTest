#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>

#include "Include/DrawMyfunc.h"
#include "../Include/LZWfunc.h"




void TACor(const char *rootname = "")
{
    TF1 *gausfit(TH1 *h, int rbU, double fac, RANGE *U);
    TF1 *profilefit(TH2 *Rt, double rbU, double rbt, RANGE t, RANGE U, char *name);
    TF1 *CH1Correction(TTree *t1, EVENT *A, RANGE t, RANGE U, int rbU, int rbt, double fac, int iter);
    void CH2Correction(TTree *t1, vector<g4EVENT *> ch, int rbU, int rbt, double fac, int iter);
    TF1 *SumCorrection(TTree *t1, vector<g4EVENT *> ch, const int chN, int rbU, int rbt, double fac, int iter)
    
    //void MygStyle();
    //MygStyle();
    gStyle->SetOptFit(111);

    char name[100];
    char buff[1024];

    //***************************************************//
    //--------------Configuration-----------------------//
    //***************************************************//
    RANGE t = {-1e-9, 2e-9};
    RANGE u = {-15e3, 0};
    RANGE posx = {-20, 20};
    RANGE posy = posx;
    int rbt = 2;
    int rbu = 16;

    // the range set After Correct
    //double L2 = -100e-12;
    //double R2 = 100e-12;
    int bint = (t.R - t.L) / 1e-12;
    int binu = (u.R - u.L) / 1;
    int binpx = (posx.R - posx.L) / 1;
    int binpy = (posy.R - posy.L) / 1;
    const int T = 1;    // the number of tiers
    const int iter = 5; //the number of iteration
    double fac = 2;
    //---------------------------------------------------//
    //***************************************************//

    cout << "Start excuate TA correction procedure =====>>>>" << endl;
    //return;

    //cout<<"<<---- Succeed excuating ---->>"<<endl;
    //return;
    //thrd = (s+1)*0.2;
    sprintf(name, "%s", rootname);
    sprintf(buff, "%s.root", name);
    TFile *f = new TFile(buff, "READ");
    TTree *tr = (TTree *)f->Get("data");
    //double UL[T],UR[T],T0L[T],T0R[T],T0[T];
    double simutime, PEtime, CFD, LED, Amp;
    int npe;
    g4EVENT event[2];
    POSITION pos;
    double Utemp = 0, Ttemp = 0;
    double T0_cor = 0;
    float TT[1000000] = {0.};

    tr->SetBranchAddress("UL", &event[0].Amp);
    tr->SetBranchAddress("UR", &event[1].Amp);
    tr->SetBranchAddress("T0L", &event[0].time);
    tr->SetBranchAddress("T0R", &event[1].time);
    tr->SetBranchAddress("T01stpeL", &event[0].PEtime);
    tr->SetBranchAddress("T01stpeR", &event[1].PEtime);
    tr->SetBranchAddress("npeL", &event[0].npe);
    tr->SetBranchAddress("npeR", &event[1].npe);
    tr->SetBranchAddress("InX", &pos.x);
    tr->SetBranchAddress("InY", &pos.y);
    
    /*

*/
    //
    //*  simu time resolution of left PMT
    double TRL = 0;
    TF1 *fsL = CH1Correction(tr, event[0], t, u, rbu, rbt, fac, iter);
    if (fsL)
        TRL = fsL->GetParameter(2);

    //
    //*  simu time resolution of right PMT
    double TRR = 0;
    TF1 *fsR = CH1Correction(tr, event[1], t, u, rbu, rbt, fac, iter);
    if (fsR)
        TRR = fsR->GetParameter(2);

    //
    //*  simu time resolution of  left - right
    double TRLMR = 0;
    TF1 *fsLMR = CH2Correction(tr, event, rbu, rbt, fac, iter);
    if (fsLMR)
        TRLMR = fsLMR->GetParameter(2);

    //
    //* simu time resolution of (left+right)/2
    double TRLPR = 0;
    TF1 *fsLPR = SumCorrection(tr, event, 2, rbu, rbt, fac, int iter);
    if (fsLPR)
        TRLPR = fsLPR->GetParameter(2);

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    TCanvas *c3 = new TCanvas("c3", "", 800, 600);

    TCanvas *cU = new TCanvas("cU", "", 800, 600);

    TH1D *hL = new TH1D("hL", ";Time (s);Counts", bint, t.L, t.R);
    TH1D *hR = new TH1D("hR", ";Time (s);Counts", bint, t.L, t.R);
    TH1D *hLPR = new TH1D("hLPR", ";Time (s);Counts", bint, t.L, t.R);
    TH1D *hLMR = new TH1D("hLMR", ";Time (s);Counts", bint, t.L, t.R);
    TH2D *hpos = new TH2D("hpos", ";x (mm);y (mm)", binpx, posx.L, posx.R, binpy, posy.L, posy.R);

    //* intrinsic timeresolution and hit position

    int N = 0;
    N = tr->GetEntries();
    cout << "the Entries is :" << N << endl;
    for (int i = 0; i < N; i++)
    {
        tr->GetEntry(i);
        hL->Fill(event[0].PEtime);
        hR->Fill(event[1].PEtime);
        hLPR->Fill((event[0].PEtime + event[1].PEtime) / 2);
        hLMR->Fill(event[0].PEtime - event[1].PEtime);
        hpos->Fill(pos.x, pos.y);
    }
    c1->cd();
    TF1 *fL = gausfit(hL, rbt, fac, &t);

    double inTR_L = 0;
    if (fL)
        inTR_L = fL->GetParameter(2);
    else
    {
        fL = new TF1("fL", "landau", t.L, t.R);
        TFitResultPtr failed = hL->Fit("R");
        if (!failed)
            inTR_L = fL->GetParameter(2);
        else
            cout << "fit left: gaus and landau fit both failed !" << endl;
    }

    sprintf(buff, "%s_intrinsic_TR_left.png", name);
    c1->SaveAs(buff);

    c1->Clear();
    TF1 *fR = gausfit(hR, rbt, fac, &t);
    double inTR_R = 0;
    if (fR)
        inTR_R = fR->GetParameter(2);
    else
    {
        fR = new TF1("fR", "landau", t.L, t.R);
        TFitResultPtr failed = hR->Fit("R");
        if (!failed)
            inTR_R = fR->GetParameter(2);
        else
            cout << "fit right: gaus and landau fit both failed !" << endl;
    }

    sprintf(buff, "%s_intrinsic_TR_right.png", name);
    c1->SaveAs(buff);

    c1->Clear();
    TF1 *fLPR = gausfit(hLPR, rbt, fac, &t);
    double inTR_LPR = 0;
    if (fLPR)
        inTR_LPR = fLPR->GetParameter(2);
    else
    {
        fLPR = new TF1("fLPR", "landau", t.L, t.R);
        TFitResultPtr failed = hLPR->Fit("R");
        if (!failed)
            inTR_LPR = fLPR->GetParameter(2);
        else
            cout << "fit (left+right)/2: gaus and landau fit both failed !" << endl;
    }

    sprintf(buff, "%s_intrinsic_TR_LPR.png", name);
    c1->SaveAs(buff);

    c1->Clear();
    TF1 *fLMR = gausfit(hLMR, rbt, fac, &t);
    double inTR_LMR = 0;
    if (fLMR)
        inTR_LMR = fLMR->GetParameter(2);
    else
    {
        fLMR = new TF1("fLMR", "landau", t.L, t.R);
        TFitResultPtr failed = hLMR->Fit("R");
        if (!failed)
            inTR_LMR = fLMR->GetParameter(2);
        else
            cout << "fit left-right: gaus and landau fit both failed !" << endl;
    }

    sprintf(buff, "%s_intrinsic_TR_LMR.png", name);
    c1->SaveAs(buff);

    c2->cd();
    hpos->Draw("colz");
    sprintf(buff, "%s_hitpostion.png", name);
    c2->SaveAs(buff);
}

TF1 *gausfit(TH1 *h, int rbU, double fac, RANGE *U)
{

    double mean = 0;
    double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser((*U).L, (*U).R);
    TF1 *fitU = new TF1("fitU", "gaus", (*U).L, (*U).R);
    mean = hU->GetBinCenter(hU->GetMaximumBin());
    fitU->SetParameter(1, mean);
    cout << mean << "\t" << sigma << endl;
    TFitResultPtr failed = hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;

    if (failed)
        return fitU = 0;

    else
    {

        if ((*U).L < mean - 8 * sigma)
            (*U).L = mean - 8 * sigma;
        if ((*U).R > mean + 8 * sigma)
            (*U).R = mean + 8 * sigma;

        hU->GetXaxis()->SetRangeUser((*U).L, (*U).R);

        return fitU;
    }
}

TF1 *profilefit(TH2 *Rt, double rbU, double rbt, RANGE t, RANGE U, char *name)
{

    TCanvas *cpfx = new TCanvas("cpfx", "cpfx", 1600, 600);
    cpfx->Clear();
    cpfx->Divide(2, 1);
    cpfx->cd(1);
    TH2 *Qt = (TH2 *)Rt->Clone();
    //TH2* Qt = (TH2*) Rt->Clone("tmp");

    Qt->Draw("colz");

    Qt->RebinX(rbU);
    Qt->RebinY(rbt);
    Qt->GetYaxis()->SetRangeUser(t.L, t.R);
    Qt->GetXaxis()->SetRangeUser(U.L, U.R);
    Qt->ProfileX();

    cpfx->cd(2);
    TH1 *Qpfx = Qt->ProfileX();
    //Qpfx->Reset();

    //Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
    Qpfx->Draw();
    Qpfx->GetYaxis()->SetRangeUser(t.L, t.R);
    Qpfx->GetXaxis()->SetRangeUser(U.L, U.R);

    TF1 *fitQt = new TF1("fitQt", "[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)", U.L, U.R);

    Qpfx->Fit(fitQt, "R");

    sprintf(buff, "%s.png", name);
    cpfx->SaveAs(buff);
    Qpfx->Reset();
    //delete Qt;
    //delete Qpfx;
    //delete c5;
    return fitQt;
}

TF1 *CH1Correction(TTree *t1, EVENT *A, RANGE t, RANGE U, int rbU, int rbt, double fac, int iter)
{
    //A channel with trigger
    //write for SIMU

    //Correction QT
    double Q;
    RANGE initial_t;
    initial_t.L = -5;
    initial_t.R = 5;
    sprintf(buff, "%s_ch1Correction.dat", name.c_str());
    output.open(buff, ios::trunc);
    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);

    TH1D *ht;
    TH2D *hAT;

    TH1D *htraw = new TH1D("htraw", "time Resolution;T (ps);Counts", 2e4, t.L, t.R);

    TH2D *hATraw = new TH2D("hATraw", "", 600, U.L, U.R, 2e4, t.L, t.R);

    TH1D *htcor = new TH1D("htcor", "time Resolution;T (ps);Counts", 2e4, initial_t.L, initial_t.R);

    TH2D *hATcor = new TH2D("hATcor", "", 600, U.L, U.R, 2e4, initial_t.L, initial_t.R);

    TF1 *fitAT = new TF1("fitAT", "0", U.L, U.R);
    TF1 *fitT;

    double TAcor = 0;
    double T[1000000] = {0};
    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;
        for (int i = 0; i < N; i++)
        {
            t1->GetEntry(i);

            Q = (*A).Amp;
            TAcor = fitAT->Eval(Q);
            //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3){
            //if((*A).Q>50&&(*A).A<1.05&&(*B).Q>50&&(*B).A<1.05&&(*MCP).Q>0.3&&(*A).rise<4&&(*B).rise<4){
            // else if((*A).A<0&&(*A).time>0){
            if ((*A).x > charcut.x.L &&
                (*A).x < charcut.x.R &&
                (*A).Q > charcut.q.L &&
                (*A).Q < charcut.q.R &&
                (*A).rise > charcut.r.L &&
                (*A).rise < charcut.r.R &&
                (*A).bl > charcut.bl.L &&
                (*A).bl < charcut.bl.R &&
                (*A).blrms > charcut.blrms.L &&
                (*A).blrms < charcut.blrms.R)
            {

                if ((*A).Amp > charcut.y.L &&
                    (*A).Amp < charcut.y.R)
                {
                    if (s == 0)
                        //T[i]=(*A).time-(*MCP).time-TAcor;
                        htraw->Fill(T[i]);
                    hATraw->Fill(Q, T[i]);
                }
                else
                {
                    T[i] = T[i] - TAcor;
                    htcor->Fill(T[i]);
                    hATcor->Fill(Q, T[i]);
                }
                //				cout<<"Q="<<Q<<" T="<<T[i]<<endl;
            }
        }
    }

    //cout<<"<====process check====>"<<endl;

    cCor->cd();
    cCor->Clear();
    if (s == 0)
    {
        ht = (TH1D *)htraw->Clone();
        hAT = (TH2D *)hATraw->Clone();
    }
    else
    {
        ht = (TH1D *)htcor->Clone();
        hAT = (TH2D *)hATcor->Clone();
        cout << "progress check" << endl;
    }
    //return;
    //t.L=-5;
    //t.R=5;
    t.L = ht->GetMean() - 5 * ht->GetRMS();
    t.R = ht->GetMean() + 5 * ht->GetRMS();
    //cout<<"iter="<<s<<", MEAN="<<ht->GetMean()<<"\t RMS="<<ht->GetRMS()<<endl;
    //cout<<"iter="<<s<<", tL="<<tL<<"\t tR="<<tR<<endl;
    //if (s==1) return;
    //twoguasfit(ht,&tL,&tR,fac,rbt);
    fitT = gausfit(ht, rbt, fac, &t);
    if (!fitT)
    {
        fitT = new TF1("fitT", "landau", t.L, t.R);
        TFitResultPtr failed = ht->Fit("R");
        if (failed)
        {
            fitT = 0;
            cout << "gaus and landau fit both failed !" << endl;
        }
    }
    sprintf(buff, "%s_TR_cor%d.png", name.c_str(), s);
    cCor->SaveAs(buff);

    if (!fitT)
        output << name.c_str() << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
    else
        output << name.c_str() << "\t" << s << "\t"
               << "0000"
               << "\t"
               << "0000" << endl;

    fitAT = profilefit(hAT, rbU, rbt * 2, t, U, buff);

    ht->Reset();
    hAT->Reset();
    htraw->Reset();
    hATraw->Reset();
    htcor->Reset();
    hATcor->Reset();
}

return fitT;
}

TF1 *CH2Correction(TTree *t1, vector<EVENT *> ch, vector<charRANGE> range, vector<CUT> cut, OPTION opt, string name, ofstream& op)
{
//* 
//** The number of channels : 2
//** Discrimination: single
//* 

    //pass your parameter
    RANGE t=range.at(0).t;
    vector<RANGE> U;
    U.pushback(range.at(0).y);
    U.pushback(range.at(1).y);

    RANGE initial_t;
    initial_t.L = t.L;
    initial_t.R = t.R;

    int iter = opt.iter;
    double fac = opt.fac;
    int rbt = opt.rbt;
    int rbU = opt.rbu;

    ofstream output;
    sprintf(buff, "%s.dat", name.c_str());
    output.open(buff, ios::trunc);
    //cout<<buff<<" has been built"<<endl;

    TCanvas *cCor = new TCanvas("cCor", "cCor", 800, 600);
    TH1D *ht = new TH1D("ht", "Time Resolution;T (ps);Counts", 1e3, t.L, t.R);

    TF1 *fitT;

    double TAcor = 0;

    double T[1000000] = {0};

    //EVENT** ch[chN];
    const int chN=2;

    TH2D *hAT[chN];
    TF1 *fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);
    for (int h = 0; h < chN; h++)
    {
        sprintf(buff, "hAT%d", h);
        hAT[h] = new TH2D(buff, buff, 200, U.at(h).L, U.at(h).R, 4e3, t.L, t.R);
    }

    //ch[0]=&A;
    //ch[1]=&B;
    //ch[2]=&MCP;
    //double L=20;
    //cout<<"the number of values of fraction of CFD: "<<L<<endl;
    //for(int j=0;j<L;j++){
    //multi CFD fraction dicrimination

    for (int s = 0; s < iter; s++)
    {

        int N = t1->GetEntries();
        cout << "Total entries is :" << N << endl;

        //correction of A1
        for (int h = 0; h < chN; h++)
        {
            t.L = initial_t.L;
            t.R = initial_t.R;
            for (int i = 0; i < N; i++)
            {
                t1->GetEntry(i);

                //Q1=(*A).Q;
                //Q2=(*B).Q;
                if (h - 1 < 0)
                    TAcor = fitAT->Eval((*ch.at(chN - 1)).Q);
                else
                    TAcor = fitAT->Eval((*ch.at(h - 1)).Q);
                int il = 0; // if all boolean conditions are true , il==chN.
                for (il = 0; il < chN; il++)
                {
#ifndef G4_FLAG
                    if (!((*ch.at(il)).Amp > cut.at(il).y.L &&
                          (*ch.at(il)).Amp < cut.at(il).y.R &&
                          (*ch.at(il)).Q > cut.at(il).q.L &&
                          (*ch.at(il)).Q < cut.at(il).q.R &&
                          (*ch.at(il)).rise > cut.at(il).r.L &&
                          (*ch.at(il)).rise < cut.at(il).r.R && 
                          (*ch.at(il)).bl > cut.at(il).bl.L &&
                          (*ch.at(il)).bl < cut.at(il).bl.R &&
                          (*ch.at(il)).blrms > cut.at(il).blrms.L &&
                          (*ch.at(il)).blrms < cut.at(il).blrms.R &&
                          (*ch.at(il)).x > cut.at(il).x.L &&
                          (*ch.at(il)).x < cut.at(il).x.R ))
                        break;
#else
                    if (!((*ch.at(il)).Amp > cut.at(il).u.L &&
                          (*ch.at(il)).Amp < cut.at(il).u.R &&
                          (*ch.at(il)).time > cut.at(il).t.L &&
                          (*ch.at(il)).time < cut.at(il).t.R &&
                          (*ch.at(il)).npe > cut.at(il).npe.L &&
                          (*ch.at(il)).npe < cut.at(il).npe.R ))
                        break;
#endif
                }
                if (il == chN)
                {

                    if (s == 0 && h == 0)
                    {
#ifdef ADD
                        // average the time of all channels ;
                        for (int iT = 0; iT < chN; iT++)
                        {
                            T[i] += (*ch.at(iT)).time;
                        }
                        T[i] = T[i] / chN;
#else 
                        T[i] = (*ch.at(0)).time-(*ch.at(1)).time;
#endif
                    }
                    else
                        T[i] = T[i] - TAcor;

                    ht->Fill(T[i]);
                    hAT[h]->Fill((*ch.at(h)).Q, T[i]);
                }
            }

            cCor->cd();
            cCor->Clear();
            ht->Draw();
            //return;
            //tRL=-5;
            //tRR=5;

            //twoguasfit(ht,&tL,&tR,fac,rbt);
            fitT = gausfit(ht, rbt, fac, &t);
            //	sprintf(buff,"%s_TR_ch%d_cor%d.png",name,h,s);
            sprintf(buff, "%s_ch%d_TR_cor%d.png", name.c_str(), h, s);
            cCor->SaveAs(buff);

            sprintf(buff, "%s_ch%d_At_pfx_cor%d", name.c_str(), h, s);
            //sprintf(buff,"%s_At_ch%d_pfx_cor%d",name,h,s);
            
            //*
            //* fit your histgram by gaus-function or landau function.
            fitT = gausfit(ht, rbt, fac, &t);
            if (!fitT)
            {
                fitT = new TF1("fitT", "landau", t.L, t.R);
                TFitResultPtr failed = ht->Fit("R");
                if (failed)
                {
                    fitT = 0;
                    cout << "gaus and landau fit both failed !" << endl;
                }
            }
            sprintf(buff, "%s_ch%d_TR_cor%d.png", name.c_str(), h, s);
            cCor->SaveAs(buff);

            //*
            //* output your fit parameter if fit successed.
            if (!fitT)
                output << name.c_str() << "\t" << h << "\t" << s << "\t" << fitT->GetParameter(2) << "\t" << fitT->GetParError(2) << endl;
            else
                output << name.c_str() << "\t" << h << "\t" << s << "\t"
                       << "0000"
                       << "\t"
                       << "0000" << endl;
            fitAT = profilefit(hAT[h], rbU, rbt * 8, t, U.at(h), buff);

            ht->Reset();
            hAT[h]->Reset();
        }
    }
    fitAT = new TF1("fitAT", "0", U.at(0).L, U.at(0).R);

    //}

    output.close();
    return fitT;
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
