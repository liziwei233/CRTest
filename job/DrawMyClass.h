#ifndef DrawMyClass_h
#define DrawMyClass_h

#include <stdio.h>
#include <vector>
#include <istream>
#include <ostream>
#include <string>
#include <cstdlib>

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TVirtualPad.h"

//using namespace std;

char buff[1024];
char line[180];
int counter = 0;

//TFile *sfile = NULL;
//TFile *prefile = NULL;

typedef struct POSITION
{
    double x;
    double y;
    double z;
} DIRECTION;

struct RANGE
{
    double L;
    double R;
};

void setgStyle()
{

    gStyle->SetFrameLineWidth(3);
    //gStyle->SetFrameBorderSize(2);
    gStyle->SetTickLength(0.04);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1111);
    gStyle->SetErrorX(1);
    gStyle->SetEndErrorSize(4);
}
void Pal2()
{
    Int_t MyPalette[100];
    Double_t Red[] = {0., 0.0, 1.0, 1.0, 1.0};
    Double_t Green[] = {0., 0.0, 0.0, 1.0, 1.0};
    Double_t Blue[] = {0., 1.0, 0.0, 0.0, 1.0};
    Double_t Length[] = {0., .25, .50, .75, 1.0};
    Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 100);
    for (int i = 0; i < 100; i++)
        MyPalette[i] = FI + i;

    gStyle->SetPalette(100, MyPalette, 0.01);
}

void Pal()
{
    TColor *col[5];
    int j = 0;
    col[j++] = gROOT->GetColor(kRed - 4);
    col[j++] = gROOT->GetColor(kOrange + 1);
    //col[j++] = gROOT->GetColor(kYellow-7);
    col[j++] = gROOT->GetColor(kSpring + 8);
    //col[j++] = gROOT->GetColor(kGreen+1);
    col[j++] = gROOT->GetColor(kCyan - 7);
    col[j++] = gROOT->GetColor(kAzure - 1);
    //col[j++] = gROOT->GetColor(kBlue-6);
    double r[5], g[5], b[5];
    float R, G, B;
    for (int i = 0; i < 5; i++)
    {
        col[i]->GetRGB(R, G, B);
        r[i] = R;
        g[i] = G;
        b[i] = B;
    }
    Int_t MyPalette[100];
    Double_t Length[] = {0., .25, .50, .75, 1.0};
    Int_t FI = TColor::CreateGradientColorTable(5, Length, (double *)r, (double *)g, (double *)b, 100);
    for (int i = 0; i < 100; i++)
        MyPalette[i] = FI + i;

    gStyle->SetPalette(100, MyPalette, 0.01);
    TColor::InvertPalette();
}

void DrawMyGraph(TGraph *datagraph, char *xtitle, char *ytitle, Size_t MSize = 1, Style_t MStyle = 28, Color_t MColor = 1, Color_t LColor = 1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t FColor = 0)
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
    datagraph->GetXaxis()->SetLabelSize(0.07);
    datagraph->GetYaxis()->SetLabelSize(0.07);
    datagraph->GetXaxis()->SetLabelOffset(0.01);
    datagraph->GetYaxis()->SetLabelOffset(0.01);
    datagraph->GetXaxis()->SetTitleFont(42);
    datagraph->GetYaxis()->SetTitleFont(42);
    //datagraph->GetXaxis()->SetTitleColor( TitleColor);
    //datagraph->GetYaxis()->SetTitleColor( TitleColor );
    datagraph->GetXaxis()->SetTitleSize(0.07);
    datagraph->GetYaxis()->SetTitleSize(0.07);
    datagraph->GetXaxis()->SetTitleOffset(1.0);
    datagraph->GetYaxis()->SetTitleOffset(1.0);
    datagraph->GetXaxis()->CenterTitle(1);
    datagraph->GetYaxis()->CenterTitle(1);
}
void DrawMyHist(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Style_t LStyle = 1, Color_t TitleColor = 1)
{
    datahist->SetLineColor(LColor);
    datahist->SetLineWidth(LWidth);
    datahist->SetLineStyle(LStyle);

    if (strlen(xtitle))
        datahist->GetXaxis()->SetTitle(xtitle);
    if (strlen(ytitle))
        datahist->GetYaxis()->SetTitle(ytitle);

    datahist->GetYaxis()->SetMaxDigits(3);
    datahist->GetXaxis()->SetAxisColor(1);
    datahist->GetYaxis()->SetAxisColor(1);
    datahist->GetXaxis()->SetLabelColor(1);
    datahist->GetYaxis()->SetLabelColor(1);
    datahist->GetXaxis()->SetLabelFont(42);
    datahist->GetYaxis()->SetLabelFont(42);
    datahist->GetXaxis()->SetLabelSize(0.07);
    datahist->GetYaxis()->SetLabelSize(0.07);
    datahist->GetXaxis()->SetLabelOffset(0.01);
    datahist->GetYaxis()->SetLabelOffset(0.01);
    datahist->GetXaxis()->SetTitleFont(42);
    datahist->GetYaxis()->SetTitleFont(42);
    datahist->GetXaxis()->SetTitleColor(TitleColor);
    datahist->GetYaxis()->SetTitleColor(TitleColor);
    datahist->GetXaxis()->SetTitleSize(0.07);
    datahist->GetYaxis()->SetTitleSize(0.07);
    datahist->GetXaxis()->SetTitleOffset(1.0);
    datahist->GetYaxis()->SetTitleOffset(1.0);
    //datahist->GetXaxis()->SetBorderSize(5);

    datahist->GetXaxis()->SetNdivisions(510);
    datahist->GetYaxis()->SetNdivisions(510);
    datahist->GetXaxis()->CenterTitle(1);
    datahist->GetYaxis()->CenterTitle(1);

     
}

void DrawMy2dHist(TH2 *datahist, char *xtitle, char *ytitle, Style_t MStyle = 8, Color_t MColor = 1, Size_t MSize = 1, Color_t TitleColor = 1)
{
    TPaveStats *stat = (TPaveStats *)datahist->FindObject("stats");
    datahist->SetStats(0);

    datahist->SetMarkerStyle(MStyle);
    datahist->SetMarkerColor(MColor);
    datahist->SetMarkerSize(MSize);

    datahist->SetLineColor(MColor);
    datahist->SetLineWidth(1);
    if (strlen(xtitle))
        datahist->GetXaxis()->SetTitle(xtitle);
    if (strlen(ytitle))
        datahist->GetYaxis()->SetTitle(ytitle);

    datahist->GetXaxis()->SetAxisColor(1);
    datahist->GetYaxis()->SetAxisColor(1);
    datahist->GetXaxis()->SetLabelColor(1);
    datahist->GetYaxis()->SetLabelColor(1);
    datahist->GetXaxis()->SetLabelFont(42);
    datahist->GetYaxis()->SetLabelFont(42);
    datahist->GetXaxis()->SetLabelSize(0.07);
    datahist->GetYaxis()->SetLabelSize(0.07);
    datahist->GetXaxis()->SetLabelOffset(0.01);
    datahist->GetYaxis()->SetLabelOffset(0.01);
    datahist->GetXaxis()->SetTitleFont(42);
    datahist->GetYaxis()->SetTitleFont(42);
    datahist->GetXaxis()->SetTitleColor(TitleColor);
    datahist->GetYaxis()->SetTitleColor(TitleColor);
    datahist->GetXaxis()->SetTitleSize(0.07);
    datahist->GetYaxis()->SetTitleSize(0.07);
    datahist->GetXaxis()->SetTitleOffset(1.0);
    datahist->GetYaxis()->SetTitleOffset(1.0);
    //datahist->GetXaxis()->SetBorderSize(5);
    datahist->GetXaxis()->SetNdivisions(510);
    datahist->GetYaxis()->SetNdivisions(510);
    datahist->GetXaxis()->CenterTitle(1);
    datahist->GetYaxis()->CenterTitle(1);
}
void DrawMyfun(TF1 *datafunc, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Style_t LStyle = 1)
{
    datafunc->SetLineColor(LColor);
    datafunc->SetLineWidth(LWidth);
    datafunc->SetLineStyle(LStyle);
    datafunc->GetXaxis()->SetTitle(xtitle);
    datafunc->GetYaxis()->SetTitle(ytitle);
    datafunc->GetXaxis()->SetAxisColor(1);
    datafunc->GetYaxis()->SetAxisColor(1);
    datafunc->GetXaxis()->SetLabelColor(1);
    datafunc->GetYaxis()->SetLabelColor(1);
    datafunc->GetXaxis()->SetLabelFont(42);
    datafunc->GetYaxis()->SetLabelFont(42);
    datafunc->GetXaxis()->SetLabelSize(0.07);
    datafunc->GetYaxis()->SetLabelSize(0.07);
    datafunc->GetXaxis()->SetLabelOffset(0.01);
    datafunc->GetYaxis()->SetLabelOffset(0.01);
    datafunc->GetXaxis()->SetTitleFont(42);
    datafunc->GetYaxis()->SetTitleFont(42);
    datafunc->GetXaxis()->SetTitleColor(1);
    datafunc->GetYaxis()->SetTitleColor(1);
    datafunc->GetXaxis()->SetTitleSize(0.07);
    datafunc->GetYaxis()->SetTitleSize(0.07);
    datafunc->GetXaxis()->SetTitleOffset(1.0);
    datafunc->GetYaxis()->SetTitleOffset(1.0);
    //datahist->GetXaxis()->SetBorderSize(5);
    datafunc->GetXaxis()->SetNdivisions(510);
    datafunc->GetYaxis()->SetNdivisions(510);
    datafunc->GetXaxis()->CenterTitle(1);
    datafunc->GetYaxis()->CenterTitle(1);
     
}
void Drawxline(float x1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t LColor = 6)
{

    //TLine* line1 = new TLine(x1,gPad->VtoPixel(gPad->GetUymin()),x1,gPad->VtoPixel(gPad->GetUymax()));
    //TLine* line2 = new TLine(x2,gPad->VtoPixel(gPad->GetUymin()),x2,gPad->VtoPixel(gPad->GetUymax()));
    gPad->Update();
    gPad->Modified();
    TLine *line1 = new TLine(x1, gPad->GetUymin(), x1, gPad->GetUymax());

    line1->SetLineWidth(LWidth);
    line1->SetLineStyle(LStyle);
    line1->SetLineColor(LColor);
    line1->Draw("same");
}
void Drawyline(float y1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t LColor = 6)
{

    //TLine* line1 = new TLine(x1,gPad->VtoPixel(gPad->GetUymin()),x1,gPad->VtoPixel(gPad->GetUymax()));
    //TLine* line2 = new TLine(x2,gPad->VtoPixel(gPad->GetUymin()),x2,gPad->VtoPixel(gPad->GetUymax()));
    gPad->Update();
    gPad->Modified();
    TLine *line1 = new TLine(gPad->GetUxmin(), y1, gPad->GetUxmax(), y1);

    line1->SetLineWidth(LWidth);
    line1->SetLineStyle(LStyle);
    line1->SetLineColor(LColor);
    line1->Draw("same");
}

TLegend *DrawMyLeg(Double_t xlow = 0.3, Double_t ylow = 0.6, Double_t xup = 0.6, Double_t yup = 0.9, Int_t textFont = 62, Double_t textSize = 0.03)
{
    TLegend *leg = new TLegend(xlow, ylow, xup, yup);
    //leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(10);
    leg->SetTextFont(textFont);
    leg->SetTextSize(textSize);
    //leg->Draw("same");
    return leg;
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
TLatex *DrawMyLatex(char *text, Double_t x = 0.65, Double_t y = 0.5, Int_t textFont = 42, Size_t textSize = 0.07, Color_t colorIndex = 1)
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
TCanvas *cdC(int n, double left = 0.15, double right = 0.1, double up = 0.1, double down = 0.15)
{
    sprintf(buff, "c%d", n);
    TCanvas *c = new TCanvas(buff, buff, 800, 600);
    c->cd();
    //gPad->SetGrid();
    //SetMyPad(gPad, 0.15, 0.05, 0.1, 0.14);
    SetMyPad(gPad, left, right, up, down);
    return c;
}
TCanvas *cd2C(int n, double left = 0.15, double right = 0.1, double up = 0.1, double down = 0.15)
{
    sprintf(buff, "c%d", n);
    TCanvas *c = new TCanvas(buff, buff, 1600, 600);
    c->Divide(2, 1);
    //c->cd();
    //gPad->SetGrid();
    //SetMyPad(gPad, 0.15, 0.05, 0.1, 0.14);
    SetMyPad(gPad, left, right, up, down);
    return c;
}
//down->up,left->right
void CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
    if (!C)
        return;
    // Setup Pad layout:
    Float_t vSpacing = 0.0;
    Float_t vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny;
    Float_t hSpacing = 0.0;
    Float_t hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx;
    Float_t vposd, vposu, vmard, vmaru, vfactor;
    Float_t hposl, hposr, hmarl, hmarr, hfactor;
    for (Int_t i = 0; i < Nx; i++)
    {
        if (i == 0)
        {
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr - hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
        }
        else if (i == Nx - 1)
        {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr - hposl;
            hmarl = 0.0;
            hmarr = rMargin / (hposr - hposl);
        }
        else
        {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hfactor = hposr - hposl;
            hmarl = 0.0;
            hmarr = 0.0;
        }
        for (Int_t j = 0; j < Ny; j++)
        {
            if (j == 0)
            {
                vposd = 0.0;
                vposu = bMargin + vStep;
                vfactor = vposu - vposd;
                vmard = bMargin / vfactor;
                vmaru = 0.0;
            }
            else if (j == Ny - 1)
            {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep + tMargin;
                vfactor = vposu - vposd;
                vmard = 0.0;
                vmaru = tMargin / (vposu - vposd);
            }
            else
            {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep;
                vfactor = vposu - vposd;
                vmard = 0.0;
                vmaru = 0.0;
            }
            C->cd(0);
            char name[16];
            sprintf(name, "pad_%i_%i", i, j);
            TPad *pad = (TPad *)gROOT->FindObject(name);
            if (pad)
                delete pad;
            pad = new TPad(name, "", hposl, vposd, hposr, vposu);
            pad->SetLeftMargin(hmarl);
            pad->SetRightMargin(hmarr);
            pad->SetBottomMargin(vmard);
            pad->SetTopMargin(vmaru);
            pad->SetFrameBorderMode(0);
            pad->SetBorderMode(0);
            pad->SetBorderSize(0);
            pad->Draw();
        }
    }
}

TCanvas *DrawMyCanvas(TH1 **h, int Nx = 2, int Ny = 2)
{
    //Example of canvas partitioning
    // Sometimes the Divide() method is not appropriate to divide a Canvas.
    // Because of the left and right margins, all the pads do not have the
    // same width and height. CanvasPartition does that properly. This
    // example also ensure that the axis labels and titles have the same
    // sizes and that the tick marks length is uniform.
    //Author:
    gStyle->SetOptStat(0);
    TCanvas *C = (TCanvas *)gROOT->FindObject("C");
    if (C)
        delete C;
    C = new TCanvas("C", "canvas", 1024, 640);
    C->SetFillStyle(4000);
    //TH1F *hFrame = (TH1F*) (*(h+1))->Clone();
    //  hFrame->Draw();
    //  return C;
    // Margins
    Float_t lMargin = 0.12;
    Float_t rMargin = 0.05;
    Float_t bMargin = 0.15;
    Float_t tMargin = 0.05;

   
    // Canvas setup
    CanvasPartition(C, Nx, Ny, lMargin, rMargin, bMargin, tMargin);
    // Dummy histogram.
    TPad *pad[Nx][Ny];
    for (Int_t i = 0; i < Nx; i++)
    {
        for (Int_t j = 0; j < Ny; j++)
        {
            C->cd(0);
            // Get the pads previosly created.
            char pname[16];
            sprintf(pname, "pad_%i_%i",  i,j);
            pad[i][j] = (TPad *)gROOT->FindObject(pname);
            pad[i][j]->Draw();
            pad[i][j]->SetFillStyle(4000);
            pad[i][j]->SetFrameFillStyle(4000);
            pad[i][j]->cd();
            // Size factors
            Float_t xFactor = pad[0][0]->GetAbsWNDC() / pad[i][j]->GetAbsWNDC();
            Float_t yFactor = pad[0][0]->GetAbsHNDC() / pad[i][j]->GetAbsHNDC();
            Float_t BW = (*(h + (Ny-1-j) * 2 + i))->GetBinWidth(0);
            char hname[16];
            sprintf(hname, "h%i", (Ny-1-j) * 2 + i);
            TH1F *hFrame = (TH1F *)(*(h + (Ny-1-j) * 2 + i))->Clone(hname);
            hFrame->Reset();
            hFrame->Draw();
            // y axis range
            hFrame->GetYaxis()->SetRangeUser(1, 1.2 * (*(h + (Ny-1-j) * 2 + i))->GetMaximum());
            // x axis range
            hFrame->GetXaxis()->SetRangeUser( (*(h + (Ny-1-j) * 2 + i))->GetXaxis()->GetXmin()+100 * BW, (*(h + (Ny-1-j) * 2 + i))->GetXaxis()->GetXmax()-100 * BW);
            // Format for y axis
            hFrame->GetYaxis()->SetLabelFont(42);
            hFrame->GetYaxis()->SetLabelSize(0.1);
            hFrame->GetYaxis()->SetLabelOffset(0.01);
            hFrame->GetYaxis()->SetTitleFont(42);
            hFrame->GetYaxis()->SetTitleSize(0.1);
            hFrame->GetYaxis()->SetTitleOffset(1.0);
            hFrame->GetYaxis()->CenterTitle();
            hFrame->GetYaxis()->SetNdivisions(505);
            // TICKS Y Axis
            hFrame->GetYaxis()->SetTickLength(xFactor * 0.04 / yFactor);
            // Format for x axis
            hFrame->GetXaxis()->SetLabelFont(42);
            hFrame->GetXaxis()->SetLabelSize(0.1);
            hFrame->GetXaxis()->SetLabelOffset(0.01);
            hFrame->GetXaxis()->SetTitleFont(42);
            hFrame->GetXaxis()->SetTitleSize(0.1);
            hFrame->GetXaxis()->SetTitleOffset(1.0);
            hFrame->GetXaxis()->CenterTitle();
            hFrame->GetXaxis()->SetNdivisions(505);
            // TICKS X Axis
            hFrame->GetXaxis()->SetTickLength(yFactor * 0.06 / xFactor);
            (*(h + (Ny-1-j) * 2 + i))->Draw("same");
        }
    }
    C->cd();
    return C;
}

double URound(double value, int digite)
{
    if (digite < 1)
    {
        return (double)(long long)(value);
    }

    long iPub = (long)pow(10, digite); // get digite squre of 10
    value *= iPub;
    long long iValue = (long long)(value + 0.5);

    return ((double)(iValue) / iPub);
}

// #define __PEAKS_C_FIT_AREAS__ 1 /* fit peaks' areas */
struct gausPAR
{
    double h; //height
    double m; //mean
    double s; //sigma
    bool operator<(const gausPAR &gP) const
    { //symbol overloading
        return m < gP.m;
    }
};

int npk = 10;

Double_t fpeaks(Double_t *x, Double_t *par)
{
    //Double_t result = par[0] + par[1]*x[0];
    Double_t result = par[0] * TMath::Gaus(x[0], par[1], par[2]);

    for (Int_t p = 0; p < npk; p++)
    {
        Double_t norm = par[3 * p + 3]; // "height" or "area"
        Double_t mean = par[3 * p + 4];
        Double_t sigma = par[3 * p + 5];
#if defined(__PEAKS_C_FIT_AREAS__)
        norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif                                                 /* defined(__PEAKS_C_FIT_AREAS__) */
        result += norm * TMath::Gaus(x[0], mean, sigma);
    }
    return result;
}

TF1 *fpeaksfit(TH1 *ha, int npeaks, double sigma, gausPAR *gPar, const char *name)
{
    gStyle->SetOptFit(111);
    double res = 2.;
    double thrd = 0.01;
    npk = TMath::Abs(npeaks);
    TCanvas *cfp1 = cdC(1);

    cfp1->cd();
    TH1 *h = (TH1 *)ha->Clone("h");
    TH1 *h2 = (TH1 *)h->Clone("h2");
    h->Draw();
    double hxL = 0;
    double hxR = 0;

    hxL = h->GetXaxis()->GetXmin();
    hxR = h->GetXaxis()->GetXmax();
    double par[300] = {0};

    double *xpeaks;

    TSpectrum *s = new TSpectrum(2 * npk, res);
    int nfound = s->Search(h, sigma, "", thrd);
    npk = 0;
    printf("Found %d candidate peaks to fit\n", nfound);
    //ha->Draw();
    //return;
    TH1 *hb = s->Background(h, 20, "same");
    if (hb)
        cfp1->Update();
    //if (np <npeaks) return;

    //estimate linear background using a fitting method
    //cfp->cd(1);
    hb->Draw();
    TF1 *fbk = new TF1("fbk", "gaus", hxL, hxR);
    fbk->SetLineColor(3);
    hb->Fit("fbk", "R");
    // Loop on all found peaks. Eliminate peaks at the background level
    //par[0] = fbk->GetParameter(0);
    //par[1] = fbk->GetParameter(1);
    //par[2] = fbk->GetParameter(2);
    par[0] = 0;
    par[1] = 0;
    par[2] = 0;
    //fbk->Draw("same");
    gPad->Clear();
    h->Draw();
    //h->Rebin();
    fbk->Draw("same");
    TF1 *fp[nfound];

    xpeaks = s->GetPositionX();

    double meanarray[100];
    double sigmaarray[100];
    for (int p = 0; p < nfound; p++)
    {
        Double_t xp = xpeaks[p];
        Int_t bin = h->GetXaxis()->FindBin(xp);
        Double_t yp = h->GetBinContent(bin);
        //if (yp - TMath::Sqrt(yp) < fbk->Eval(xp))
        //    continue;
        sprintf(buff, "fp%d", p);

        fp[p] = new TF1(buff, "gaus", hxL, hxR);
        fp[p]->SetLineColor(3);
        h->Fit(fp[p], "+q", "", xp - res * sigma, xp + res * sigma);

        //fp[p]->Draw("same");
        if (fp[p]->GetParameter(1) > hxL && fp[p]->GetParameter(1) < hxR && fp[p]->GetParameter(0) > 0)
        {
            par[3 * npk + 3] = fp[p]->GetParameter(0); // "height"
            par[3 * npk + 4] = fp[p]->GetParameter(1); // "mean"
            par[3 * npk + 5] = fp[p]->GetParameter(2); // "sigma"
            meanarray[npk] = fp[p]->GetParameter(1);
            sigmaarray[npk] = fp[p]->GetParameter(2);
#if defined(__PEAKS_C_FIT_AREAS__)
            par[3 * npk + 3] *= par[3 * npk + 5] * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif                                                                            /* defined(__PEAKS_C_FIT_AREAS__) */
            npk++;
        }
    }
    double RL = TMath::MinElement(npk, meanarray) - 10 * TMath::MaxElement(npk, sigmaarray);
    double RR = TMath::MaxElement(npk, meanarray) + 10 * TMath::MaxElement(npk, sigmaarray);
    h->GetXaxis()->SetRangeUser(RL, RR);
    printf("Found %d useful peaks to fit\n", npk);

    sprintf(buff, "%s_1.png", name);
    //sprintf(buff, "%s_1.png", name.c_str());
    cfp1->SaveAs(buff);

    //printf("Now fitting: Be patient\n");
    //return fp[1];
    TCanvas *cfp2 = cdC(2);
    cfp2->cd();
    h2->Draw();
    //h2->Rebin();

    cfp2->Update();
    //TF1 *fit = new TF1("fit",fpeaks,0,1000,2+3*npeaks);
    //LZWfunc *myfunc = new LZWfunc();
    TF1 *fit = new TF1("fit", fpeaks, RL, RR, 3 + 3 * npk);
    // We may have more than the default 25 parameters
    TVirtualFitter::Fitter(h2);
    fit->SetParameters(par);
    fit->SetNpx(10000);
    //TPaveStats* stat = (TPaveStats*)h2->FindObject("stats");
    //h->SetStats(0);
    h2->GetXaxis()->SetRangeUser(RL, RR);
    h2->Fit("fit");
    //return fit;
    //gausPAR gPar[100];

    //
    // the first of three pars are bcakground
    for (int i = 0; i < npk; i++)

    {
        gPar[i].h = fit->GetParameter(3 * i + 3);
        gPar[i].m = fit->GetParameter(3 * i + 4);
        gPar[i].s = fit->GetParameter(3 * i + 5);
    }
    sort(gPar, gPar + npk);

    TPaveStats *stat = (TPaveStats *)h2->GetListOfFunctions()->FindObject("stats");
    //stat->SetOptFit(1100);
    stat->SetOptFit(1111);
    h2->GetXaxis()->SetRangeUser(gPar[0].m - 10 * gPar[0].s, gPar[0].m + 10 * gPar[0].s);
    cfp2->Update();
    sprintf(buff, "%s_2.png", name);
    cfp2->SaveAs(buff);
    //fit->Draw("same");
    return fit;
}

TF1 *gausfit(TH1 *h, double sigma, double facleft, double facright, int rbU, double UL, double UR)
{
    double mean = 0;
    //double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser(UL, UR);
    mean = hU->GetBinCenter(h->GetMaximumBin());
    //sigma = hU->GetRMS();
    TF1 *fitU = new TF1("fitU", "gaus", mean - facleft * sigma, mean + facright * sigma);
    hU->GetXaxis()->SetRangeUser(mean - facleft * sigma, mean + facright * sigma);
    cout << mean << "\t" << sigma << endl;

    // fitU->SetParLimits(0, 0,hU->GetMaximum()*1.1);
    fitU->SetParameter(1, mean);
    hU->Fit(fitU, "Q");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;
    //return NULL;
    TFitResultPtr failed = hU->Fit(fitU, "Q", "", mean - facleft * sigma, mean + facright * sigma);
    //failed =1 means fit failed
    if (failed)
        //return hU = NULL;
        return fitU = NULL;

    else
    {

        if (UL < mean - 20 * sigma)
            UL = mean - 20 * sigma;
        if (UR > mean + 20 * sigma)
            UR = mean + 20 * sigma;

        hU->GetXaxis()->SetRangeUser(UL, UR);
        //return hU;
        return fitU;
    }
}

/*
TH1 *gausfit(TH1 *h, double facleft, double facright, int rbU, double UL, double UR)
{
    double mean = 0;
    double sigma = 0;
    TH1 *hU = (TH1 *)h->Clone();
    hU->Draw();
    hU->Rebin(rbU);
    hU->GetXaxis()->SetRangeUser(UL, UR);
    TF1 *fitU = new TF1("fitU", "gaus", UL, UR);
    mean = hU->GetBinCenter(hU->GetMaximumBin());

    fitU->SetParLimits(0, 0,hU->GetMaximum()*1.1);
    fitU->SetParameter(1, mean);
    cout << mean << "\t" << sigma << endl;
    hU->Fit(fitU, "R");
    mean = fitU->GetParameter(1);
    sigma = fitU->GetParameter(2);

    cout << mean << "\t" << sigma << endl;
    //return NULL;
    TFitResultPtr failed = hU->Fit(fitU, "", "", mean - facleft * sigma, mean + facright * sigma);
    //failed =1 means fit failed
    if (failed)
        return hU = NULL;

    else
    {

        if (UL < mean - 20 * sigma)
            UL = mean - 20 * sigma;
        if (UR > mean + 20 * sigma)
            UR = mean + 20 * sigma;

        hU->GetXaxis()->SetRangeUser(UL, UR);
        return hU;
    }
}
*/
TF1 *twogausfit(TH1 *ht, double fac, double rangefac, int rbt, double tL, double tR)
{
    //First fit for ensuring the rangement of histgram;
    TH1 *h = (TH1 *)ht->Clone();
    h->Rebin(rbt);
    double mean = h->GetBinCenter(h->GetMaximumBin());
    double sigma = h->GetRMS();
    TF1 *fit = new TF1("fit", "gaus", mean - sigma, mean + sigma);
    h->GetXaxis()->SetRangeUser(mean - sigma, mean + sigma);
    fit->SetParameter(1, mean);
    //fit->SetParameter(2,sigma);
    h->Fit(fit);
    mean = fit->GetParameter(1);
    sigma = TMath::Abs(fit->GetParameter(2));
    //return h;
    if (tL < mean - 5 * sigma || sigma > 1)
    {
        tL = mean - 10 * sigma;
        tR = mean + 10 * sigma;
    }
    cout << h->GetName() << "\t" << tL << "\t" << tR << endl;

    h->GetXaxis()->SetRangeUser(tL, tR);

    TF1 *fit2 = new TF1("fit2", "gaus(0)+gaus(3)", tL, tR);
    fit2->SetParNames("C_{TR}", "#mu_{TR}", "#sigma_{TR}", "C_{bkgnd}", "#mu_{bkgnd}", "#sigma_{bkgnd}");
    fit2->SetParameter(1, mean);
    fit2->SetParameter(2, sigma);
    fit2->SetParLimits(3, 0, fit->GetParameter(0) * fac);
    fit2->SetParameter(4, mean - sigma);
    fit2->SetParLimits(4, mean - 3 * sigma, mean + 5 * sigma);
    fit2->SetParameter(5, 2 * sigma);
    fit2->SetParLimits(5, sigma, 3 * sigma);

    //h->Fit(fit2);
    h->Fit(fit2, "", "", mean - rangefac * sigma, mean + rangefac * sigma);
    TF1 *fit_tr = new TF1("fit_tr", "gaus", tL, tR);
    fit_tr->SetParameter(0, fit2->GetParameter(0));
    fit_tr->SetParameter(1, fit2->GetParameter(1));
    fit_tr->SetParameter(2, fit2->GetParameter(2));
    TF1 *fit_bg = new TF1("fit_bg", "gaus", tL, tR);
    fit_bg->SetParameter(0, fit2->GetParameter(3));
    fit_bg->SetParameter(1, fit2->GetParameter(4));
    fit_bg->SetParameter(2, fit2->GetParameter(5));
    fit_tr->SetLineColor(3);
    fit_tr->SetLineStyle(7);
    fit_bg->SetLineColor(3);
    fit_bg->SetLineStyle(7);
    fit_tr->Draw("same");
    fit_bg->Draw("same");

    if (fit2->GetParameter(3) >= fit2->GetParameter(0) && TMath::Abs(fit2->GetParameter(5)) >= 0.05)
    {
        mean = fit2->GetParameter(4);
        sigma = TMath::Abs(fit2->GetParameter(5));
    }
    else if (TMath::Abs(fit2->GetParameter(2)) >= 0.05)
    {
        mean = fit2->GetParameter(1);
        sigma = TMath::Abs(fit2->GetParameter(2));
    }

    if (tL < mean - 10 * sigma)
    {
        tL = mean - 10 * sigma;
    }
    if (tR > mean + 10 * sigma)
    {
        tR = mean + 10 * sigma;
    }

    cout << h->GetName() << "\t" << tL << "\t" << tR << endl;
    h->GetXaxis()->SetRangeUser(tL, tR);
    //return h;
    return fit_tr;
}

TF1 *profilefit(TH2 *Rt, double rbU, double rbt, double tL, double tR, double UL, double UR, char *name)
{

    TCanvas *cAT = new TCanvas("cAT", "cAT", 800, 600);
    TCanvas *cpfx = new TCanvas("cpfx", "cpfx", 800, 600);
    cAT->Clear();
    cpfx->Clear();

    cAT->cd();
    TH2 *Qt = (TH2 *)Rt->Clone();
    //TH2* Qt = (TH2*) Rt->Clone("tmp");
    Qt->Draw("colz");

    Qt->RebinX(rbU * 2);
    Qt->RebinY(rbt * 2);
    Qt->GetYaxis()->SetRangeUser(tL, tR);
    Qt->GetXaxis()->SetRangeUser(UL, UR);
    Qt->ProfileX();

    cpfx->cd();
    TH1 *Qpfx = Qt->ProfileX();
    //Qpfx->Reset();

    //Qpfx=(TH1*)gDirectory->Get("Qt_pfx");
    Qpfx->Draw();
    Qpfx->GetYaxis()->SetRangeUser(tL, tR);
    Qpfx->GetYaxis()->SetTitle("Timediff (ps)");
    Qpfx->GetXaxis()->SetRangeUser(UL, UR);
    Qpfx->GetXaxis()->SetTitle("TOT (ps)");

    //TF1 *fitQt = new TF1("fitQt", "[0]+[1]/TMath::Sqrt(abs(x))+[2]/abs(x)+[3]/abs(x)/TMath::Sqrt(abs(x))+[4]/abs(x)/abs(x)", UL, UR);
    TF1 *fitQt = new TF1("fitQt", "pol5", UL, UR);

    fitQt->SetNpx(1000000);
    Qpfx->Fit(fitQt, "R");

    sprintf(buff, "%s_pfx.png", name);
    cpfx->SaveAs(buff);
    sprintf(buff, "%s_ATrelationship.png", name);
    cAT->SaveAs(buff);
    //sfile->WriteTObject(Qpfx);
    //sfile->WriteTObject(Qt);

    Qpfx->Reset();
    //delete Qt;
    //delete Qpfx;
    //delete c5;
    return fitQt;
}
TGraph *DrawMyGraph(float *x, float *y, const int N)
{
    //TGraph();
    TGraph *g1 = new TGraph(N, x, y);
    TCanvas *c1 = cdC(0);
    g1->Draw("AP");
    return g1;
    //TGraph* g2 = new TGraph(n,x1,y1cor);
}

#endif // #ifdef MyClass_cxx