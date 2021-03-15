#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TRandom.h>
#include "TString.h"
#include <time.h>
#include "DrawMyClass.h"
TRandom3 r;
TH1D *hr = new TH1D("hr", "hr", 2000, -0.1e-9, 0.1e-9);
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
bool SearchforID(int id, vector<int> &vector)
{
    for (int i = 0; i < vector.size(); i++)
    {
        if (vector.at(i) == id)
            return true;
    }
    return false;
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
    double Trecept = 500e-12; //waiting the photons hit
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
        TFitResultPtr p3 = g->Fit(fitU, "QR");
        //cout<<"p3=\t"<<p3<<endl;
        if (p3)
        {
            TF1 *fit2 = new TF1("fit2", "pol2", U_RL, U_RR);
            TFitResultPtr p2 = g->Fit(fit2, "QR");
            //cout<<"p2=\t"<<p2<<endl;

            if (p2)
            {
                TF1 *fit1 = new TF1("fit1", "pol1", U_RL, U_RR);
                TFitResultPtr p1 = g->Fit(fit1, "QR");
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

    void GetTrackerAngle(int N, double *TrackerX, double *TrackerY, double *TrackerZ, double *theta, double *phi, double targetX, double *targetY, double *targetZ)
    {
        double *p[3];
        p[0] = TrackerX;
        p[1] = TrackerY;
        p[2] = TrackerZ;

        double exp[2] = {0};
        double delta[3] = {0};
        TGraphErrors *g;

        TVector3 v1;
        TVector3 vXaxis(1, 0, 0);
        for (int i = 0; i < 2; i++)
        {

            g = new TGraphErrors();
            for (int j = 0; j < N; j++)
            {
                //cout<<"x"<<p[0][j]<<",yz"<<i<<","<<p[1+i][j]<<endl;
                g->SetPoint(j, p[0][j], p[1 + i][j]);
            }
            // fit
            g->Fit("pol1", "q");
            g->Draw();
            exp[i] = g->GetFunction("pol1")->Eval(targetX);
            //delta[i + 1] = p[1 + i][0] - p[1 + i][N - 1];
            delta[1 + i] = g->GetFunction("pol1")->Eval(p[0][N - 1]) - g->GetFunction("pol1")->Eval(p[0][0]);
        }
        *targetY = exp[0];
        *targetZ = exp[1];
        // x =(1,0,0), line=()
        delta[0] = p[0][N - 1] - p[0][0];
        v1.SetX(delta[0]);
        v1.SetY(delta[1]);
        v1.SetZ(delta[2]);
        //*theta = v1.Angle(vXaxis);
        *theta = TMath::ACos(-1 * v1.x() / v1.Mag());
        *phi = TMath::ACos(v1.z() / TMath::Sqrt(v1.z() * v1.z() + v1.y() * v1.y()));
        if (v1.y() < 0)
            *phi = -1 * (*phi);
        //cout<<"theta = "<< *theta<<", phi"<<*phi<<endl;
        /*
           cout<<"theta 1 = "<< *theta<<endl;
         *theta = TMath::ACos(delta[0]/TMath::Sqrt(delta[1]*delta[1]+delta[2]*delta[2]+delta[0]*delta[0])); //Unit: rad
         cout<<"theta 2 = "<< *theta<<endl;
         *theta = TMath::ACos(v1.x()/v1.Mag());
         cout<<"theta 3 ="<< *theta<<endl;

         cout<<"Phi"<< *phi<<endl;
         cout<<v1.z()<<endl;
         cout<<v1.y()<<endl;
         cout<<"z="<<v1.Mag()*sin(*theta)*cos(*phi)<<endl;
         cout<<"y="<<v1.Mag()*sin(*theta)*sin(*phi)<<endl;
         */
    }
    double remainder(double A, int r)
    {
        int n = int(A) / r;
        return A - n * r;
    }
    double RebuildTrack(double np, double Inx, double Iny, double Inz, double Px, double Py, double Pz, int ID = 1)
    {
        /*
           double Ox, Oy, Oz;
           Oz = A1z;
           Ox = (A1z - Inz) / pz * px + Inx;
        //cout<<Ox<<endl;
        Oy = (A1z - Inz) / pz * py + Iny;
        return sqrt(Ox * Ox + Oy * Oy + Oz * Oz);
        //cout << "A1z: " << A1z<<"\t"<<Inz<<"\t"<<pz<<"\t"<<px<<"\t"<<Inx << endl;
        */
        //cout << "Input: " << Inx << "\t" << Iny << "\t" << Inz << "\t" << Px << "\t" << Py << "\t" << Pz << endl;
        double Ax0, Ay0, Az0;
        double px, py, pz;
        double swap;
        if (ID == 1 || ID == 3)
        {
            Az0 = -1 * (ID - 2) * 90;
            Ax0 = 0;
            Ay0 = 0;
            px = Px;
            py = Py;
            pz = -1 * (ID - 2) * Pz;
        }
        if (ID == 2 || ID == 4)
        {
            Az0 = -1 * (ID - 3) * 90;
            Ax0 = 0;
            Ay0 = 0;
            px = Px;
            py = Pz;
            pz = -1 * (ID - 3) * Py;
            swap = Iny;
            Iny = Inz;
            Inz = swap;
        }
        double thetaC = TMath::ACos(1 / np);
        double theta = TMath::ACos(pz); //angle between momentum direction of mu and z axis;

        double Ax, Ay, Az;
        Az = Az0;
        double dx, dy, dz;
        dz = Az0 - Inz;
        double Rdxdy;
        Rdxdy = dz * TMath::Tan(theta - thetaC);
        dx = Rdxdy / sqrt(px * px + py * py) * px;
        dy = Rdxdy / sqrt(px * px + py * py) * py;
        Ax = Inx + dx;
        Ay = Iny + dy;
        return sqrt(Ax * Ax + Ay * Ay + Az * Az);

        //return -1;
    }

    // why the result is only correct in the case of ID==3
    double RebuildTOP(double np, double Inx, double Iny, double Inz, double Px, double Py, double Pz, int ID = 1)
    {
        //cout << "Input: " << Inx << "\t" << Iny << "\t" << Inz << "\t" << px << "\t" << py << "\t" << pz << endl;
        double A1x, A1y, A1z;
        double px, py, pz;
        if (ID == 1 || ID == 3)
        {
            A1z = -1 * (ID - 2) * 90;
            A1x = 0;
            A1y = 0;
            px = Px;
            py = Py;
            pz = Pz;
        }
        if (ID == 2 || ID == 4)
        {
            A1z = -1 * (ID - 3) * 90;
            A1x = 0;
            A1y = 0;
            px = Px;
            py = Pz;
            pz = Py;
        }

        double TOP;
        double v = 298 / np; //velocity of photon,Unit:mm/ns

        double dx, dy;
        double dz = A1z - Inz;
        //double Amaxx, Amaxy;
        double theta = TMath::ACos(-1 * pz); //angle between momentum direction of mu and z axis;
        double thetaC = TMath::ACos(1 / np);

        double Ax, Ay, Az;

        double Rdxdy = dz * TMath::Tan(theta - thetaC);
        //cout << "dz= " << dz << endl;
        //cout << "theta= " << theta << endl;
        //cout << "thetaC= " << thetaC << endl;
        //cout << "Rdxdy = " << Rdxdy << endl;
        if (Rdxdy > 600)
        {
            cout << "No  solution! Rdxdy > 600 " << endl;
            return -1;
        }
        else if (TMath::Abs(theta - 0.5 * TMath::Pi()) < 0.1)
        {
            Az = A1z;
            double Rdxdy = dz * TMath::Tan(theta - thetaC);
            dx = Rdxdy / sqrt(px * px + py * py) * px;
            dy = Rdxdy / sqrt(px * px + py * py) * py;
            Ax = Inx + dx;
            Ay = Iny + dy;
            //cout << "Ax,Ay,Az = " << Ax << "\t" << Ay << "\t" << Az << endl;
            if (remainder(Ax, 10) * remainder(Ax, 10) + Ay * Ay < 5 * 5)
            {

                return TOP = sqrt(Ax * Ax + Ay * Ay + Az * Az) / v;
            }
            else
            {
                cout << "No  solution! No hit within PMT sensitive area" << endl;
                return TOP = sqrt(Ax * Ax + Ay * Ay + Az * Az) / v;
                //return -1;
            }
        }
        else
        {

            double Ox, Oy, Oz;
            Oz = A1z;
            //cout << "A1z: " << A1z<<"\t"<<Inz<<"\t"<<pz<<"\t"<<px<<"\t"<<Inx << endl;
            Ox = (A1z - Inz) / pz * px + Inx;
            //cout<<Ox<<endl;
            Oy = (A1z - Inz) / pz * py + Iny;

            double Rxy = sqrt(Ox * Ox + Oy * Oy); //r of out point of mu in xy-plane
            double a = (dz) / TMath::Cos(theta - thetaC) * TMath::Sin(thetaC);
            double b = a / TMath::Cos(theta);
            //Aminx = (sqrt(Ox*Ox+Oy*Oy)-b)/sqrt(Ox*Ox+Oy*Oy)*Ox;
            //Aminy = (sqrt(Ox*Ox+Oy*Oy)-b)/sqrt(Ox*Ox+Oy*Oy)*Oy;
            //Amaxx = (sqrt(Ox*Ox+Oy*Oy)+b)/sqrt(Ox*Ox+Oy*Oy)*Ox;
            //Amaxy = (sqrt(Ox*Ox+Oy*Oy)+b)/sqrt(Ox*Ox+Oy*Oy)*Oy;
            double Oxprim = remainder(Ox, 10);
            double n = int(Ox / 10); //the number of reflection
            //cout << "A1z: " << A1z<< endl;
            //cout << "Input: " << Inx << "\t" << Iny << "\t" << Inz << "\t" << px << "\t" << py << "\t" << pz << endl;
            //cout << "Rxy = " <<Rxy << endl;
            //cout << "Ox,Oy,Oz = " << Ox << "\t" << Oy << "\t" << Oz << endl;
            if ((-Ox) * (-Ox) / a / a + (-Oy) * (-Oy) / b / b <= 1)
            {
                cout << "The solution exists!" << endl;
                Ax = 0 + int(Ox / 10) * 10;
                Ay = 0 + int(Oy / 10) * 10;
                //cout << "Ax,Ay,Az = " << Ax << "\t" << Ay << "\t" << Az << endl;
                return TOP = sqrt(Ax * Ax + Ay * Ay + Az * Az) / v;
            }
            else
            {

                if (remainder(Ox, 10) * remainder(Ox, 10) + Oy * Oy < (5 + b) * (5 + b))
                {
                    cout << "No  solution! No hit within PMT sensitive area" << endl;
                    return TOP = sqrt(Ox * Ox + Oy * Oy + Oz * Oz) / v;
                    //return -1;
                }
                else
                {

                    if (Oxprim < 0)
                        Ax = -1 * sqrt(5 * 5 * Oxprim * Oxprim / (Oy * Oy + Oxprim * Oxprim));
                    if (Oxprim > 0)
                        Ax = +1 * sqrt(5 * 5 * Oxprim * Oxprim / (Oy * Oy + Oxprim * Oxprim));
                    Ay = Oxprim / Oxprim * Oy;
                    if ((Ax - Oxprim) * (Ax - Oxprim) / a / a + (Ay - Oy) * (Ay - Oy) / b / b <= 1)
                    {
                        cout << "The solution exists!" << endl;
                        Ax = Ax + 10 * n;
                        //cout << "Ax,Ay,Az = " << Ax << "\t" << Ay << "\t" << Az << endl;
                        return TOP = sqrt(Ax * Ax + Ay * Ay + Az * Az) / v;
                    }
                    else
                    {
                        cout << "Uncertainty situation" << endl;
                        return TOP = sqrt(Ox * Ox + Oy * Oy + Oz * Oz) / v;
                        //return -1;
                    }
                }
            }

            /*
               if (remainder(Ox, 10) * remainder(Ox, 10) + Oy * Oy > (5 + b) * (5 + b))
               {
               cout << "Ox,Oy,Oz = " << Ox << "\t" << Oy << "\t" << Oz << endl;
               if (Rxy - b < 0)
               {
               Ax = 0 + int(Ox / 10) * 10;
               Ay = 0 + int(Oy / 10) * 10;
               }
               if (Rxy - b > 0)
               {
               Ax = (Rxy - b) / Rxy * Ox;
               Ay = (Rxy - b) / Rxy * Oy;
               }
               return TOP = (Ax * Ax + Ay * Ay + Az * Az) / v;
               }
               else
               {
               cout << "No  solution! " << endl;
               return -1;
               }
               */
        }
    }
    void GetFileList(TString filePath, TString filePattern, vector<TString> &fList)
    {
        char line[1000];
        fList.clear();
        if (filePath.EndsWith("/"))
            filePath.Remove(filePath.Length() - 1, 1);
        FILE *fp = gSystem->OpenPipe("ls " + filePath + "/*" + filePattern, "r");
        if (!fp)
        {
            cout << "----> NO data files exists in " << filePath << "!" << endl;
            return;
        }

        while (fgets(line, sizeof(line), fp))
        {
            TString s(line);
            if (s.Index(filePattern) == -1)
                continue;
            fList.push_back(s.ReplaceAll("\n", ""));
        }
    };
    TString GetFilename(TString fileName)
    {
        cout << "The fileName is:" << fileName << endl;
            TString prename;
                prename=fileName.Copy().Remove(0, fileName.Last('/')+1);
                    prename = prename.Remove(prename.Last('.'),prename.Length());
                        cout << "Dat filename is:" << prename << endl;
                            return prename;
                            };
    //char path[1000] = "/Users/liziwei/learning/CRTest/build/angleop";
    //char path[1000] = "/Users/liziwei/learning/CRTest/build";
    char path[1000] = "/mnt/c/Subsys/work/CRTest/build";
    void CalculateEff(TString filepath = "../build")
    {
        clock_t start, finish;
        double totaltime;
        start = clock();

        char buff[1024];

        const int T0N = 1; //the number of T0
        const int R3809_sensorN = 4;
        const int R3809_sensor_ch_N = 1;
        const int T0_sensorN = T0N * R3809_sensorN * R3809_sensor_ch_N; //the number of channels of T0
        const int FTOFN = 1;                                            //the number of FTOF
        const int R10754_sensorN = 8;
        const int R10754_sensor_ch_N = 16;
        const int FTOF_sensorN = FTOFN * R10754_sensorN * R10754_sensor_ch_N; //the number of channels of FTOF
        const int TrackerN = 4;                                               //the number of Trackers
        const int TriggerN = 2;                                               //the number of Triggers
        TH1I *hTrackereff = new TH1I("hTrackereff", "", TrackerN + 1, 0, TrackerN + 1);
        TH1I *hMMeff = new TH1I("hMMeff", "", TrackerN, 0, TrackerN);
        TH1I *hlackMMeff = new TH1I("hlackMMeff", "", TrackerN, 0, TrackerN);
        TH1I *hT0eff = new TH1I("hT0eff", "", T0_sensorN + 1, 0, T0_sensorN + 1);
        TH1I *hT0CHeff = new TH1I("hT0CHeff", "", T0_sensorN, 0, T0_sensorN);
        TH1I *hFTOFeff = new TH1I("hFTOFeff", "", FTOF_sensorN + 1, 0, FTOF_sensorN + 1);
        TH1I *hFTOFCHeff = new TH1I("hFTOFCHeff", "", FTOF_sensorN, 0, FTOF_sensorN);
        TChain *t1 = new TChain("Run");
        vector<TString> rootlist;
        GetFileList(filepath,"root",rootlist);
        TString rootname=GetFilename(rootlist.at(0));
        for(int i=0;i<rootlist.size();i++){
        t1->Add(rootlist.at(i));
        cout<<"Add rootfile to My chain: "<<rootlist.at(i)<<endl;
        }
        //return;

        vector<int> *T0chID = new vector<int>;
        vector<int> *FTOFchID = new vector<int>;
        vector<double> *mu_hittime = new vector<double>;
        vector<double> *mu_hitx = new vector<double>; //the position of incident event
        vector<double> *mu_hity = new vector<double>;
        vector<double> *mu_hitz = new vector<double>;
        vector<int> *mu_hitID = new vector<int>;
        vector<double> *mu_px = new vector<double>; //the direction of incident event
        vector<double> *mu_py = new vector<double>;
        vector<double> *mu_pz = new vector<double>;
        vector<double> *mu_energy = new vector<double>;
        //TFile *f1 = new TFile(buff, "READ");
        //TTree *t1 = (TTree *)f1->Get("Run");

        t1->SetBranchAddress("R107.id", &FTOFchID);
        t1->SetBranchAddress("R380.id", &T0chID);

        t1->SetBranchAddress("mu.t", &mu_hittime);
        t1->SetBranchAddress("mu.x", &mu_hitx);
        t1->SetBranchAddress("mu.y", &mu_hity);
        t1->SetBranchAddress("mu.z", &mu_hitz);
        t1->SetBranchAddress("mu.DetID", &mu_hitID);

        t1->SetBranchAddress("mu.px", &mu_px);
        t1->SetBranchAddress("mu.py", &mu_py);
        t1->SetBranchAddress("mu.pz", &mu_pz);
        t1->SetBranchAddress("mu.E", &mu_energy);

        int triggercounter = 0;
        int trackercounter = 0;
        int T0counter = 0;
        int FTOFcounter = 0;

        vector<int> triggervec;
        vector<int> trackervec;
        vector<int> T0vec;
        vector<int> FTOFvec;
        int N = t1->GetEntries();
        cout << "Entries = " << N << endl;
        for (int iEvent = 0; iEvent < N; iEvent++)
            //for (int iEvent = 0; iEvent < 10; iEvent++)
        {
            if (iEvent % 10000 == 0)
                cout << "The Entry No: " << iEvent << endl;

            //-----------initial----------------------//
            T0chID->clear();
            FTOFchID->clear();
            mu_hittime->clear();
            mu_hitx->clear();
            mu_hity->clear();
            mu_hitz->clear();
            mu_hitID->clear();
            mu_px->clear();
            mu_py->clear();
            mu_pz->clear();
            mu_energy->clear();

            triggervec.clear();
            trackervec.clear();
            T0vec.clear();
            FTOFvec.clear();

            t1->GetEntry(iEvent);
            // create vector to store id

            for (int ihitdet = 0; ihitdet < mu_hitID->size(); ihitdet++)
            {
                int theID = mu_hitID->at(ihitdet);
                if (theID == 300 || theID == 301)
                    triggervec.push_back(theID);
                if (theID >= 0 && theID < TrackerN)
                    {
                        trackervec.push_back(theID);
                        hMMeff->Fill(theID);
                    }
                if (theID == 100)
                {
                    for (int iT0hit = 0; iT0hit < T0chID->size(); iT0hit++)
                    {
                        int theT0chID = T0chID->at(iT0hit);
                        if (T0vec.empty())
                            T0vec.push_back(theT0chID);
                        else if (!SearchforID(theT0chID, T0vec))
                            T0vec.push_back(theT0chID);
                    }
                }
                if (theID == 200)
                {
                    for (int iFTOFhit = 0; iFTOFhit < FTOFchID->size(); iFTOFhit++)
                    {
                        int theFTOFchID = FTOFchID->at(iFTOFhit);
                        if (FTOFvec.empty())
                            FTOFvec.push_back(theFTOFchID);
                        else if (!SearchforID(theFTOFchID, FTOFvec))
                            FTOFvec.push_back(theFTOFchID);
                    }
                }
            }
            if (triggervec.size() > 1)
            {

                triggercounter++;
                if (trackervec.size() > 0)
                {

                    trackercounter++;
                    hTrackereff->Fill(trackervec.size());
                    if (trackervec.size() == 3)
                        for (int j = 0; j < TrackerN; j++)
                        {

                            if (!SearchforID(j, trackervec))
                                hlackMMeff->Fill(j);
                        }
                }
                if (T0vec.size() > 0)
                {

                    T0counter++;
                    hT0eff->Fill(T0vec.size());
                    for (int j = 0; j < T0vec.size(); j++)
                    {
                        hT0CHeff->Fill(T0vec.at(j));
                    }
                }
                if (FTOFvec.size() > 0)
                {

                    FTOFcounter++;
                    hFTOFeff->Fill(FTOFvec.size());
                    for (int j = 0; j < FTOFvec.size(); j++)
                    {
                        hFTOFCHeff->Fill(FTOFvec.at(j));
                    }
                }
            }
        }

        cout << "Trigger coincidence: " << triggercounter << endl
            << "tracker counter: " << trackercounter << endl
            << "T0 counter: " << T0counter << endl
            << "FTOF counter: " << FTOFcounter << endl;
        TLatex *la;
        TCanvas *c;
        int ccounter=0;


        //** draw tracker eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hTrackereff, "Numbers of MM", "Counts",2, 3);
        hTrackereff->Draw();
        double Trackereff[TrackerN];   
        for (int i = 0; i < TrackerN; i++)
        {

            Trackereff[i] = hTrackereff->GetBinContent(i + 2) / trackercounter;
            sprintf(buff, "Numbers of MM=%d,Eff=%.1f%%", i + 1, Trackereff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        TString canvasname;
        canvasname = filepath + "/"+rootname+"Trackereff.png";
        //sprintf(buff,"%s/%sTrackereff.png",filepath,rootname);
        c->SaveAs(canvasname);

        //** draw MM eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hMMeff, "ID of MM", "Counts",2, 3);
        hMMeff->Draw();
        double MMeff[TrackerN];   
        for (int i = 0; i < TrackerN; i++)
        {

            MMeff[i] = hMMeff->GetBinContent(i + 1) / trackercounter;
            sprintf(buff, "ID of MM=%d,Eff=%.1f%%", i + 1, MMeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
       // sprintf(buff,"%s/%sMMeff.png",filepath,rootname);
        canvasname = filepath+"/"+rootname+"MMeff.png";
        c->SaveAs(canvasname);

        //** draw MM true eff
        //*
        int allfires=allfires = hTrackereff->GetBinContent(TrackerN + 1);
        c = cdC(ccounter++);
        DrawMyHist(hlackMMeff, "ID of MM", "Counts",2, 3);
        hlackMMeff->Draw();
        double MMtrueeff[TrackerN];   
        for (int i = 0; i < TrackerN; i++)
        {

            MMtrueeff[i] = allfires/(allfires+hlackMMeff->GetBinContent(i + 1));
            sprintf(buff, "ID of MM=%d,true Eff=%.1f%%", i + 1, MMtrueeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        //sprintf(buff,"%s/%sMMtrueeff.png",filepath,rootname);
        canvasname = filepath+"/"+rootname+"MMtrueeff.png";
        c->SaveAs(canvasname);

        //** draw T0 eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hT0eff, "Numbers of firedch", "Counts",2, 3);
        hT0eff->Draw();
        double T0eff[T0_sensorN];   
        for (int i = 0; i < T0_sensorN; i++)
        {

            T0eff[i] = hT0eff->GetBinContent(i + 2) / T0counter;
            sprintf(buff, "Numbers of firedch=%d,Eff=%.1f%%", i + 1, T0eff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        canvasname = filepath+"/"+rootname+"T0eff.png";
        c->SaveAs(canvasname);

        //** draw T0ch eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hT0CHeff, "ID of T0", "Counts",2, 3);
        hT0CHeff->Draw();
        double T0CHeff[T0_sensorN];   
        for (int i = 0; i < T0_sensorN; i++)
        {

            T0CHeff[i] = hT0CHeff->GetBinContent(i + 1) / T0counter;
            sprintf(buff, "ID of T0=%d,Eff=%.1f%%", i + 1, T0CHeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        canvasname = filepath+"/"+rootname+"T0CHeff.png";
        c->SaveAs(canvasname);

        //** draw FTOF eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hFTOFeff, "Numbers of firedch", "Counts",2, 3);
        hFTOFeff->Draw();
        double FTOFeff[FTOF_sensorN];   
        for (int i = 0; i < FTOF_sensorN; i++)
        {

            FTOFeff[i] = hFTOFeff->GetBinContent(i + 2) / FTOFcounter;
            sprintf(buff, "Numbers of firedch=%d,Eff=%.1f%%", i + 1, FTOFeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        canvasname = filepath+"/"+rootname+"FTOFeff.png";
        c->SaveAs(canvasname);

        //** draw FTOFch eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hFTOFCHeff, "ID of FTOF", "Counts",2, 3);
        hFTOFCHeff->Draw();
        double FTOFCHeff[FTOF_sensorN];   
        for (int i = 0; i < FTOF_sensorN; i++)
        {

            FTOFCHeff[i] = hFTOFCHeff->GetBinContent(i + 1) / FTOFcounter;
            sprintf(buff, "ID of FTOF=%d,Eff=%.1f%%", i + 1, FTOFCHeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        canvasname = filepath+"/"+rootname+"FTOFCHeff.png";
        c->SaveAs(canvasname);

    }
    void CalculateTR(const char *rootname = "CRY3", double fac = 0.2, const char *ParType = "CFD", unsigned long processN = 1)
    {
        //void Outputfun_MCP(const char *rootname="",double fac = -30, const char* ParType="FIX"){

        cout << "fac=" << fac << ",dicriminate:" << ParType << endl;

        /*===========================================
         * ============Procedure timing start========
         * =========================================*/
        clock_t start, finish;
        double totaltime;
        start = clock();

        /*
           TLegend *DrawMyLeg(Double_t xlow = 0.2, Double_t ylow = 0.2, Double_t xup = 0.5, Double_t yup = 0.5, Int_t textFont = 62, Size_t textSize = 0.05);

           TLatex *DrawMyLatex(char *text, Double_t x = 0.65, Double_t y = 0.5, Int_t textFont = 62, Size_t textSize = 0.05, Color_t colorIndex = 2);

           void DrawMyGraph(TGraph * datagraph, char *xtitle, char *ytitle, Size_t MSize = 1, Style_t MStyle = 28, Color_t MColor = 1, Color_t LColor = 1, Width_t LWidth = 1, Style_t LStyle = 1, Color_t FColor = 16);

           void DrawMyHist1(TH1 * datahist, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Color_t TitleColor = 1);

           void SetMyPad(TVirtualPad * pad, float left, float right, float top, float bottom);

           void DrawMyPad(TVirtualPad * pad, const char *xname, const char *yname, float x1, float x2, float y1, float y2);

        //void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor=1, Width_t LWidth=3, Color_t TitleColor=1);
        */
        //gStyle->SetOptStat(0);
        //gStyle->SetOptFit(0);
        gStyle->SetOptTitle(0);

        //TRandom3 r;
        //r.SetSeed(0);
        char name[1024];
        char buff[1024];

        const int PMTN = 4; //the number of PMT copyvolume
        const int DetN = 1;
        const int T = PMTN * DetN; //the number of PMT copyvolume
        const int TrackerN = 2;    //the number of Trackers

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

        double ttssigma = 20e-12;
        const int range = 400; // 25ps/sample
        Double_t thrd = -30;   //Umax = -28.94mV
        double Rate = 0;

        double possigma = 3000e-3; //50 um

        bool flag = 0;
        int index = 0;
        double keypoint = 0;
        double T0_1stpe[T] = {0}; // hit time of first pe
        double xT0[T] = {0};      //time stamp of waveform
        double TOP[T] = {0};      // TOP of photon
        int CHID[T] = {0};
        int NPE[T] = {0};
        double U[T] = {0}; //Amplitude of waveform

        int DetID[DetN] = {0};    //incident time
        double DetT0[DetN] = {0}; //incident time
        double InX[DetN] = {0};   //incident position
        double InY[DetN] = {0};
        double InZ[DetN] = {0};

        double InE[DetN] = {0};

        double InPX; //incident direction
        double InPY;
        double InPZ;

        double Intheta; //incident direction
        double Inphi;

        double TrackerX[TrackerN] = {0}; //Tracker impact pos
        double TrackerY[TrackerN] = {0}; //Tracker impact pos
        double TrackerZ[TrackerN] = {0}; //Tracker impact pos
        double Caltheta;                 //incident direction
        double Calphi;
        double CalX[DetN] = {0};
        double CalY[DetN] = {0};
        double CalZ[DetN] = {0};

        const int certain = 2;

        vector<vector<double>> vecTOP(T);
        vector<vector<double>> par(T);
        vector<vector<double>> tts(T);

        vector<double> *hitT = new vector<double>;
        vector<double> *TOPT = new vector<double>;
        vector<int> *ID = new vector<int>;

        vector<double> *DethitT = new vector<double>;

        vector<double> *IncidX = new vector<double>; //the position of incident event
        vector<double> *IncidY = new vector<double>;
        vector<double> *IncidZ = new vector<double>;
        vector<double> *IncidPX = new vector<double>; //the direction of incident event
        vector<double> *IncidPY = new vector<double>;
        vector<double> *IncidPZ = new vector<double>;

        vector<double> *IncidE = new vector<double>;
        vector<int> *IncidID = new vector<int>;

        //count = new vector<int>;
        int N = 0, hitN = 0, pN = 0;
        int effN = 0; // effective Events;
        double temp = 0.;
        Double_t x[range] = {};
        Double_t y[T][range] = {};

        sprintf(name, "%s/%s", path, rootname);
        sprintf(buff, "%s.root", name);

        TFile *f1 = new TFile(buff, "READ");
        TTree *t1 = (TTree *)f1->Get("Run");

        t1->SetBranchAddress("PmtS.t", &hitT);
        t1->SetBranchAddress("PmtS.TOP", &TOPT);
        t1->SetBranchAddress("PmtS.id", &ID);

        t1->SetBranchAddress("mu.t", &DethitT);
        t1->SetBranchAddress("mu.x", &IncidX);
        t1->SetBranchAddress("mu.y", &IncidY);
        t1->SetBranchAddress("mu.z", &IncidZ);
        t1->SetBranchAddress("mu.px", &IncidPX);
        t1->SetBranchAddress("mu.py", &IncidPY);
        t1->SetBranchAddress("mu.pz", &IncidPZ);

        t1->SetBranchAddress("mu.E", &IncidE);
        t1->SetBranchAddress("mu.DetID", &IncidID);

        //sprintf(name,"Thrd_%g",abs(thrd));

        sprintf(buff, "%sdata.root", name);

        TFile *f2 = new TFile(buff, "RECREATE");
        TTree *t2 = new TTree("data", "restore analysed data  from G4");
        sprintf(buff, "CHID[%d]/I", T); // CHID
        t2->Branch("CHID", &CHID, buff);
        sprintf(buff, "NPE[%d]/I", T); // NPE
        t2->Branch("NPE", &NPE, buff);
        sprintf(buff, "U[%d]/D", T);
        t2->Branch("U", &U, buff);     //-mV
        sprintf(buff, "xT0[%d]/D", T); // ns
        t2->Branch("xT0", &xT0, buff);
        sprintf(buff, "T0_1stpe[%d]/D", T); //ns
        t2->Branch("T0_1stpe", &T0_1stpe, buff);
        sprintf(buff, "TOP[%d]/D", T); //ns
        t2->Branch("TOP", &TOP, buff);

        sprintf(buff, "DetID[%d]/I", DetN); // DetectorID
        t2->Branch("DetID", &DetID, buff);
        sprintf(buff, "DetT0[%d]/D", DetN); // ns
        t2->Branch("DetT0", &DetT0, buff);

        sprintf(buff, "InX[%d]/D", DetN);
        t2->Branch("InX", &InX, buff); //mm
        sprintf(buff, "InY[%d]/D", DetN);
        t2->Branch("InY", &InY, buff); //mm
        sprintf(buff, "InZ[%d]/D", DetN);
        t2->Branch("InZ", &InZ, buff); //mm
        sprintf(buff, "InE[%d]/D", DetN);
        t2->Branch("InE", &InE, buff); //mm

        sprintf(buff, "InPX[%d]/D", 1);
        t2->Branch("InPX", &InPX, buff); //mm
        sprintf(buff, "InPY[%d]/D", 1);
        t2->Branch("InPY", &InPY, buff); //mm
        sprintf(buff, "InPZ[%d]/D", 1);
        t2->Branch("InPZ", &InPZ, buff); //mm

        sprintf(buff, "Intheta[%d]/D", 1);
        t2->Branch("Intheta", &Intheta, buff); //mm
        sprintf(buff, "Inphi[%d]/D", 1);
        t2->Branch("Inphi", &Inphi, buff); //mm

        sprintf(buff, "CalX[%d]/D", DetN);
        t2->Branch("CalX", &CalX, buff); //mm
        sprintf(buff, "CalY[%d]/D", DetN);
        t2->Branch("CalY", &CalY, buff); //mm
        sprintf(buff, "CalZ[%d]/D", DetN);
        t2->Branch("CalZ", &CalZ, buff); //mm

        sprintf(buff, "Caltheta[%d]/D", 1);
        t2->Branch("Caltheta", &Caltheta, buff); //mm
        sprintf(buff, "Calphi[%d]/D", 1);
        t2->Branch("Calphi", &Calphi, buff); //mm
        /*
           sprintf(buff, "InPX[%d]/D", T);
           t2->Branch("InPX", &InPX, buff); //mm

           sprintf(buff, "InPY[%d]/D", T);
           t2->Branch("InPY", &InPY, buff); //mm
           sprintf(buff, "InPZ[%d]/D", T);
           t2->Branch("InPZ", &InPZ, buff); //mm
           */

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
            h[iT] = new TH1D(buff, buff, binNum, RL * 1e9, RR * 1e9);
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
        //N=1;
        for (int iEvent = 0; iEvent < N; iEvent++)
            //for (int iEvent = 0; iEvent < 10; iEvent++)
        {
            if (iEvent % 10000 == 0)
                cout << "The Entry No: " << iEvent << endl;

            //-----------initial----------------------//
            hitT->clear();
            TOPT->clear();
            ID->clear();
            DethitT->clear();
            IncidX->clear();
            IncidY->clear();
            IncidZ->clear();
            IncidPX->clear();
            IncidPY->clear();
            IncidPZ->clear();
            IncidE->clear();
            IncidID->clear();
            vector<vector<double>>(T).swap(vecTOP);
            vector<vector<double>>(T).swap(par);
            vector<vector<double>>(T).swap(tts);

            //for(int i = certain; i < certain+1; i++){
            //par[i]=r.Gaus(2.4e-9,0.5e-9);
            //par[i]=4e-9;
            t1->GetEntry(iEvent);
            if (!IncidX->empty() && !hitT->empty())
            {
                int Trackers = 0;

                // Pick up which layer fires
                for (int iLayer = 0; iLayer < IncidX->size(); iLayer++)
                {
                    if (IncidID->at(iLayer) == 0 || IncidID->at(iLayer) == 3)
                        Trackers++;
                }
                if (Trackers >= 2)

                    //if (IncidX->size() == TrackerN + DetN)
                {
                    effN++;

                    hitN = hitT->size();
                    pN = IncidX->size();
                    //cout << "hitN = " << hitN << endl;
                    //cout << "pN = " << pN << endl;
                    //myFun = new TF1("myFun",outputfunc,RL,RR,temp);

                    //Record the impact position and direction of muon
                    int Trackercounter = 0;
                    int Detcounter = 0;
                    //cout<<"counter = "<< pN <<endl;

                    for (int i = 0; i < pN; i++)
                    {
                        if (IncidID->at(i) >= 100)
                        {

                            DetT0[Detcounter] = (*DethitT)[i];
                            DetID[Detcounter] = IncidID->at(i);
                            InE[Detcounter] = (*IncidE)[i];
                            InX[Detcounter] = (*IncidX)[i];
                            InY[Detcounter] = (*IncidY)[i];
                            InZ[Detcounter] = (*IncidZ)[i];
                            if (IncidID->at(i) == 100)
                            {

                                InPX = (*IncidPX)[i];
                                InPY = (*IncidPY)[i];
                                InPZ = (*IncidPZ)[i];
                                Intheta = TMath::ACos(-1 * InPX);

                                Inphi = TMath::ACos(InPZ / TMath::Sqrt(InPY * InPY + InPZ * InPZ));
                                if (InPY < 0)
                                    Inphi = -1 * Inphi;
                            }
                            //Inphi[Detcounter] = TMath::ATan(InPY[Detcounter] / InPZ[Detcounter]);
                            //cout << "SIMU Y ,Z " << InY[Detcounter] << ", " << InZ[Detcounter] << endl;
                            //cout << "theta SIMU = " << Intheta << endl;
                            //cout << "Phi SIMU = " << Inphi << endl;

                            Detcounter++;
                        }
                        else if (IncidID->at(i) == 1 || IncidID->at(i) == 2)
                        {
                            TrackerX[Trackercounter] = IncidX->at(i);
                            TrackerY[Trackercounter] = IncidY->at(i) + r.Gaus(0, possigma);
                            TrackerZ[Trackercounter] = IncidZ->at(i) + r.Gaus(0, possigma);
                            Trackercounter++;
                        }
                    }
                    for (int i = 0; i < Detcounter; i++)
                    {

                        CalX[i] = InX[i];
                        GetTrackerAngle(Trackercounter, TrackerX, TrackerY, TrackerZ, &Caltheta, &Calphi, CalX[i], &CalY[i], &CalZ[i]);
                    }
                    for (int iT = 0; iT < T; iT++)
                    {
                        //return;
                        CHID[iT] = iT;
                        for (int k = 0; k < hitN; k++)
                        {
                            //cout<< T[][k] <<endl;
                            temp = (*hitT)[k] * 1e-9;
                            if ((*ID)[k] == iT)
                            {
                                //r.SetSeed(time(NULL)*processN+k);
                                vecTOP[iT].push_back((*TOPT)[k]);
                                par[iT].push_back(temp);
                                tts[iT].push_back(r.Gaus(0, ttssigma));
                                if (iEvent == N - 1)
                                    h[iT]->Fill(temp * 1e9);
                            }
                        }
                        NPE[iT] = par[iT].size();
                        if (!vecTOP[iT].empty())
                            TOP[iT] = TMath::MinElement(NPE[iT], &vecTOP[iT][0]);
                        //parR[k]=8.3e-9;
                        //cout<<" [+] par "<<k<<"\t"<<parR.at(k)<<endl;
                        //cout<<"par"<<k<<" = "<<par[k]<<endl;
                        sort(par[iT].begin(), par[iT].end());
                        //cout<<">>> progress check <<<"<<endl;
                        if (!par[iT].empty())
                            T0_1stpe[iT] = par[iT].at(0) * 1e9;
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
                            y[iT][j] = outputfunc(x[j], par[iT], tts[iT]);
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
            }
        }
        f1->Close();
        if (effN < 1)
        {
            cout << " No effective Event, Check Your data please !" << endl;

            return;
        }
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

        ymin = TMath::MinElement(T, U) * 1.2;
        ymax = -0.2 * ymin;
        DrawMyPad(gPad, "Time (ns)", "Amp (mV)", RL * 1e9, RR * 1e9, ymin, ymax);

        TLegend *leg;
        leg = DrawMyLeg(0.6, 0.2, 0.8, 0.45);
        //return;
        for (int iT = 0; iT < T; iT++)
        {
            DrawMyGraph(g[iT], "Time (ns)", "Amp (mV)", 1, 28, 1, iT + 1, 2, 1);
            g[iT]->Draw("L");
            sprintf(buff, "PMT%d", iT);
            leg->AddEntry(g[iT], buff, "l");
        }
        leg->Draw();

        c->cd(2);
        SetMyPad(gPad, 0.12, 0.05, 0.05, 0.12);
        TLatex *l[T];
        for (int iT = 0; iT < T; iT++)
        {

            DrawMyHist(h[iT], "Time (ns)", "Counts", iT + 1, 2);
            if (iT == 0)
                h[iT]->Draw();
            else
                h[iT]->Draw("SAMES");
            sprintf(buff, "pmt%d hitN=%g", iT, h[iT]->GetSum());
            l[iT] = DrawMyLatex(buff, 0.65, 0.6 - 0.1 * iT, 62, 0.05, iT + 1);
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

        //hSig[0]->FillRandom("gaus",1000);
        //hSig[0]->Draw();
        //return ;

        TCanvas *c1 = new TCanvas("c1", "", 800, 600);
        SetMyPad(gPad, 0.12, 0.05, 0.05, 0.12);

        c1->cd();
        for (int iT = 0; iT < T; iT++)
        {

            DrawMyHist(hSig[iT], "Time (ns)", "Counts", iT + 1, 2);
            if (iT == 0)
                hSig[iT]->Draw();
            else
                hSig[iT]->Draw("SAMES");
            hSig[iT]->Rebin(1);
            hSig[iT]->Fit(fSig[iT], "QR");

            sprintf(buff, "pmt%d TR=%0.2f", iT, fSig[iT]->GetParameter(2));
            l[iT] = DrawMyLatex(buff, 0.65, 0.6 - 0.1 * iT, 62, 0.05, iT + 1);
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
        vector<double>().swap(*TOPT);
        vector<double>().swap(*DethitT);
        vector<int>().swap(*ID);
        vector<double>().swap(*IncidX);
        vector<double>().swap(*IncidY);
        vector<double>().swap(*IncidZ);
        vector<double>().swap(*IncidPX);
        vector<double>().swap(*IncidPY);
        vector<double>().swap(*IncidPZ);
        vector<double>().swap(*IncidE);
        vector<int>().swap(*IncidID);

        delete hitT;
        delete TOPT;
        delete DethitT;
        delete ID;
        delete IncidID;
        delete IncidX;
        delete IncidY;
        delete IncidZ;
        delete IncidPX;
        delete IncidPY;
        delete IncidPZ;
        delete IncidE;

        /*=======================================================*
         * ================Procedure timing end==================*
         * ======================================================*/
        finish = clock();
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "\nthe whole time through the procedure is " << totaltime << "s!!" << endl; //delete count;
    }

    void GetTimeRes(const char *rootname = "CRY3data")
    {
        //void MygStyle();
        //MygStyle();
        gStyle->SetOptFit(111);

        sprintf(buff, "%s/%s.dat", path, rootname);
        ofstream out(buff);
        cout << "===> Create your date file: " << buff << endl;

        char name[100];
        char buff[1024];

        //***************************************************//
        //--------------Configuration-----------------------//
        //***************************************************//
        double tL = -2.;
        double tR = 2.;
        double AL = 0;
        double AR = 300;

        double uL = -1e3;
        double uR = 0;

        int rbA = 10;
        int rbt = 4;
        int rbu = 1;

        // the range set After Correct
        //double L2 = -100e-12;
        //double R2 = 100e-12;
        int bint = (tR - tL) / 0.5e-3;
        int binu = (uR - uL) / 30;
        const int PMTN = 4; //the number of PMT copyvolume
        const int DetN = 1;
        const int T = PMTN * DetN; //the number of PMT copyvolume

        //---------------------------------------------------//
        //***************************************************//

        //double RL,RR;
        //double U_RL,U_RR;
        //int binT = (init_tR-init_tL)/0.5e-12;
        //int binU = (init_UR-init_UL)/1;

        cout << "Start =====>>>>" << endl;
        //return;
        //cout<<"<<---- Succeed excuating ---->>"<<endl;
        //return;
        //thrd = (s+1)*0.2;
        sprintf(name, "%s/%s", path, rootname);

        sprintf(buff, "%s.root", name);
        TFile *f = new TFile(buff, "READ");
        TTree *t = (TTree *)f->Get("data");
        cout << "Read rootfile: " << buff << endl;

        double T0_1stpe[T] = {0}; // hit time of first pe
        double xT0[T] = {0};      //time stamp of waveform
        double TOP[T] = {0};      // TOP of photon
        int CHID[T] = {0};
        int NPE[T] = {0};
        double U[T] = {0}; //Amplitude of waveform

        int DetID[DetN] = {0};    //incident time
        double DetT0[DetN] = {0}; //incident time
        double InX[DetN] = {0};   //incident position
        double InY[DetN] = {0};
        double InZ[DetN] = {0};

        double InE[DetN] = {0};

        double InPX; //incident direction
        double InPY;
        double InPZ;

        double Intheta; //incident direction
        double Inphi;

        double Caltheta; //incident direction
        double Calphi;
        double CalX[DetN] = {0};
        double CalY[DetN] = {0};
        double CalZ[DetN] = {0};

        t->SetBranchAddress("CHID", CHID);
        t->SetBranchAddress("NPE", NPE);
        t->SetBranchAddress("U", U); //-mV
        t->SetBranchAddress("xT0", xT0);
        t->SetBranchAddress("T0_1stpe", T0_1stpe);
        t->SetBranchAddress("TOP", TOP);

        t->SetBranchAddress("DetID", DetID);
        t->SetBranchAddress("DetT0", DetT0);

        t->SetBranchAddress("InE", InE);    //mm
        t->SetBranchAddress("InX", InX);    //mm
        t->SetBranchAddress("InY", InY);    //mm
        t->SetBranchAddress("InZ", InZ);    //mm
        t->SetBranchAddress("InPX", &InPX); //momentum
        t->SetBranchAddress("InPY", &InPY); //momentum
        t->SetBranchAddress("InPZ", &InPZ); //momentum
        t->SetBranchAddress("Intheta", &Intheta);
        t->SetBranchAddress("Inphi", &Inphi);

        t->SetBranchAddress("CalX", CalX); //mm
        t->SetBranchAddress("CalY", CalY); //mm
        t->SetBranchAddress("CalZ", CalZ); //mm
        t->SetBranchAddress("Caltheta", &Caltheta);
        t->SetBranchAddress("Calphi", &Calphi);

        TH1D *ht = new TH1D("ht", ";Time (ns);Counts", bint, tL, tR);
        TH1D *htcor = new TH1D("htcor", ";Time (ns);Counts", bint, tL, tR);
        TH1D *htheta = new TH1D("htheta", ";#theta (rad);Counts", 200, 0, 1.5);
        TH1D *hphi = new TH1D("hphi", ";#phi (rad);Counts", 200, -3.15, 3.15);
        TH1D *hdtheta = new TH1D("hdtheta", ";#Delta#theta (rad);Counts", 2000, -1, 1);
        TH1D *hdphi = new TH1D("hdphi", ";#Delta#phi (rad);Counts", 2000, -1, 1);
        TH1D *hNPMT = new TH1D("hNPMT", ";Number of Fired PMT;Counts", 5, 0, 5);
        TH1D *hPMTID = new TH1D("hPMTID", ";PMTID;Counts", 4, 0, 4);

        TH2D *hpos = new TH2D("hpos", ";Y (mm); Z (mm)", 200, -95, 95, 200, -95, 95);
        TH1D *hdy = new TH1D("hdy", ";#Delta#Y (mm);Counts", 2000, -10, 10);
        TH1D *hdz = new TH1D("hdy", ";#Delta#Z (mm);Counts", 2000, -10, 10);

        TH2D *hNPE = new TH2D("hNPE", ";PMTID; NPE", 8, 0, 8, 16, -1, 15);
        TH2D *h2dPMT = new TH2D("h2dPMT", ";Number of Fired PMT;TR (ns)", 5, 0, 5, bint, tL, tR);

        TH1D *hA = new TH1D("hA", ";reTrack (mm);Counts", 3e3, AL, AR);

        TH2D *hAT = new TH2D("hAT", ";reTrack;TR (ns)", 3e3, AL, AR, bint, tL, tR);

        TH1D *htheta_cut = (TH1D *)htheta->Clone("htheta_cut");
        TH1D *hphi_cut = (TH1D *)hphi->Clone("hphi_cut");
        TH2D *hpos_cut = (TH2D *)hpos->Clone("hpos_cut");

        int N = t->GetEntries();
        cout << "Total trigger events is: " << N << endl;
        double Timestamp = 0;
        double Timestampcor = 0;
        double TimeSum = 0;
        double PMTtime = 0;
        double PMTtimecor = 0;
        double reTOP = 0;
        double reTrack = 0;
        double reTrackSum = 0;
        double meanreTrack = 0;
        double CalPx, CalPy, CalPz;
        vector<double> TT;
        vector<double> AA;
        TT.reserve(20e4);
        AA.reserve(20e4);
        //for (int i = 0; i < N; i++)
        for (int i = 0; i < N; i++)
        {
            t->GetEntry(i);
            hdtheta->Fill(Caltheta - Intheta);
            hdphi->Fill(Calphi - Inphi);
            htheta->Fill(Caltheta);
            hphi->Fill(Calphi);
            hpos->Fill(InY[0], InZ[0]);
            hdy->Fill(InY[0] - CalY[0]);
            hdz->Fill(InZ[0] - CalZ[0]);
            CalPx = -1 * TMath::Cos(Caltheta);
            CalPz = TMath::Sin(Caltheta) * TMath::Cos(Calphi);
            CalPy = TMath::Sin(Caltheta) * TMath::Sin(Calphi);
            //cout<<CalPx<<"\t"<<CalPy<<"\t"<<CalPz<<"\t"<<endl;
            //cout<<InPX<<"\t"<<InPY<<"\t"<<InPZ<<"\t"<<endl;

            for (int s = 0; s < T; s++)
            {

                hNPE->Fill(CHID[s], NPE[s]);
            }
            Timestamp = 0;
            Timestampcor = 0;
            PMTtime = 0;
            PMTtimecor = 0;
            reTrack = 0;
            reTrackSum = 0;
            //cout<<"U & xT0: "<<U[0]<<"\t"<<xT0[0]<<endl;
            int PMTcounter = 0;
            /*
               for (int k = 0; k < PMTN; k++)
               {
               if (U[PMTN + k] < 0 && xT0[PMTN + k] != 0 && U[k] < 0 && xT0[k] != 0)
               {
               PMTtime += xT0[PMTN + k] - xT0[k];
               PMTcounter++;
               }
               }
               */
            for (int k = 0; k < PMTN; k++)
            {
                //int k = 3;
                if (U[k] < 0 && xT0[k] > 0)
                {
                    //reTOP = RebuildTOP(1.5, InX[0], InY[0], InZ[0], InPX, InPY, InPZ, k + 1);
                    reTrack = RebuildTrack(1.5, CalX[0], CalY[0], CalZ[0], CalPx, CalPy, CalPz, k + 1);
                    reTrackSum += reTrack;
                    //cout << "delta TOP" << reTrack / 298 * 1.5 - TOP[0] << endl;

                    PMTtime += xT0[k] - DetT0[0];
                    PMTtimecor += xT0[k] - DetT0[0] - reTrack / 298 * 1.5;
                    //PMTtime = xT0[0] - DetT0[0];

                    PMTcounter++;
                    hPMTID->Fill(CHID[k]);
                }
            }
            //if (PMTcounter >= 1)    return;
            /*
               for (int j = 0; j < DetN; j++)
               {
               int PMTcounter = 0;
               DetTime[j] = 0;
               TimeSum = 0;
               for (int k = 0; k < PMTN; k++)
               {
            //hNPE->Fill(CHID[j * PMTN + k],NPE[j * PMTN + k]);
            //if (1)
            if (U[j * PMTN + k] < 0 && xT0[j * PMTN + k] != 0)
            {

            TimeSum += xT0[j * PMTN + k];
            PMTcounter++;
            //cout << "DetTime " << PMTcounter << "\t," << TimeSum << endl;
            }
            }
            //if (InY[0]*InY[0]+InZ[0]*InZ[0]<5*5)
            //if (PMTcounter == 1)
            if (1)
            {
            DetTime[j] = xT0[j * PMTN + 3];
            //DetTime[j] = TimeSum / PMTcounter;
            //cout << "DetTime & PMTcounter: " << PMTcounter << "& " << DetTime[j] << endl;
            }
            }
            */
            //FlyTime = (InX[0] - InX[1]) / cos(Caltheta) / 1.e3 / 3.e8 * 1.e9;
            //Timestamp = PMTtime / PMTcounter - FlyTime;
            Timestamp = PMTtime / PMTcounter;
            Timestampcor = PMTtimecor / PMTcounter;
            meanreTrack = reTrackSum / PMTcounter;
            hNPMT->Fill(PMTcounter);
            h2dPMT->Fill(PMTcounter, Timestamp);

            //cout<<"FlyTime"<<FlyTime<<endl;
            if (PMTcounter == 4)
            {
                TT.push_back(Timestampcor);
                //AA.push_back(reTrack);
                AA.push_back(meanreTrack);

                //cout << "DetTime[1]" << DetTime[1] << "\t" << DetTime[0] << endl;
                //cout << "timestamp: " << Timestamp << endl;
                //cout << "Timestamp: " << Timestamp << endl;
                //ht->Fill(Timestamp);
                htcor->Fill(Timestamp);
                htheta_cut->Fill(Caltheta);
                hphi_cut->Fill(Calphi);
                hpos_cut->Fill(InY[0], InZ[0]);
            }
            /*
               if (DetTime[1] > 0 && DetTime[0] > 0)
               {

               Timestamp = DetTime[1] - DetTime[0] + FlyTime;
            //cout << "DetTime[1]" << DetTime[1] << "\t" << DetTime[0] << endl;
            //cout << "timestamp: " << Timestamp << endl;
            ht->Fill(Timestamp);
            htheta_cut->Fill(Caltheta[0]);
            hphi_cut->Fill(Calphi[0]);
            hpos_cut->Fill(InY[0],InZ[0]);
            }
            */
        }
        TCanvas *c;
        int CNum = 0;
        TF1 *fit;
        TLegend *leg;
        TLatex *la;

        char pngprefix[100];
        sprintf(pngprefix, "%s/%s", path, rootname);

        //
        // ---------draw efficiency of PMT--------//
        //
        double PMTEff[4];
        c = cdC(CNum++);
        DrawMyHist(hPMTID, "", "", 1, 2);
        //ht->Rebin(2);
        hPMTID->Draw();
        hPMTID->GetXaxis()->SetNdivisions(505);
        for (int i = 0; i < 4; i++)
        {

            PMTEff[i] = hPMTID->GetBinContent(hPMTID->FindBin(i)) / N;
            sprintf(buff, "NPMT=%d,Eff=%.1f%%", i, PMTEff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.4 + 0.1 * i, 42, 0.06);
            la->Draw("same");
        }
        //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sPMTID.png", pngprefix);
        c->SaveAs(buff);
        //return;

        //
        // ---------draw Time Res --------//
        //
        c = cdC(CNum++);
        DrawMyHist(ht, "", "", 1, 2);
        //ht->Rebin(2);
        htcor->Draw();
        //return;
        //fit = gausfit(htcor, 200e-3, 3, 3, rbt*8, tL, tR);
        fit = gausfit(htcor, 200e-3, 3, 3, rbt * 2, tL, tR);
        //fitT = gausfit(ht, 0.1, 3, 3, rbt*2, tL, tR);
        if (!fit)
        {
            cout << "the htfit hist is NULL " << endl;
            return;
        }
        sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sTimeRes.png", pngprefix);
        c->SaveAs(buff);
        //return;

        //
        // ---------draw track vs TR --------//
        //
        c = cdC(CNum++);
        const int iter = 5; //the number of iteration

        //ht->Reset();
        TH1F *htfit;
        TH1F *hAfit;
        TF1 *fitT;
        TF1 *fitA;
        TF1 *fitAT = new TF1("fitAT", "0", AL, AR);
        double Treserve[200000];
        double TAcor = 0;
        double Amean = 0;
        double Asigma = 0;
        double Tmean = 0;
        double Tsigma = 0;
        for (int s = 0; s < iter; s++)
        {
            if (s != 0)
            {
                if (Amean - 3 * Asigma <= 0)

                    //fitAT = profilefit(hAT, rbA *2, rbt * 10, Tmean - 3 * Tsigma, Tmean + 3 * Tsigma, 0.1, Amean + 3 * Asigma, buff);
                    fitAT = profilefit(hAT, rbA, rbt * 8, Tmean - 3 * Tsigma, Tmean + 3 * Tsigma, 0.1, Amean + 3 * Asigma, buff);
                else
                    //fitAT = profilefit(hAT, rbA * 2, rbt * 10, Tmean - 3 * Tsigma, Tmean + 3 * Tsigma, Amean - 3 * Asigma, Amean + 3 * Asigma, buff);
                    fitAT = profilefit(hAT, rbA, rbt * 8, Tmean - 3 * Tsigma, Tmean + 3 * Tsigma, Amean - 3 * Asigma, Amean + 3 * Asigma, buff);
                if (!fitAT)
                {
                    cout << " the profilefit is failed! " << endl;
                    return;
                }
                ht->Reset();
                hAT->Reset();
            }
            for (int i = 0; i < TT.size(); i++)
            {
                if (s == 0)
                {
                    cout << "AA & TT: " << AA[i] << "\t" << TT[i] << endl;
                    hA->Fill(AA[i]);
                    Treserve[i] = TT[i];
                    //tL = 23.5;
                    //tR = 25;
                    tL = -2;
                    tR = 2;
                    //tL = 43.5;
                    //tR = 46;
                }
                else
                {

                    TAcor = fitAT->Eval(AA[i]);
                    Treserve[i] = Treserve[i] - TAcor;
                    tL = -1;
                    tR = 1;
                }
                ht->Fill(Treserve[i]);
                hAT->Fill(AA[i], Treserve[i]);
            }
            c->cd();
            if (s == 0)
            {
                fitA = gausfit(hA, 10, 3, 3, rbA, AL, AR);
                //return;
                if (!fitA)
                {
                    cout << "the hAfit hist is NULL " << endl;
                    return;
                }
                DrawMyHist(hA, "", "", 1, 3);
                //fitA = (TF1 *)hAfit->GetFunction("fitU");
                Amean = fitA->GetParameter(1);
                Asigma = fitA->GetParameter(2);
                cout << "Amean=" << Amean << ",\tAsigma=" << Asigma << endl;
                sprintf(buff, "%s_reTrack.png", pngprefix);
                c->SaveAs(buff);
            }
            c->Clear();
            //fitT = gausfit(ht, 0.1, 3, 3, rbt * 8, tL, tR);
            fitT = gausfit(ht, 0.1, 3, 3, rbt * 4, tL, tR);
            if (!fitT)
            {
                cout << "the htfit hist is NULL " << endl;
                return;
            }
            DrawMyHist(ht, "", "", 1, 3);
            //fitT = (TF1 *)htfit->GetFunction("fitU");
            Tmean = fitT->GetParameter(1);
            Tsigma = fitT->GetParameter(2);

            ht->SetNdivisions(505);
            sprintf(buff, "#sigma=%.0fps", Tsigma * 1e3);
            TLatex *l = DrawMyLatex(buff, 0.2, 0.5);
            l->Draw();
            cout << "Tmean=" << Tmean << ",\tTsigma=" << Tsigma << endl;
            sprintf(buff, "%s_TR_cor%d.png", pngprefix, s);
            c->SaveAs(buff);
        }
        c = cdC(CNum++);
        DrawMy2dHist(hAT, "", "", 1, 2);
        //ht->Rebin(2);
        hAT->Draw("colz");
        //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%shAT.png", pngprefix);
        c->SaveAs(buff);
        //return;

        //
        // ---------draw Number of fired PMT--------//
        //
        double Eff[5];
        c = cdC(CNum++);
        DrawMyHist(hNPMT, "", "", 1, 2);
        //ht->Rebin(2);
        hNPMT->Draw();
        hNPMT->GetXaxis()->SetNdivisions(505);
        for (int i = 0; i < 5; i++)
        {

            Eff[i] = hNPMT->GetBinContent(hNPMT->FindBin(i)) / hNPMT->Integral();
            sprintf(buff, "NPMT=%d,Eff=%.1f%%", i, Eff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.4 + 0.1 * i, 42, 0.06);
            la->Draw("same");
        }
        //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sNPMT.png", pngprefix);
        c->SaveAs(buff);
        //return;

        //
        // ---------draw NPMT vs TR --------//
        //
        c = cdC(CNum++);
        DrawMy2dHist(h2dPMT, "", "", 1, 2);
        //ht->Rebin(2);
        h2dPMT->Draw("colz");
        h2dPMT->GetXaxis()->SetNdivisions(505);
        //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%s2dPMT.png", pngprefix);
        c->SaveAs(buff);
        //return;

        //
        // ---------draw delta theta --------//
        //
        c = cdC(CNum++);
        DrawMyHist(hdtheta, "", "", 1, 2);
        //ht->Rebin(2);
        hdtheta->GetXaxis()->SetNdivisions(505);
        hdtheta->Draw();
        fit = gausfit(hdtheta, 30e-3, 3, 3, 1, -0.1, 0.1);
        sprintf(buff, "#sigma_{#theta}=%.0fmrad", fit->GetParameter(2) * 1e3);
        la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sdeltatheta.png", pngprefix);
        c->SaveAs(buff);

        //
        // ---------draw delta phi --------//
        //
        c = cdC(CNum++);
        DrawMyHist(hdphi, "", "", 1, 2);
        //ht->Rebin(2);
        hdphi->GetXaxis()->SetNdivisions(505);
        hdphi->Draw();
        fit = gausfit(hdphi, 100e-3, 3, 3, 2, -0.2, 0.2);
        sprintf(buff, "#sigma_{#phi}=%.0fmrad", fit->GetParameter(2) * 1e3);
        la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sdeltaphi.png", pngprefix);
        c->SaveAs(buff);

        //
        // ---------draw NPE --------//
        //
        c = cdC(CNum++);
        SetMyPad(gPad, 0.15, 0.15, 0.05, 0.15);
        DrawMy2dHist(hNPE, "", "");
        hNPE->Draw("colz");
        sprintf(buff, "%sNPE.png", pngprefix);
        c->SaveAs(buff);

        //
        // ---------draw theta --------//
        //
        c = cdC(CNum++);
        leg = DrawMyLeg(0.6, 0.7, 0.8, 0.9);
        DrawMyHist(htheta, "", "", 1, 2);
        htheta->Draw();
        leg->AddEntry(htheta, "No cut", "l");
        DrawMyHist(htheta_cut, "", "", 2, 2);
        htheta_cut->Draw("SAME");
        leg->AddEntry(htheta_cut, "All NPE>0", "l");
        leg->Draw();
        sprintf(buff, "%stheta.png", pngprefix);
        c->SaveAs(buff);
        //
        // ---------draw phi --------//
        //
        c = cdC(CNum++);
        DrawMyHist(hphi, "", "", 1, 2);
        hphi->Draw();
        //hphi->GetYaxis()->SetRangeUser(0, 15e3);
        DrawMyHist(hphi_cut, "", "", 2, 2);
        //hphi_cut->Draw("");
        hphi_cut->Draw("SAME");
        leg->Draw();
        sprintf(buff, "%sphi.png", pngprefix);
        c->SaveAs(buff);

        //
        // ---------draw hit pos--------//
        //
        leg = DrawMyLeg(0.6, 0.85, 0.9, 0.9);
        leg->SetNColumns(2);
        c = cdC(CNum++, 800, 800);
        DrawMy2dHist(hpos, "", "", 24, 1, 1);
        hpos->Draw();
        leg->AddEntry(hpos, "No cut", "l");

        DrawMy2dHist(hpos_cut, "", "", 3, 2, 1);
        hpos_cut->Draw("SAME");
        leg->AddEntry(hpos_cut, "All NPE>0", "p");
        leg->Draw();
        sprintf(buff, "E=%.2f", hpos_cut->Integral() / hpos->Integral());
        la = DrawMyLatex(buff, 0.2, 0.85);
        sprintf(buff, "%shitpos.png", pngprefix);
        c->SaveAs(buff);

        //
        // ---------draw delta y --------//
        //
        c = cdC(CNum++);
        DrawMyHist(hdy, "", "", 1, 2);
        hdy->Rebin(20);
        hdy->GetXaxis()->SetNdivisions(505);
        hdy->Draw();
        fit = gausfit(hdy, 3, 3, 3, 1, -10, 10);
        sprintf(buff, "#sigma_{y}=%.0fmm", fit->GetParameter(2));
        la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sdeltay.png", pngprefix);
        c->SaveAs(buff);

        //
        // ---------draw delta z --------//
        //
        c = cdC(CNum++);
        DrawMyHist(hdz, "", "", 1, 2);
        hdz->Rebin(20);
        hdz->GetXaxis()->SetNdivisions(505);
        hdz->Draw();
        fit = gausfit(hdz, 3, 3, 3, 1, 10, 10);
        sprintf(buff, "#sigma_{z}=%.0fmm", fit->GetParameter(2));
        la = DrawMyLatex(buff, 0.2, 0.4);
        sprintf(buff, "%sdeltaz.png", pngprefix);
        c->SaveAs(buff);
    }

    void PrintResults(const char *rootname = "x2y2z1data")
    {

        //void MygStyle();
        //MygStyle();
        //gStyle->SetOptFit(111);

        sprintf(buff, "%s/%s.dat", path, rootname);
        ofstream out(buff);
        cout << "===> Create your date file: " << buff << endl;

        char name[100];
        char buff[1024];

        //***************************************************//
        //--------------Configuration-----------------------//
        //***************************************************//
        double tL = 0.;
        double tR = 3.;

        double uL = -1e3;
        double uR = 0;

        int rbt = 2;
        int rbu = 1;

        // the range set After Correct
        //double L2 = -100e-12;
        //double R2 = 100e-12;
        int bint = (tR - tL) / 5e-3;
        int binu = (uR - uL) / 30;
        const int T = 4;
        const int iter = 5; //the number of iteration
        //---------------------------------------------------//
        //***************************************************//

        //double RL,RR;
        //double U_RL,U_RR;
        //int binT = (init_tR-init_tL)/0.5e-12;
        //int binU = (init_UR-init_UL)/1;

        cout << "Start =====>>>>" << endl;
        //return;
        //cout<<"<<---- Succeed excuating ---->>"<<endl;
        //return;
        //thrd = (s+1)*0.2;
        sprintf(name, "%s/%s", path, rootname);

        sprintf(buff, "%s.root", name);
        TFile *f = new TFile(buff, "READ");
        TTree *t = (TTree *)f->Get("data");

        double xT0[T] = {0}; //time stamp of waveform
        int CHID[T] = {0};
        int NPE[T] = {0};
        double U[T] = {0}; //Amplitude of waveform
        double InE[T] = {0};
        double InX[T] = {0}; //incident position
        double InY[T] = {0};
        double InZ[T] = {0};
        double Intheta[T] = {0};
        double Inphi[T] = {0};

        double CalX[T] = {0}; //incident position
        double CalY[T] = {0};
        double CalZ[T] = {0};
        double Caltheta[T] = {0};
        double Calphi[T] = {0};

        //double InPX[T] = {0}; //incident direction
        //double InPY[T] = {0};
        //double InPZ[T] = {0};

        double T0_1stpe[T] = {0}; // hit time of first pe

        double Utemp = 0, Ttemp = 0;
        double T0_cor = 0;
        float TT[1000000] = {0.};

        t->SetBranchAddress("CHID", CHID);
        t->SetBranchAddress("NPE", NPE);
        t->SetBranchAddress("U", U); //-mV
        t->SetBranchAddress("xT0", xT0);
        t->SetBranchAddress("T0_1stpe", T0_1stpe);

        t->SetBranchAddress("InE", InE); //mm
        t->SetBranchAddress("InX", InX); //mm
        t->SetBranchAddress("InY", InY); //mm
        t->SetBranchAddress("InZ", InZ); //mm
        t->SetBranchAddress("Intheta", Intheta);
        t->SetBranchAddress("Inphi", Inphi);

        t->SetBranchAddress("CalX", CalX); //mm
        t->SetBranchAddress("CalY", CalY); //mm
        t->SetBranchAddress("CalZ", CalZ); //mm
        t->SetBranchAddress("Caltheta", Caltheta);
        t->SetBranchAddress("Calphi", Calphi);

        //t->SetBranchAddress("InPX", InPX); //mm
        //t->SetBranchAddress("InPY", InPY); //mm
        //t->SetBranchAddress("InPZ", InPZ); //mm

        TH1D *ht[4];
        TH1D *h1st[4];
        TH1D *hNPE[4];
        TH1D *hU[4];
        for (int iT = 0; iT < T; iT++)
        {
            sprintf(buff, "ht%d", iT);
            ht[iT] = new TH1D(buff, ";Time (ns);Counts", bint, tL, tR);
            sprintf(buff, "h1st%d", iT);
            h1st[iT] = new TH1D(buff, ";Time (ns);Counts", bint, tL, tR);
            sprintf(buff, "hU%d", iT);
            hU[iT] = new TH1D(buff, ";Amplitude (mV);Counts", binu, uL, uR);
            sprintf(buff, "hNPE%d", iT);
            hNPE[iT] = new TH1D(buff, ";NPE;Counts", 20, 0, 30);
        }
        //t->Draw("T0>>h");

        int N = 0;
        N = t->GetEntries();
        cout << "the Entries is :" << N << endl;
        for (int i = 0; i < N; i++)
        {
            t->GetEntry(i);

            for (int iT = 0; iT < T; iT++)
            {

                if (U[iT] < 0 && xT0[iT] != 0)
                {

                    ht[iT]->Fill(xT0[iT]);
                    h1st[iT]->Fill(T0_1stpe[iT]);
                    hU[iT]->Fill(U[iT]);
                    hNPE[iT]->Fill(NPE[iT]);
                }
            }
        }
        int CCounter = 0;
        TLatex *la;
        TF1 *ff;
        TCanvas *c1;
        double ANPE[T] = {0};
        double ANPEErr[T] = {0};
        double TRSigma[T] = {0};
        double TRSigmaErr[T] = {0};
        double TRMean[T] = {0};
        double TRMeanErr[T] = {0};

        for (int i = 0; i < T; i++)
        {
            // --- Draw NPE distribution ---
            c1 = cdC(CCounter++);
            DrawMyHist(hNPE[i], "Number of Photons", "Counts");
            hNPE[i]->Draw();
            ANPE[i] = -1;
            ANPEErr[i] = -999;
            ANPE[i] = hNPE[i]->GetMean();
            ANPEErr[i] = hNPE[i]->GetRMS();
            /*
               ff = gausfit(hNPE[i], 5, 2, 2, 1, 0, 30);
               if (ff)
               {
               ANPE[i] = ff->GetParameter(1);
               ANPEErr[i] = ff->GetParError(1);
               }
               */
            sprintf(buff, "Npe=%.0f", ANPE[i]);
            la = DrawMyLatex(buff, 0.55, 0.25);
            sprintf(buff, "%s/%sCH%dPhotonNum.png", path, rootname, i);
            c1->SaveAs(buff);

            // ---Time distribution ---//
            c1 = cdC(CCounter++);
            DrawMyHist(ht[i], "TimeRes (ns)", "Counts");
            ht[i]->Draw();
            ff = gausfit(ht[i], 20e-3, 2, 2, 1, 0., 3);
            TRSigma[i] = -999;
            TRSigmaErr[i] = -999;
            if (ff)
            {
                TRMean[i] = ff->GetParameter(1) * 1e3;
                TRMeanErr[i] = ff->GetParError(1) * 1e3;
                TRSigma[i] = ff->GetParameter(2) * 1e3;
                TRSigmaErr[i] = ff->GetParError(2) * 1e3;
            }
            sprintf(buff, "Sigma=%.0f#pm%.0fps", TRSigma[i], TRSigmaErr[i]);
            la = DrawMyLatex(buff, 0.55, 0.25);
            sprintf(buff, "Mean=%.0f#pm%.0fps", TRMean[i], TRMeanErr[i]);
            la = DrawMyLatex(buff, 0.55, 0.35);

            sprintf(buff, "%s/%sCH%dTR.png", path, rootname, i);
            c1->SaveAs(buff);
        }
        out << InX[0] << "\t"
            << InY[0] << "\t"
            << InZ[0] << "\t"
            << Intheta[0] << "\t"
            << Inphi[0] << endl;
        for (int i = 0; i < T; i++)
        {

            out << ANPE[i] << "\t"
                << ANPEErr[i] << "\t"
                << TRSigma[i] << "\t"
                << TRSigmaErr[i] << endl;
        }
        return;
        //c1 = cdC(CCounter++);
        //RL = h->GetMean();
        //cout<<"RL = "<<RL<<endl;
        //return;
        //c1->Divide(2,2);

        // DrawMyCanvas((TH1**)&ht[0],2, 2);
        //TH1F *hFrame = (TH1F*) (*(&ht[0]+1))->Clone();
        //hFrame->Draw();
        //return;
        /*
           for (int i = 0; i < 4; i++)
           {
           c1 = cdC(CCounter++);
           ht[i]->Draw();
           DrawMyHist(ht[i], "", "");
           sprintf(buff, "Mean=%.0fps", ht[i]->GetMean()*1e3);
           DrawMyLatex(buff, 0.6, 0.5,42,0.06,2);
           sprintf(buff, "RMS=%.0fps", ht[i]->GetRMS()*1e3);
           DrawMyLatex(buff, 0.6, 0.4,42,0.06,2);
           sprintf(buff, "%s_Time%d.png", name,i);
           c1->SaveAs(buff);
           }
           */
        c1 = cdC(CCounter++);
        c1->Divide(2, 2, 0.01, 0.01);
        for (int i = 0; i < T; i++)
        {
            c1->cd(i + 1);
            hNPE[i]->Draw();
            DrawMyHist(hNPE[i], "", "");
            hNPE[i]->GetXaxis()->SetNdivisions(505);
            gPad->SetLeftMargin(0.2);
            gPad->SetRightMargin(0.01);
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.2);
            sprintf(buff, "CH%d", i + 1);
            DrawMyLatex(buff, 0.5, 0.9, 42, 0.06, 2);
            sprintf(buff, "Mean=%.0f PE", hNPE[i]->GetMean());
            DrawMyLatex(buff, 0.6, 0.5, 42, 0.06, 2);
        }
        sprintf(buff, "%s_NPE.png", name);
        c1->SaveAs(buff);

        c1 = cdC(CCounter++);
        c1->Divide(2, 2, 0.01, 0.01);
        for (int i = 0; i < 4; i++)
        {
            c1->cd(i + 1);
            ht[i]->Draw();
            DrawMyHist(ht[i], "", "");
            ff = gausfit(ht[i], 20e-3, 2, 2, 1, 0., 3);
            TRSigma[i] = -999;
            TRSigmaErr[i] = -999;
            if (ff)
            {
                TRMean[i] = ff->GetParameter(1) * 1e3;
                TRMeanErr[i] = ff->GetParError(1) * 1e3;
                TRSigma[i] = ff->GetParameter(2) * 1e3;
                TRSigmaErr[i] = ff->GetParError(2) * 1e3;
            }
            ht[i]->GetXaxis()->SetNdivisions(505);
            gPad->SetLeftMargin(0.2);
            gPad->SetRightMargin(0.01);
            gPad->SetTopMargin(0.01);
            gPad->SetBottomMargin(0.2);

            sprintf(buff, "CH%d", i + 1);
            DrawMyLatex(buff, 0.5, 0.9, 42, 0.06, 2);

            sprintf(buff, "Sigma=%.0f#pm%.0fps", TRSigma[i], TRSigmaErr[i]);
            la = DrawMyLatex(buff, 0.55, 0.25);
            sprintf(buff, "Mean=%.0f#pm%.0fps", TRMean[i], TRMeanErr[i]);
            la = DrawMyLatex(buff, 0.55, 0.35);
        }
        sprintf(buff, "%s_Time.png", name);
        c1->SaveAs(buff);

        /* 
        //
        // -- print 1st pe arrival time distrubutions --
        //
        c1 = cdC(CCounter++);
        c1->Divide(2, 2, 0.01, 0.01);
        for (int i = 0; i < 4; i++)
        {
        c1->cd(i + 1);
        h1st[i]->Draw();
        DrawMyHist(h1st[i], "", "");
        h1st[i]->GetXaxis()->SetNdivisions(505);
        gPad->SetLeftMargin(0.2);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.2);
        sprintf(buff, "Mean=%.0fps", h1st[i]->GetMean() * 1e3);
        DrawMyLatex(buff, 0.6, 0.5, 42, 0.06, 2);
        sprintf(buff, "RMS=%.0fps", h1st[i]->GetRMS() * 1e3);
        DrawMyLatex(buff, 0.6, 0.4, 42, 0.06, 2);
        }
        sprintf(buff, "%s_1stPETime.png", name);
        c1->SaveAs(buff);

        //
        // -- print AMP distributions --
        //
        c1 = cdC(CCounter++);
        c1->Divide(2, 2, 0.01, 0.01);
        for (int i = 0; i < 4; i++)
        {
        c1->cd(i + 1);
        hU[i]->Draw();
        DrawMyHist(hU[i], "", "");
        hU[i]->GetXaxis()->SetNdivisions(505);
        gPad->SetLeftMargin(0.2);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.2);
        sprintf(buff, "Mean=%.0fps", hU[i]->GetMean());
        DrawMyLatex(buff, 0.6, 0.5, 42, 0.06, 2);
        }
        sprintf(buff, "%s_Amp.png", name);
        c1->SaveAs(buff);
        */
    }

    void Output(const char *rootname = "x2y2z1_model2blackwrap")
    {
        gStyle->SetOptFit(111);
        char buff[1024];
        sprintf(buff, "%s/%sdata.root", path, rootname);
        cout << "Check File :" << buff << endl;
        if (gSystem->AccessPathName(buff))
        {
            cout << "==> the file isn't exist!" << endl;
            CalculateTR(rootname, 0.2, "CFD", 1);
        }
        sprintf(buff, "%sdata", rootname);
        PrintResults(buff);
    }

    void readdatresults()
    {
        gStyle->SetOptTitle(1);
        int theta[] = {5, 10, 15, 20, 25, 30, 35};
        int phi[] = {45, 60, 75, 90, 105, 120, 135};
        const int T = 4;
        double NPE[T] = {0};
        double NPEErr[T] = {0};
        double TRSigma[T] = {0};
        double TRSigmaErr[T] = {0};
        double InX[T] = {0}; //incident position
        double InY[T] = {0};
        double InZ[T] = {0};
        double InPX[T] = {0}; //incident direction
        double InPY[T] = {0};
        double InPZ[T] = {0};

        char buff[200];
        char name[200];
        TH2D *hTR[T];
        TH2D *hNPE[T];
        ifstream input;
        int xNum = sizeof(phi) / sizeof(phi[0]);
        double xMin = phi[0] - 0.5 * (phi[1] - phi[0]);
        double xMax = phi[xNum - 1] + 0.5 * (phi[1] - phi[0]);

        int yNum = sizeof(theta) / sizeof(theta[0]);
        double yMin = theta[0] - 0.5 * (theta[1] - theta[0]);
        double yMax = theta[yNum - 1] + 0.5 * (theta[1] - theta[0]);

        for (int i = 0; i < T; i++)
        {
            sprintf(buff, "hTR%d", i);
            sprintf(name, "TR(ps),CH=%d", i);
            hTR[i] = new TH2D(buff, name, xNum, xMin, xMax, yNum, yMin, yMax);
            sprintf(buff, "hNPE%d", i);
            sprintf(name, "NPE,CH=%d", i);
            hNPE[i] = new TH2D(buff, name, xNum, xMin, xMax, yNum, yMin, yMax);
            ifstream input;
            for (int j = 0; j < xNum; j++)
            {
                for (int k = 0; k < yNum; k++)
                {

                    sprintf(buff, "%s/model1_theta%dphi%ddata.dat", path, theta[j], phi[k]);
                    cout << " ==> Read data file: " << buff << endl;
                    if (gSystem->AccessPathName(buff))
                    {
                        cout << "the file isn't exist!" << endl;
                        return;
                    }
                    input.open(buff);
                    input >> InX[0] >> InY[0] >> InZ[0] >> InPX[0] >> InPY[0] >> InPZ[0];
                    for (int ich = 0; ich < T; ich++)
                    {

                        input >> NPE[ich] >> NPEErr[ich] >> TRSigma[ich] >> TRSigmaErr[ich];
                    }

                    cout << NPE[i] << "\t" << TRSigma[i] << endl;
                    //hTR[i]->Fill(phi[k],theta[j]);
                    //cout<<reTrkLTRSigma<<endl;

                    TRSigma[i] = URound(TRSigma[i], 1);
                    hTR[i]->SetBinContent(k + 1, j + 1, TRSigma[i]);

                    cout << " TR " << i << "\t:" << TRSigma[i] << endl;

                    NPE[i] = URound(NPE[i], 0);
                    hNPE[i]->SetBinContent(k + 1, j + 1, NPE[i]);
                    input.close();
                }
            }
        }

        //return;
        TCanvas *c;
        int cNum = 1;

        //char prefix[100] = "nogroundpi";
        for (int i = 0; i < T; i++)
        {

            c = cdC(cNum++);
            DrawMy2dHist(hTR[i], "#phi (degree)", "#theta (degree)");
            hTR[i]->Draw("colz");
            //TExec *ex1 = new TExec("ex1", "Pal();");
            TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(55, 0, 0.45);");
            ex1->Draw();

            hTR[i]->Draw("colz TEXT SAME");
            hTR[i]->SetMarkerSize(2.8);
            sprintf(buff, "%s/CH%dTR.png", path, i);
            c->Update();
            c->Modified();
            c->SaveAs(buff);

            c = cdC(cNum++);
            DrawMy2dHist(hNPE[i], "#phi (degree)", "#theta (degree)");
            hNPE[i]->Draw("colz");
            TExec *ex2 = new TExec("ex2", "gStyle->SetPalette(55, 0, 0.45);");
            ex2->Draw();
            hNPE[i]->Draw("colz text SAME");
            hNPE[i]->SetMarkerSize(2.8);
            sprintf(buff, "%s/CH%dNPE.png", path, i);
            c->SaveAs(buff);
        }
    }

    /*
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

    void DrawMyHist1(TH1 *datahist, char *xtitle, char *ytitle, Color_t LColor = 1, Width_t LWidth = 3, Color_t TitleColor = 1)
    {
        datahist->SetLineColor(LColor);
        datahist->SetLineWidth(LWidth);
        datahist->GetXaxis()->SetTitle(xtitle);
        datahist->GetYaxis()->SetTitle(ytitle);
        datahist->GetXaxis()->SetAxisColor(1);
        datahist->GetYaxis()->SetAxisColor(1);
        datahist->GetXaxis()->SetLabelColor(1);
        datahist->GetYaxis()->SetLabelColor(1);
        datahist->GetXaxis()->SetLabelFont(42);
        datahist->GetYaxis()->SetLabelFont(42);
        datahist->GetXaxis()->SetLabelSize(0.06);
        datahist->GetYaxis()->SetLabelSize(0.06);
        datahist->GetXaxis()->SetLabelOffset(0.01);
        datahist->GetYaxis()->SetLabelOffset(0.01);
        datahist->GetXaxis()->SetTitleFont(42);
        datahist->GetYaxis()->SetTitleFont(42);
        datahist->GetXaxis()->SetTitleColor(TitleColor);
        datahist->GetYaxis()->SetTitleColor(TitleColor);
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
