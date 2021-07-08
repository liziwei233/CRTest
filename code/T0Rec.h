#ifndef T0Rec_h
#define T0Rec_h

#include "CRdata.h"
#include "LUT.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMatrixD.h"

//#include "T0Rec.h"

using namespace TMath;
using namespace std;
class T0Rec : public CRdata
{
public:
  T0Rec(string in = "data.root",
        string out = "T0Rec.root",
        TTree *tree = 0);
  ~T0Rec();
  void InitialOutput();
  void InitialRefcounter();
  enum lightpath
  {
    side0,
    side1,
    side2,
    side3
  };
  lightpath Mirror(lightpath);

  /// ReConstruction
  void MirrorTransformation(lightpath);
  void RecHitPos();
  void RecTOF(double);
  void Reconstruction(double,double);

  /// Set Hit
  void Clear();
  void SetTrackHit();
  void SetSimuPhotonHit(int i);
  void SetPhotonHit(int i);
  void Simu2Lab(double &x, double &y, double &z);
  void GetPhotonDirection();
  double GetRemainder(double, double);
  void Loop();
  static bool P(double x, double y){return x>y;};
  void myswap(double &a, double &b);

  /// LUT
  LUT LUT_n;

  /// Hit
  double gammaDx[4];
  double gammaDy[4];
  double gammaDz[4];
  double Dx;
  double Dy;
  double Dz;
  double Px;
  double Py;
  double Pz;
  double FL;
  int fID;
  int ChX;
  int ChY;
  double T;
  double TTS, T0;

  /// result of Reconstruction
  double HitX, HitY;
  std::vector<double> RecDeltaX;
  std::vector<double> RecDeltaY;
  std::vector<double> RecDircX;
  std::vector<double> RecDircY;
  std::vector<double> RecDircZ;
  std::vector<double> RecPropLength;
  std::vector<double> RecT0detTime;
  vector<double> Weight;
  std::vector<std::vector<double>> RecT0detTime_Vec;
  std::vector<std::vector<double>> Weight_Vec;
  double HyT0detTime;
  double BestPropLength;
  double BestT0detTime;
  vector<double> RBT0vec;
  vector<double> LOPvec;
  vector<double> TOPvec;
  vector<double> trueTOPvec;
  vector<int> chidvec;
  vector<int> Xvec;
  vector<int> Yvec;
  vector<double> TT;
  vector<double> AA;
  double RBT0[4];
  double LOP[4];
  double TOP[4];
  double trueTOP[4];

  vector<double> *pTOFvec = &RBT0vec;
  vector<double> *pLOPvec = &LOPvec;
  vector<double> *pTOPvec = &TOPvec;
  vector<double> *ptrueTOPvec = &trueTOPvec;
  vector<int> *pchidvec = &chidvec;
  vector<int> *pXvec = &Xvec;
  vector<int> *pYvec = &Yvec;

  double SimuT0dettime;
  double MeanT0time;
  double MeanTL;
  int side0ctr;
  int side1ctr;
  int side2ctr;
  int side3ctr;
  /// parameters
  double T0mediumL;
  double T0mediumt;
  double SensorR;

  double A, B, C;
  double NpMax, Np, Ng;
  double Mass_Mu;
  int MaxReflection;
  double PrimaryMomentum;

  /// random number
  TRandom3 *fRandom;
  TTree *ot;
  ofstream op;
  
  // histgram
  TH1D* hRBT0;
};

T0Rec::T0Rec(string in, string out, TTree *tree) : CRdata(in, out, tree), LUT_n("LUT/n.txt", 400, 2)
{
  InitialOutput();
  Mass_Mu = 105.65837; //MeV
  PrimaryMomentum = 2; //GeV
  NpMax = 1.48779, Np = 1.46979, Ng = 1.51404;
  MaxReflection = 2;

  T0mediumL = 180;
  T0mediumt = 10;
  SensorR = 5;

  cout << "[+] - Rebuild parameters have been set!" << endl;
}
T0Rec::~T0Rec()
{
}
void T0Rec::InitialOutput()
{
  hRBT0 = new TH1D("hRBT0","",200,0,3);
  op.open("Rectoflog.dat");
  ot = new TTree("data", "restore T0 rebuilded data  from G4");
  ot->Branch("MeanT0time", &MeanT0time, "MeanT0time/D");
  ot->Branch("HyT0detTime", &HyT0detTime, "HyT0detTime/D");
  ot->Branch("SimuT0dettime", &SimuT0dettime, "SimuT0dettime/D");
  ot->Branch("Px", &Px, "Px/D");
  ot->Branch("Py", &Py, "Py/D");
  ot->Branch("Pz", &Pz, "Pz/D");

  ot->Branch("RBT0",RBT0,"RBT0[4]/D");
  ot->Branch("LOP",LOP,"LOP[4]/D");
  ot->Branch("TOP",TOP,"TOP[4]/D");
  ot->Branch("trueTOP",trueTOP,"trueTOP[4]/D");

/*
  ot->Branch("RBT0", pTOFvec);
  ot->Branch("LOP", pLOPvec);
  ot->Branch("TOP", pTOPvec);
  ot->Branch("trueTOP", ptrueTOPvec);
*/
  ot->Branch("chid", pchidvec);
  ot->Branch("Xvec", pXvec);
  ot->Branch("Yvec", pYvec);

  ot->Branch("CRRBtheta", &CRRBtheta);
  ot->Branch("CRE", &CRE);
}
T0Rec::lightpath T0Rec::Mirror(T0Rec::lightpath a)
{
  int mirrorside = (a + 2) % 4;
  return (T0Rec::lightpath)mirrorside;
}

void T0Rec::MirrorTransformation(T0Rec::lightpath a)
{

  if (a == fID)
    return;
  if (a == side0)
  {
    side0ctr++;
    double x1 = T0mediumL / 2;
    double y1 = T0mediumL / 2;
    double x2 = T0mediumL / 2;
    double y2 = -T0mediumL / 2;

    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    double tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    double tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;

    Dx = -Dx;
    Dy = Dy;
  }
  if (a == side1)
  {
    side1ctr++;
    double x1 = T0mediumL / 2;
    double y1 = -T0mediumL / 2;
    double x2 = -T0mediumL / 2;
    double y2 = -T0mediumL / 2;

    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    double tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    double tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;

    Dx = Dx;
    Dy = -Dy;
  }
  if (a == side2)
  {
    side2ctr++;
    double x1 = -T0mediumL / 2;
    double y1 = T0mediumL / 2;
    double x2 = -T0mediumL / 2;
    double y2 = -T0mediumL / 2;

    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    double tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    double tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;

    Dx = -Dx;
    Dy = Dy;
  }

  if (a == side3)
  {
    side3ctr++;
    double x1 = T0mediumL / 2;
    double y1 = T0mediumL / 2;
    double x2 = -T0mediumL / 2;
    double y2 = T0mediumL / 2;

    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    double tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    double tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;

    Dx = Dx;
    Dy = -Dy;
  }
}

void T0Rec::Clear()
{
  RecDeltaX.clear();
  RecDeltaY.clear();
  RecDircX.clear();
  RecDircY.clear();
  RecDircZ.clear();
  RecPropLength.clear();
  RecT0detTime.clear();
  Weight.clear();
}
void T0Rec::Simu2Lab(double &x, double &y, double &z)
{
  double Sdata[] = {x, y, z};
  double Mdata[] = {0, 0, -1, 1, 0, 0, 0, -1, 0};
  //double Mdata[] = {0, 0, 1, 0, 1, 0, -1,0, 0};
  TMatrixD M(3, 3, Mdata);
  TMatrixD S(1, 3, Sdata);
  TMatrixD L = S * M;
  x = TMatrixDRow(L, 0)(0);
  y = TMatrixDRow(L, 0)(1);
  z = TMatrixDRow(L, 0)(2);
}

void T0Rec::myswap(double &a, double &b)
{
    double temp = 0;
    temp = a;
    a = b;
    b = temp;
}
double T0Rec::GetRemainder(double D, double d)
{
  int n = D / d;
  return D - n * d;
}
void T0Rec::InitialRefcounter()
{
  side0ctr = 0;
  side1ctr = 0;
  side2ctr = 0;
  side3ctr = 0;
}
void T0Rec::SetTrackHit()
{
  //T0 = fRandom->Gaus(0,0.04);
  // TODO: simulation to lab
  /*
  Dx = CRRBpx;
  Dy = CRRBpy;
  Dz = CRRBpz;
  Simu2Lab(Dx,Dy,Dz);
  Px = FTOFdetRBx;
  Py = FTOFdetRBy;
  Pz = FTOFdetRBz;
  Simu2Lab(Px,Py,Pz);
  */
  Dx = CRpx;
  Dy = CRpy;
  Dz = CRpz;
  Px = T0detx - T0mediumt / 2;
  Py = T0dety - T0mediumt / 2 / Dx * Dy;
  Pz = T0detz - T0mediumt / 2 / Dx * Dz;
  //double T0Px = T0detx -T0mediumt/2;
  //double T0Py = T0dety -T0mediumt/2/Dx*Dy;
  //double T0Pz = T0detz -T0mediumt/2/Dx*Dz;
  //FL = sqrt(Power(Px - T0Px, 2) + Power(Py - T0Py, 2) + Power(Pz - T0Pz, 2));

  Simu2Lab(Dx, Dy, Dz);
  Simu2Lab(Px, Py, Pz);
  double DieL = Sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
  Dx = Dx / DieL;
  Dy = Dy / DieL;
  Dz = Dz / DieL;

  //cout << "FL="<<FL<<endl;
  InitialRefcounter();
}
void T0Rec::GetPhotonDirection()
{
  vector<double> gammadx[4];
  vector<double> gammady[4];
  vector<double> gammadz[4];
  vector<double> gammat[4];
  for(int i =0 ; i<T0photont.size();i++)
  {
     gammadx[T0photonid[i]].push_back(T0photonpx[i]);
     gammady[T0photonid[i]].push_back(T0photonpy[i]);
     gammadz[T0photonid[i]].push_back(T0photonpz[i]);
     gammat[T0photonid[i]].push_back(T0photont[i]);
     
  }
  for(int i =0; i<4;i++){
    if(gammat[i].size()<1) continue;
    for (int kk = 0; kk < gammat[i].size()-1; kk++)
        {
            for (int jj = 0; jj < gammat[i].size()-1 - kk; jj++)
            {
                if (gammat[i][jj] > gammat[i][jj+1])
                {
                    
                    myswap(gammat[i][jj], gammat[i][jj+1]);
                    myswap(gammadx[i][jj], gammadx[i][jj+1]);
                    myswap(gammady[i][jj], gammady[i][jj+1]);
                    myswap(gammadz[i][jj], gammadz[i][jj+1]);
                    
                }
            }
        }
      gammaDx[i] = gammadx[i][0];
      gammaDy[i] = gammady[i][0];
      gammaDz[i] = gammadz[i][0];
      Simu2Lab( gammaDx[i], gammaDy[i], gammaDz[i]);
      gammaDx[i] = gammaDx[i]/gammaDy[i];
      gammaDy[i] = gammaDy[i]/gammaDy[i];
      gammaDz[i] = gammaDz[i]/gammaDy[i];
  }
}
void T0Rec::SetSimuPhotonHit(int i)
{
  //fID = FTOFeleid[i];
  //T = FTOFelefittime1[i] - T0dett; // T= TOF(T0-FTOF) + TOP
  if(T0eleid[i]==1) fID=0;
  else if(T0eleid[i]==2) fID=3;
  else if(T0eleid[i]==3) fID=2;
  else if(T0eleid[i]==0) fID=1;
  //fID = T0photonid[i];
  T = T0fasttime[i]; // T= TOF(T0-FTOF) + TOP
  /* for lab */
  ChX = TMath::Power(-1, (fID / 2 + 1)) * (fID % 2 - 1);
  ChY = TMath::Power(-1, (fID / 2 + 1)) * (fID % 2);

  /* for simu */
  //ChX = TMath::Power(-1, (fID / 2 )) * (fID % 2 );
  //ChY = TMath::Power(-1, (fID / 2 )) * (fID % 2-1);
  op << "T= " << T << endl;
  op << "FID= " << fID << ",(" << ChX << "," << ChY << ")" << endl;
  op << "photon direction: ("<<gammaDx[T0eleid[i]]<<","<<gammaDy[T0eleid[i]]<<","<<gammaDz[T0eleid[i]]<<")"<<endl;
  //T = GlobalTime->at(i)+TTS-T0;
  //cout << fID << ",ChX= " << ChX << ", ChY= " << ChY << endl;
}
void T0Rec::SetPhotonHit(int i)
{
  //fID = FTOFeleid[i];
  //T = FTOFelefittime1[i] - T0dett; // T= TOF(T0-FTOF) + TOP
 if(T0eleid[i]==1) fID=3;
  else if(T0eleid[i]==3) fID=1;
  else fID = T0eleid[i];
  //fID = T0photonid[i];
  T = T0elefittime1[i]+0.3; // T= TOF(T0-FTOF) + TOP

  ChX = Power(-1, (fID / 2 + 1)) * (fID % 2 - 1);
  ChY = Power(-1, (fID / 2 + 1)) * (fID % 2);
  op << "T= " << T << endl;
  op << "FID= " << fID << ",(" << ChX << "," << ChY << ")" << endl;
  //T = GlobalTime->at(i)+TTS-T0;
  //cout << fID << ",ChX= " << ChX << ", ChY= " << ChY << endl;
}
void T0Rec::RecHitPos()
{
  HitX = ChX * 90;
  HitY = ChY * 90;
  op << "Hit pos (" << HitX << "," << HitY << ")" << endl;
}

void T0Rec::RecTOF(double beta)
{

  int i = RecDeltaX.size();
  double CosThetaC = 1. / Np / beta;
  RecDeltaX.push_back(HitX - Px);
  RecDeltaY.push_back(HitY - Py);
  RecDircX.push_back(RecDeltaX[i] / Abs(RecDeltaY[i]));
  RecDircY.push_back(RecDeltaY[i] / Abs(RecDeltaY[i]));
  RecDircZ.push_back(0);
  RecPropLength.push_back(0);
  RecT0detTime.push_back(0);

  //cout << "Hit (x,y)=(" << HitX << "," << HitY << ")" << endl;
  //cout << " (Px, Py)=(" << Px << "," << Py << ")" << endl;
  //cout << "RecDelta =(" << RecDeltaX[i] << "," << RecDeltaY[i] << ")" << endl;
  A = Dz * Dz - CosThetaC * CosThetaC;
  B = 2. * Dz * (Dx * RecDircX[i] + Dy * RecDircY[i]);
  C = Power(Dx * RecDircX[i] + Dy * RecDircY[i], 2) - CosThetaC * CosThetaC * (Power(RecDircX[i], 2) + Power(RecDircY[i], 2));

  double Delta = B * B - 4. * A * C;
  if (Delta >= 0)
  {
    double Z1 = (-B + Sqrt(Delta)) / 2. / A;
    double Z2 = (-B - Sqrt(Delta)) / 2. / A;
    double L_Z1 = Abs(RecDeltaY[i] / RecDircY[i] * Z1);
    double L_Z2 = Abs(RecDeltaY[i] / RecDircY[i] * Z2);
    int tag1 = 0, tag2 = 0;
    double costheta1 = Dx * RecDircX[i] + Dy * RecDircY[i] + Dz * Z1;
    double totalref1 = (Power(RecDircX[i], 2) + Power(RecDircY[i], 2)) / (Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Z1 * Z1) * NpMax * NpMax / 0.95;
    double costheta2 = Dx * RecDircX[i] + Dy * RecDircY[i] + Dz * Z2;
    double totalref2 = (Power(RecDircX[i], 2) + Power(RecDircY[i], 2)) / (Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Z2 * Z2) * NpMax * NpMax / 0.95;
    double R1 = GetRemainder(L_Z1, 2 * T0mediumt);
    double R2 = GetRemainder(L_Z2, 2 * T0mediumt);
    //op << "costheta :" << costheta1 << "\t" << costheta2 << endl;
    //op << "totalref :" << totalref1 << "\t" << totalref2 << endl;
    //op << "Remainder :" << R1 << "\t" << R2 << endl;
    if (costheta1 > 0 && totalref1 >= 1)
      tag1 = 1;
    if (costheta2 > 0 && totalref2 >= 1)
      tag2 = 1;
    if (tag1 || tag2)
    {
      if (tag1 && tag2)
      {

        RecDircZ[i] = Abs(Z1) < Abs(Z2) ? Z1 : Z2;
      }
      if (tag1 && !tag2)
      {

        RecDircZ[i] = Z1;
      }
      if (tag2 && !tag1)
      {

        RecDircZ[i] = Z2;
      }
      RecPropLength[i] = Abs(RecDeltaY[i] * Sqrt(Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Power(RecDircZ[i], 2)) / RecDircY[i]);
      RecT0detTime[i] = T - RecPropLength[i] / TMath::C() * Ng * 1e6;
    }
    op << side0ctr << "\t" << side1ctr << "\t" << side2ctr << "\t" << side3ctr << "\t" << RecDircX[i] << "\t" << RecDircY[i] << "\t" << RecDircZ[i] << "\t" << RecDeltaX[i] << "\t" << RecDeltaY[i] << "\t" << Abs(RecDeltaY[i] / RecDircY[i] * RecDircZ[i]) << "\t" << RecPropLength[i] << "\t" << RecT0detTime[i] << "\t" << SimuT0dettime << endl;
    hRBT0->Fill(RecT0detTime[i]);
    //op << side0ctr << "\t" << side1ctr << "\t" << side2ctr << "\t" << side3ctr << "\t" << tag1 << "\t" << RecDircX[i] << "\t" << RecDircY[i] << "\t" << Z1 << "\t" << RecDeltaX[i] << "\t" << RecDeltaY[i] << "\t" << L_Z1 << "\t" << RecPropLength[i] << "\t" << RecT0detTime[i] << "\t" << SimuT0dettime << endl;
    //op << side0ctr << "\t" << side1ctr << "\t" << side2ctr << "\t" << side3ctr << "\t" << tag2 << "\t" << RecDircX[i] << "\t" << RecDircY[i] << "\t" << Z2 << "\t" << RecDeltaX[i] << "\t" << RecDeltaY[i] << "\t" << L_Z2 << "\t" << RecPropLength[i] << "\t" << RecT0detTime[i] << "\t" << SimuT0dettime << endl;
  }
}
void T0Rec::Reconstruction(double mass,double Ek)
{
  Clear();
  //double beta = 1. / Sqrt(1 + Power(mass / Ek, 2));
  double beta = Ek/Sqrt(Ek*Ek+mass*mass);
  //HyFlightTime = FL / TMath::C() / beta * 1e6;
  //cout << "[+] | - Start Reconstruction " << endl;
  //cout<< "beta= "<<beta<<endl;
  RecHitPos();
  //RecTOF(beta);

  // *direct
  RecTOF(beta);
  Weight.push_back(0.69);
  // *mirror 1 time
  for (int j = 0; j < 4; j++)
  {

    if (j == fID)
      continue;
    MirrorTransformation((enum lightpath)j);
    RecTOF(beta);
    Weight.push_back(0.199);
    SetTrackHit();
  }

  /*
  // *mirror 2 time
  for (int j = 0; j < 4; j++)
  {
    if (j == fID || j == (fID + 2) % 4)
      continue;
    MirrorTransformation((enum lightpath)j);
    MirrorTransformation(Mirror((enum lightpath)j));
    RecTOF(beta);
    Weight.push_back(0.068);
    SetTrackHit();
  }
*/
  // *mirror 2 time
  for (int j = 0; j < 4; j++)
  {
    if (j == fID)
      continue;
    for (int s = 0; s < 4; s++)
    {
      if (s == j || s == fID)
        continue;
      MirrorTransformation((enum lightpath)j);
      MirrorTransformation((enum lightpath)s);
      RecTOF(beta);
      Weight.push_back(0.068);
      SetTrackHit();
    }
  }
  // *mirror 3 time
  for (int j = 0; j < 4; j++)
  {
    if (j == fID || j == (fID + 2) % 4)
      continue;
    MirrorTransformation((enum lightpath)j);
    MirrorTransformation(Mirror((enum lightpath)j));
    MirrorTransformation((enum lightpath)j);
    RecTOF(beta);
    Weight.push_back(0.03);
    SetTrackHit();
  }

  MirrorTransformation((enum lightpath)((fID + 1) % 4));
  MirrorTransformation((enum lightpath)((fID + 2) % 4));
  MirrorTransformation((enum lightpath)((fID + 3) % 4));
  RecTOF(beta);
  Weight.push_back(0.03);
  SetTrackHit();
  MirrorTransformation((enum lightpath)((fID + 3) % 4));
  MirrorTransformation((enum lightpath)((fID + 2) % 4));
  MirrorTransformation((enum lightpath)((fID + 1) % 4));
  RecTOF(beta);
  Weight.push_back(0.03);
  SetTrackHit();

  // *mirror 4 time
  for (int j = 0; j < 4; j++)
  {
    if (j == fID || j == (fID + 2) % 4)
      continue;
    MirrorTransformation((enum lightpath)j);
    MirrorTransformation(Mirror((enum lightpath)j));
    MirrorTransformation((enum lightpath)j);
    MirrorTransformation(Mirror((enum lightpath)j));
    RecTOF(beta);
    Weight.push_back(0.03);
    SetTrackHit();
  }

  /*
  int tag = 0;
  for (int id = 0; id < (int)RecT0detTime.size(); id++)
  {
    //if(Abs(RecT0detTime[id]-HyFlightTime) < Abs(RecT0detTime[tag]-HyFlightTime)) tag = id;
    if (Abs(RecT0detTime[id] - HyFlightTime) < Abs(RecT0detTime[tag] - HyFlightTime))
      tag = id;
  }
  BestPropLength = RecPropLength[tag];
  BestT0detTime = RecT0detTime[tag];
*/

  //cout << "BestPropLength= " << BestPropLength<< endl
  //     << "BestT0detTime= " << BestT0detTime << endl;
}

void T0Rec::Loop()
{
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();

  //double FT[16*sensorN];
  //double TL[16*sensorN];
  TT.clear();
  AA.clear();
  cout << "[+] | - nentries = " << nentries << endl;
  int exctr = 0;
  for (Long64_t tentry = 0; tentry < nentries&&exctr<1e3; tentry++)
  {
    fChain->GetEntry(tentry);
    if (T0detRBy == -999 ||CRRBtheta < 0)
      continue;

    RBT0vec.clear();
    LOPvec.clear();
    TOPvec.clear();
    chidvec.clear();
    trueTOPvec.clear();
    Xvec.clear();
    Yvec.clear();
    SimuT0dettime = T0dett;
    SetTrackHit();
    op << "Track diraction (Dx,Dy,Dz) = (" << Dx << "," << Dy << "," << Dz << ")" << endl;
    op << "Track hit position (Px,Py) = (" << Px << "," << Py << ")" << endl;
    //int N = T0photonid.size();
   op << "Entry$= " << tentry << endl;
   GetPhotonDirection();
    double sum = 0;
    double wsum = 0;
    HyT0detTime = 0;
    RecT0detTime_Vec.clear();
    Weight_Vec.clear();
    int N=4;
    for (int i = 0; i < N; i++)
    {
      if (T0eletot[i] < 0)
        continue;

      SetSimuPhotonHit(i);
      //SetPhotonHit(i);
      Reconstruction(Mass_Mu / 1e3,CRE);
      sort(RecT0detTime.begin(), RecT0detTime.end(),P);
      RecT0detTime_Vec.push_back(RecT0detTime);
      Weight_Vec.push_back(Weight);
      chidvec.push_back(fID);
      for (int j = 0; j < RecT0detTime_Vec.back().size(); j++)
      {
        if (RecT0detTime_Vec.back()[j] < 0)
          continue;
        sum += RecT0detTime_Vec.back()[j] * Weight_Vec.back()[j];
        wsum += Weight_Vec.back()[j];
      }
    };
    if (RecT0detTime_Vec.size() < 3)
      continue;
    HyT0detTime = sum / wsum;

    double HyT0detTime_vec[N];
HyT0detTime = SimuT0dettime;
    int ctr = 0;
    /*
    for (int i = 0; i < RecT0detTime_Vec.size(); i++)
    {
      HyT0detTime_vec[i] = RecT0detTime_Vec[i][0];
    }
    HyT0detTime = Mean(RecT0detTime_Vec.size(), HyT0detTime_vec);

    for (int s = 0; s < 5; s++)
    {

      for (int i = 0; i <  RecT0detTime_Vec.size(); i++)
      {
        int tag = 0;
        bool flag = 0;
        for (int j = 0; j < RecT0detTime_Vec[i].size(); j++)
        {
          if (Abs(RecT0detTime_Vec[i][j] - HyT0detTime) < Abs(RecT0detTime_Vec[i][tag] - HyT0detTime))
          {

            tag = j;
            flag = 1;
          }
        }
        if (flag)
        {

          HyT0detTime_vec[i] = RecT0detTime_Vec[i][tag];
          HyT0detTime = Mean(RecT0detTime_Vec.size(), HyT0detTime_vec);
        }
      }
    }
    */
    for (int i = 0; i < RecT0detTime_Vec.size(); i++)
    {
      int tag = 0;
      for (int j = 0; j < RecT0detTime_Vec[i].size(); j++)
      {
        if (Abs(RecT0detTime_Vec[i][j] - HyT0detTime) < Abs(RecT0detTime_Vec[i][tag] - HyT0detTime))
        {

          tag = j;
        }
      }
      BestPropLength = RecPropLength[tag];
      BestT0detTime = RecT0detTime_Vec[i][tag];
      op << T0photonid[i] << "\t" << BestPropLength << "\t" << BestT0detTime << "\t" << SimuT0dettime << "\t" << HyT0detTime << endl;
      RBT0[chidvec[i]] = BestT0detTime;
      LOP[chidvec[i]] = BestPropLength;
      TOP[chidvec[i]] = BestPropLength / TMath::C() * Ng * 1e6;
      trueTOP[chidvec[i]] = FTOFphotonTOP[i];
      RBT0vec.push_back(BestT0detTime);
      LOPvec.push_back(BestPropLength);
      TOPvec.push_back(BestPropLength / TMath::C() * Ng * 1e6);
      trueTOPvec.push_back(FTOFphotonTOP[i]);
      
      Xvec.push_back(ChX);
      Yvec.push_back(ChY);
    }

    int cut = 2;
    if (RBT0vec.size() > cut)
    {
      sort(RBT0vec.begin(), RBT0vec.end());
      double sum = 0, count = 0;
      double sum2 = 0;
      for (int i = cut / 2; i < RBT0vec.size() - cut / 2; i++)
      {
        if (TOPvec.at(i) <= 0.1)
          continue;
        if (abs(RBT0vec.at(i) - HyT0detTime) > 0.3)
          continue;
        sum += RBT0vec.at(i);
        count += 1;
        sum2 += LOPvec.at(i);
      }
      if (count != 0)
      {
        MeanT0time = sum / count;
        MeanTL = sum2 / count;
        TT.push_back(MeanT0time);
        AA.push_back(MeanTL);
        ot->Fill();
        exctr++;
      }
    }
  }
  hRBT0->Draw();
  fout->WriteTObject(ot);
}

#endif