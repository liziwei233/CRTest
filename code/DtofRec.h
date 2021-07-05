#ifndef DtofRec_h
#define DtofRec_h

#include "CRdata.h"
#include "LUT.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMatrixD.h"
//#include "DtofRec.h"

using namespace TMath;
using namespace std;
class DtofRec : public CRdata
{
public:
  DtofRec(string in = "data.root",
          string out = "Rec.root",
          TTree *tree = 0);
  ~DtofRec();
  void InitialOutput();
  void InitialRefcounter();
  enum lightpath
  {
    direct,
    left,
    right,
    front
  };
  lightpath Mirror(lightpath);

  /// ReConstruction
  void MirrorTransformation(lightpath, bool);
  void RecHitPos();
  void RecTOF(double);
  void Reconstruction(double);

  /// Set Hit
  void Clear();
  void SetTrackHit();
  void SetPhotonHit(int i);
  void Simu2Lab(double &x, double &y, double &z);
  double GetRemainder(double, double);
  void Loop();
  /// LUT
  LUT LUT_n;

  /// Hit
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
  std::vector<double> RecFlightTime;
  double HyFlightTime;
  double BestPropLength;
  double BestFlightTime;
  vector<double> TOFvec;
  vector<double> LOPvec;
  vector<double> TOPvec;
  vector<double> trueTOPvec;
  vector<int> chidvec;
  vector<int> Xvec;
  vector<int> Yvec;
  vector<double> TT;
  vector<double> AA;

  vector<double>* pTOFvec = &TOFvec;
  vector<double>* pLOPvec = &LOPvec;
  vector<double>* pTOPvec = &TOPvec;
  vector<double>* ptrueTOPvec = &trueTOPvec;
  vector<int>* pchidvec = &chidvec;
  vector<int>* pXvec = &Xvec;
  vector<int>* pYvec = &Yvec;

  
  double SimuFT;
  double RefFT;
  double MeanFT;
  double MeanTL;
  int Dctr;
  int Lctr;
  int Rctr;
  int Fctr;
  /// parameters
  double QuartzY1, QuartzY2, QuartzX, QuartzZ;
  double MCPPMToffset;
  double R10754_sensor_L;
  double R10754_sensor_edge_L;
  double R10754_anode_L;
  double R10754_anode_interval_L;
  double R10754_interval_L;
  int sensorN;
  double A, B, C;
  double NpMax, Np, Ng;
  double Mass_Mu;
  double interval, PhotonDetW, CathodeW; // mm
  double R, Gap, QuartzR1, QuartzR2;     // mm
  int SectorNu, PhotonDetNu;
  int MaxReflection;
  double PrimaryMomentum;

  /// random number
  TRandom3 *fRandom;
  TTree *ot;
  ofstream op;
};

DtofRec::DtofRec(string in, string out, TTree *tree) : CRdata(in, out, tree), LUT_n("LUT/n.txt", 400, 2)
{
  InitialOutput();
  Mass_Mu = 105.65837; //MeV
  PrimaryMomentum = 2; //GeV
  NpMax = 1.48779, Np = 1.46979, Ng = 1.51404;
  MaxReflection = 2;

  SectorNu = 12;
  sensorN = 6;
  QuartzY1 = 295;
  QuartzY2 = 533;
  QuartzX = 10;
  QuartzZ = 454;
  MCPPMToffset = 18;
  R10754_sensor_L = 27.6;
  R10754_sensor_edge_L = 2.5;
  R10754_anode_L = 5.28;
  R10754_anode_interval_L = 0.3;
  R10754_interval_L = 5;
  cout << "[+] - Rebuild parameters have been set!" << endl;
}

DtofRec::~DtofRec()
{
}

void DtofRec::InitialOutput()
{
  op.open("Rectoflog.dat");
  ot = new TTree("data", "restore Dtof rebuilded data  from G4");
  ot->Branch("MeanFT", &MeanFT, "MeanFT/D");
  ot->Branch("MeanTL", &MeanTL, "MeanTL/D");
  ot->Branch("SimuFT", &SimuFT, "SimuFT/D");
  ot->Branch("RefFT", &HyFlightTime, "RefFT/D");
  ot->Branch("Px", &Px, "Px/D");
  ot->Branch("Py", &Py, "Py/D");
  ot->Branch("Pz", &Pz, "Pz/D");

  ot->Branch("TOF",pTOFvec);
  ot->Branch("LOP",pLOPvec);
  ot->Branch("TOP",pTOPvec);
  ot->Branch("trueTOP",ptrueTOPvec);
  ot->Branch("chid",pchidvec);
  ot->Branch("Xvec",pXvec);
  ot->Branch("Yvec",pYvec);

  ot->Branch("CRRBtheta", &CRRBtheta);
  ot->Branch("CRE", &CRE);
}
DtofRec::lightpath DtofRec::Mirror(DtofRec::lightpath a)
{
  if (a == left)
    return right;
  if (a == right)
    return left;
  if (a == direct)
    return direct;
}

void DtofRec::MirrorTransformation(DtofRec::lightpath a, bool s)
{
  if (a == direct )
  {
    Dctr++;
    
  }
  if (a == left)
  {
    Lctr++;
    // dot and  line cross
    double angle, angleD, x1, y1, x2, y2, tmpX, tmpY;
    angle = Pi() / 2. - Pi() / SectorNu;
    //x1 = (QuartzR1 * Sin(angle) + interval / Cos(angle)) * Cos(Pi() / 4.) - QuartzR1 * Cos(angle) * Sin(Pi() / 4.);
    //y1 = (QuartzR1 * Sin(angle) + interval / Cos(angle)) * Sin(Pi() / 4.) + QuartzR1 * Cos(angle) * Cos(Pi() / 4.);
    //x2 = (QuartzR2 * Sin(angle) + interval / Cos(angle)) * Cos(Pi() / 4.) - QuartzR2 * Cos(angle) * Sin(Pi() / 4.);
    //y2 = (QuartzR2 * Sin(angle) + interval / Cos(angle)) * Sin(Pi() / 4.) + QuartzR2 * Cos(angle) * Cos(Pi() / 4.);
    x1 = QuartzY1 / 2;
    y1 = QuartzZ / 2;
    x2 = QuartzY2 / 2;
    y2 = -QuartzZ / 2;

    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;
    //cout << "after left mirror, Px,Py=(" << Px << "," << Py << ")" << endl;

    angle = Pi() / 4. + Pi() / SectorNu;
    //cout << "Dx,Dy=(" << Dx << "," << Dy << ")" << endl;
    if (Dy > 0)
      angleD = ATan(Dy / Dx);
    if (Dy < 0)
      angleD = ATan(Dy / Dx);
    tmpX = Sqrt(Dx * Dx + Dy * Dy) * Cos(2. * angle + angleD - Pi());
    tmpY = Sqrt(Dx * Dx + Dy * Dy) * Sin(2. * angle + angleD - Pi());
    Dx = tmpX;
    Dy = tmpY;
    //cout << "after left mirror, Dx,Dy=(" << Dx << "," << Dy << ")" << endl;
  }
  //* TODO
  if (a == right)
  {
    Rctr++;
    double angle, angleD, x1, y1, x2, y2, tmpX, tmpY;
    angle = Pi() / 2. - Pi() / SectorNu;
    //x1 = (QuartzR1 * Sin(angle) + interval / Cos(angle)) * Cos(Pi() / 4.) + QuartzR1 * Cos(angle) * Sin(Pi() / 4.);
    //y1 = (QuartzR1 * Sin(angle) + interval / Cos(angle)) * Sin(Pi() / 4.) - QuartzR1 * Cos(angle) * Cos(Pi() / 4.);
    //x2 = (QuartzR2 * Sin(angle) + interval / Cos(angle)) * Cos(Pi() / 4.) + QuartzR2 * Cos(angle) * Sin(Pi() / 4.);
    //y2 = (QuartzR2 * Sin(angle) + interval / Cos(angle)) * Sin(Pi() / 4.) - QuartzR2 * Cos(angle) * Cos(Pi() / 4.);
    x1 = -QuartzY1 / 2;
    y1 = QuartzZ / 2;
    x2 = -QuartzY2 / 2;
    y2 = -QuartzZ / 2;
    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;
    //cout << "after right mirror, Px,Py=(" << Px << "," << Py << ")" << endl;

    angle = Pi() / 4. - Pi() / SectorNu;
    //cout << "Dx,Dy=(" << Dx << "," << Dy << ")" << endl;
    if (Dy > 0)
      angleD = ATan(Dy / Dx);
    if (Dy < 0)
      angleD = ATan(Dy / Dx);
    tmpX = Sqrt(Dx * Dx + Dy * Dy) * Cos(2. * angle + angleD - Pi());
    tmpY = Sqrt(Dx * Dx + Dy * Dy) * Sin(2. * angle + angleD - Pi());
    Dx = tmpX;
    Dy = tmpY;
    //cout << "after right mirror, Dx,Dy=(" << Dx << "," << Dy << ")" << endl;
  }
  if (a == front)
  {
    Fctr++;
    double x1, y1, x2, y2, tmpX, tmpY;
    x1 = QuartzY1 / 2;
    y1 = -QuartzZ / 2;
    x2 = -QuartzY2 / 2;
    y2 = -QuartzZ / 2;

    A = y1 - y2;
    B = x2 - x1;
    C = x1 * y2 - x2 * y1;
    tmpX = ((B * B - A * A) * Px - 2. * A * B * Py - 2. * A * C) / (A * A + B * B);
    tmpY = ((A * A - B * B) * Py - 2. * A * B * Px - 2. * B * C) / (A * A + B * B);
    //cout << "Px,Py=(" << Px << "," << Py << ")" << endl;
    Px = tmpX;
    Py = tmpY;

    Dx = Dx;
    Dy = -Dy;
  }
  if (s)
  {
    MirrorTransformation(front, 0);
  }
}

void DtofRec::Clear()
{
  RecDeltaX.clear();
  RecDeltaY.clear();
  RecDircX.clear();
  RecDircY.clear();
  RecDircZ.clear();
  RecPropLength.clear();
  RecFlightTime.clear();
}
void DtofRec::Simu2Lab(double &x, double &y, double &z)
{
  double Sdata[] = {x, y, z};
  double Mdata[] = {0, 0, -1, 1, 0, 0, 0, -1, 0};
  TMatrixD M(3, 3, Mdata);
  TMatrixD S(1, 3, Sdata);
  TMatrixD L = S * M;
  x = TMatrixDRow(L, 0)(0);
  y = TMatrixDRow(L, 0)(1);
  z = TMatrixDRow(L, 0)(2);
}
double DtofRec::GetRemainder(double D, double d){
  int n = D/d;
  return D-n*d;
}
void DtofRec::InitialRefcounter()
{
  Dctr=0;
  Lctr=0;
  Rctr=0;
  Fctr=0;
}
void DtofRec::SetTrackHit()
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
  Px = FTOFdetx-QuartzX/2;
  Py = FTOFdety-QuartzX/2/Dx*Dy;
  Pz = FTOFdetz-QuartzX/2/Dx*Dz;
  double T0Px = T0detx -QuartzX/2;
  double T0Py = T0dety -QuartzX/2/Dx*Dy;
  double T0Pz = T0detz -QuartzX/2/Dx*Dz;
  FL = sqrt(Power(Px - T0Px, 2) + Power(Py - T0Py, 2) + Power(Pz - T0Pz, 2));

  Simu2Lab(Dx, Dy, Dz);
  Simu2Lab(Px, Py, Pz);
  double DieL = Sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
  Dx = Dx / DieL;
  Dy = Dy / DieL;
  Dz = Dz / DieL;
  //cout << "Track diraction (Dx,Dy,Dz) = (" << Dx << "," << Dy << "," << Dz << ")" << endl;
  // cout << "Track hit position (Px,Py) = (" << Px << "," << Py << ")" << endl;
  //cout << "FL="<<FL<<endl;
  InitialRefcounter();
}

void DtofRec::SetPhotonHit(int i)
{
  //fID = FTOFeleid[i];
  //T = FTOFelefittime1[i] - T0dett; // T= TOF(T0-FTOF) + TOP
  fID = FTOFphotonid[i];
  T = FTOFphotont[i] - T0dett; // T= TOF(T0-FTOF) + TOP
  ChX = fID / 4;
  ChY = fID % 4;
  //op<<"FID= "<<fID<<",("<<ChX<<","<<ChY<<")"<<endl;
  //T = GlobalTime->at(i)+TTS-T0;
  //cout << fID << ",ChX= " << ChX << ", ChY= " << ChY << endl;
}

void DtofRec::RecHitPos()
{

  int PhotonDetID = ChX / 4;
  double posX0 = R10754_sensor_L / 2 - R10754_sensor_edge_L - R10754_anode_interval_L - R10754_anode_L / 2 - (R10754_anode_interval_L + R10754_anode_L) * (ChX % 4);
  double posY0 = R10754_sensor_L / 2 - R10754_sensor_edge_L - R10754_anode_interval_L - R10754_anode_L / 2 - (R10754_anode_interval_L + R10754_anode_L) * ChY;

  //cout << "anode pos(x,y):" << posX0 << "\t" << posY0 << endl;
  HitX = posX0 + (sensorN / 2 - 0.5) * (R10754_sensor_L + R10754_interval_L) - (R10754_sensor_L + R10754_interval_L) * PhotonDetID;
  HitY = posY0 - QuartzZ / 2 + R10754_sensor_L / 2 + MCPPMToffset;
  op << "Hit pos (" << HitX << "," << HitY << ")" << endl;
  //cout << "Hit pos (" << HitX << "," << HitY << ")" << endl;
  /*
  double posX0 = (R + interval / Cos(Pi() / 2. - Pi() / SectorNu)) * Cos(Pi() / 4.);
  double posY0 = (R + interval / Cos(Pi() / 2. - Pi() / SectorNu)) * Sin(Pi() / 4.);
  double PDW = -(PhotonDetW + Gap) * (PhotonDetNu - 1) / 2. + (PhotonDetW + Gap) * PhotonDetID;
  HitX = posX0 + (PDW + CathodeW * (ChX % 4 - 1.5)) * Sin(Pi() / 4.) + CathodeW * (ChY - 1.5) * Cos(Pi() / 4.);
  HitY = posY0 - (PDW + CathodeW * (ChX % 4 - 1.5)) * Cos(Pi() / 4.) + CathodeW * (ChY - 1.5) * Sin(Pi() / 4.);
  */
}
void DtofRec::RecTOF(double beta)
{
  //MirrorTransformation(a, num);
  
  int i = RecDeltaX.size();
  double CosThetaC = 1. / Np / beta;
  RecDeltaX.push_back(HitX - Px);
  RecDeltaY.push_back(HitY - Py);
  RecDircX.push_back(RecDeltaX[i] / Abs(RecDeltaY[i]));
  RecDircY.push_back(RecDeltaY[i] / Abs(RecDeltaY[i]));
  RecDircZ.push_back(0);
  RecPropLength.push_back(0);
  RecFlightTime.push_back(0);

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
    double L_Z1 = Abs(RecDeltaY[i] / RecDircY[i]*Z1);
    double L_Z2 = Abs(RecDeltaY[i] / RecDircY[i]*Z2);
    int tag1 = 0, tag2 = 0;
    double costheta1 =Dx * RecDircX[i] + Dy * RecDircY[i] + Dz * Z1;
    double totalref1 = (Power(RecDircX[i], 2) + Power(RecDircY[i], 2)) / (Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Z1 * Z1)*NpMax*NpMax/0.95;
    double costheta2 =Dx * RecDircX[i] + Dy * RecDircY[i] + Dz * Z2;
    double totalref2 = (Power(RecDircX[i], 2) + Power(RecDircY[i], 2)) / (Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Z2 * Z2)*NpMax*NpMax/0.95;
    double R1 = GetRemainder(L_Z1,2*QuartzX);
    double R2 = GetRemainder(L_Z2,2*QuartzX);
    //op<<"costheta :"<<costheta1<<"\t"<<costheta2<<endl;
    //op<<"totalref :"<<totalref1<<"\t"<<totalref2<<endl;
    //op<<"Remainder :"<<R1<<"\t"<<R2<<endl;
    if (costheta1 > 0 
    &&  totalref1 >= 1 
    )
      tag1 = 1;
    if (costheta2 > 0 
    && totalref2 >= 1 
    )
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
      RecFlightTime[i] = T - RecPropLength[i] / TMath::C() * Ng * 1e6;
    }
  //op<<Dctr<<"\t"<<Lctr<<"\t"<<Rctr<<"\t"<<Fctr<<"\t"<<tag1<<"\t"<<RecDircX[i]<<"\t"<<RecDircY[i]<<"\t"<<Z1<<"\t"<<RecDeltaX[i]<<"\t"<<RecDeltaY[i]<<"\t"<<L_Z1<<"\t"<< RecPropLength[i]<<"\t"<< RecFlightTime[i]<<"\t"<<SimuFT<<"\t"<<HyFlightTime<<endl;
  //op<<Dctr<<"\t"<<Lctr<<"\t"<<Rctr<<"\t"<<Fctr<<"\t"<<tag2<<"\t"<<RecDircX[i]<<"\t"<<RecDircY[i]<<"\t"<<Z2<<"\t"<<RecDeltaX[i]<<"\t"<<RecDeltaY[i]<<"\t"<<L_Z2<<"\t"<< RecPropLength[i]<<"\t"<< RecFlightTime[i]<<"\t"<<SimuFT<<"\t"<<HyFlightTime<<endl;
  }
  //cout << "Rec photon Diraction=(" << RecDircX.back() << "," << RecDircY.back() << "," << RecDircZ.back() << ")" << endl
  //<< "RecPropLength=" << RecPropLength.back() << endl
  //<< "Rec TOP =" <<RecPropLength.back() / TMath::C() * Ng * 1e6<< endl;

  //if (num % 2 == 1)
  //  MirrorTransformation(a, num);
  //else
  //  MirrorTransformation(Mirror(a), num);
}

void DtofRec::Reconstruction(double mass)
{
  Clear();
  double beta = 1. / Sqrt(1 + Power(mass / PrimaryMomentum, 2));
  HyFlightTime = FL / TMath::C() / beta * 1e6;
  //cout << "[+] | - Start Reconstruction " << endl;
  //cout<< "beta= "<<beta<<endl;
  RecHitPos();
  //RecTOF(beta);
  

    for (int j = 0; j < 3; j++)
    {
      for (int s = 0; s < 2; s++)
      {
        MirrorTransformation((enum lightpath)j, s);
        RecTOF(beta);
        SetTrackHit();
      }
    }

    // mirror 2 time
    for (int j = 1; j < 3; j++)
    {
      for (int s = 0; s < 2; s++)
      {
        MirrorTransformation((enum lightpath)j, 0);
        MirrorTransformation(Mirror((enum lightpath)j), s);
        RecTOF(beta);
        SetTrackHit();
      }
    }
  

  int tag = 0;
  for (int id = 0; id < (int)RecFlightTime.size(); id++)
  {
    //if(Abs(RecFlightTime[id]-HyFlightTime) < Abs(RecFlightTime[tag]-HyFlightTime)) tag = id;
    if (Abs(RecFlightTime[id] - HyFlightTime) < Abs(RecFlightTime[tag] - HyFlightTime))
      tag = id;
  }
  BestPropLength = RecPropLength[tag];
  BestFlightTime = RecFlightTime[tag];
  //op<< BestPropLength<<"\t"<< BestFlightTime<<"\t"<<SimuFT<<"\t"<<HyFlightTime<<endl;

  //cout << "BestPropLength= " << BestPropLength<< endl
  //     << "BestFlightTime= " << BestFlightTime << endl;
}

void DtofRec::Loop()
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
  for (Long64_t tentry = 0; tentry < nentries; tentry++)
  {
    fChain->GetEntry(tentry);
    if (CRRBtheta < 0)
      continue;
    
    TOFvec.clear();
    LOPvec.clear();
    TOPvec.clear();
    chidvec.clear();
    trueTOPvec.clear();
    Xvec.clear();
    Yvec.clear();
    SimuFT = FTOFdett - T0dett;
    SetTrackHit();
    /*
    for (int i = 0; i < 16 * sensorN; i++)
    {

      if (FTOFeletot[i] < 0)
        continue;
      SetPhotonHit(i);
      Reconstruction(Mass_Mu / 1e3);
      FT[i] = BestFlightTime;
      TL[i] = BestPropLength;
      TOP[i] = BestPropLength / TMath::C() * Ng * 1e6;
      TOFvec.push_back(BestFlightTime);
      LOPvec.push_back(BestPropLength);
      chidvec.push_back(FTOFelechid[i]);
    };
    */
    int N = FTOFphotonid.size();
      //op<<"Entry$= "<<tentry<<endl;
    for (int i = 0; i < N; i++)
    {
      //if (FTOFeletot[i] < 0)
      //  continue;
        
      SetPhotonHit(i);
      Reconstruction(Mass_Mu / 1e3);
      TOFvec.push_back(BestFlightTime);
      LOPvec.push_back(BestPropLength);
      TOPvec.push_back(BestPropLength / TMath::C() * Ng * 1e6);
      trueTOPvec.push_back(FTOFphotonTOP[i]);
      chidvec.push_back(fID);
      Xvec.push_back(ChX);
      Yvec.push_back(ChY);
    };
    int cut = 0;
    if (TOFvec.size() > cut)
    {
      sort(TOFvec.begin(), TOFvec.end());
      double sum = 0, count = 0;
      double sum2 = 0;
      for (int i = cut / 2; i < TOFvec.size() - cut / 2; i++)
      {
        if (TOPvec.at(i) <=0.1)
          continue;
        if(abs(TOFvec.at(i)-HyFlightTime)>0.3) continue;
        sum += TOFvec.at(i);
        count += 1;
        sum2 += LOPvec.at(i);
      }
      if (count != 0)
      {
        MeanFT = sum / count;
        MeanTL = sum2 / count;
        TT.push_back(MeanFT);
        AA.push_back(MeanTL);
        ot->Fill();
        exctr++;
      }
    }
  }
  fout->WriteTObject(ot);
}

#endif
