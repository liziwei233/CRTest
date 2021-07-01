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
  enum lightpath
  {
    direct,
    left,
    right,
    front
  };
  lightpath Mirror(lightpath);

  /// ReConstruction
  void MirrorTransformation(lightpath, int, int);
  void RecHitPos();
  void RecTOF(double);
  void Reconstruction(double);

  /// Set Hit
  void Clear();
  void SetTrackHit();
  void SetPhotonHit(int i);
  void Simu2Lab(double &x, double &y, double &z);

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
  double BestPropLength;
  double BestFlightTime;
  vector<double> TT;
  vector<double> AA;

  double FT[128];
  double TL[128];
  double TOP[128];
  double HitXpos[128];
  double HitYpos[128];
  double SimuFT;
  double MeanFT;
  double MeanTL;

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
  ot = new TTree("data", "restore Dtof rebuilded data  from G4");
  ot->Branch("FlightTime", FT, "FlightTime[128]/D");
  ot->Branch("MeanFT", &MeanFT, "MeanFT/D");
  ot->Branch("TrackLength", TL, "TrackLength[128]/D");
  ot->Branch("TOP", TOP, "TOP[128]/D");
  ot->Branch("MeanTL", &MeanTL, "MeanTL/D");
  ot->Branch("HitX", HitXpos, "HitX[128]/D");
  ot->Branch("HitY", HitYpos, "HitY[128]/D");
  ot->Branch("SimuFT", &SimuFT, "SimuFT/D");
  ot->Branch("CRRBtheta", &CRRBtheta);
  ot->Branch("CRE", &CRE);
  ot->Branch("FTOFphotonid", &FTOFphotonid);
  ot->Branch("FTOFphotonTOP", &FTOFphotonTOP);
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

void DtofRec::MirrorTransformation(DtofRec::lightpath a, int i, int s)
{
  if (a == direct)
    return;
  if (a == left)
  {
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
  if (i > 1)
  {
    int lp = (a + 1 + s) % 4;
    if (lp == 0)
      lp++;
    MirrorTransformation((enum lightpath)lp, i - 1, s);
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
  Simu2Lab(Dx, Dy, Dz);
  Px = FTOFdetx;
  Py = FTOFdety;
  Pz = FTOFdetz;
  Simu2Lab(Px, Py, Pz);

  FL = sqrt(Power(FTOFdetRBx - T0detRBx, 2) + Power(FTOFdetRBy - T0detRBy, 2) + Power(FTOFdetRBz - T0detRBz, 2));
  double DieL = Sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
  Dx = Dx / DieL;
  Dy = Dy / DieL;
  Dz = Dz / DieL;
  //cout << "Track diraction (Dx,Dy,Dz) = (" << Dx << "," << Dy << "," << Dz << ")" << endl;
  // cout << "Track hit position (Px,Py) = (" << Px << "," << Py << ")" << endl;
  //cout << "FL="<<FL<<endl;
}

void DtofRec::SetPhotonHit(int i)
{
  //TTS = fRandom->Gaus(0,0.07);
  fID = FTOFeleid[i];
  ChX = fID / 4;
  ChY = fID % 4;
  T = FTOFelefittime1[i] - T0dett; // T= TOF(T0-FTOF) + TOP
  //ChX = ChannelX->at(i);
  //ChY = ChannelY->at(i);
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
  //cout << "Hit pos (" << HitX << "," << HitY << ")" << endl;
  HitXpos[fID] = HitX;
  HitYpos[fID] = HitY;
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
    double Z1 = (-B + Sqrt(Delta)) / 2. / A, Z2 = (-B - Sqrt(Delta)) / 2. / A;
    int tag1 = 0, tag2 = 0;
    if (Dx * RecDircX[i] + Dy * RecDircY[i] + Dz * Z1 < 0 &&
        (Power(RecDircX[i], 2) + Power(RecDircY[i], 2)) / (Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Z1 * Z1) >= 0.95 / NpMax / NpMax)
      tag1 = 1;
    if (Dx * RecDircX[i] + Dy * RecDircY[i] + Dz * Z2 < 0 &&
        (Power(RecDircX[i], 2) + Power(RecDircY[i], 2)) / (Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Z2 * Z2) >= 0.95 / NpMax / NpMax)
      tag2 = 1;
    if (tag1 || tag2)
    {
      if (tag1 && tag2)
        RecDircZ[i] = Abs(Z1) < Abs(Z2) ? Z1 : Z2;
      if (tag1 && !tag2)
        RecDircZ[i] = Z1;
      if (tag2 && !tag1)
        RecDircZ[i] = Z2;
      RecPropLength[i] = Abs(RecDeltaY[i] * Sqrt(Power(RecDircX[i], 2) + Power(RecDircY[i], 2) + Power(RecDircZ[i], 2)) / RecDircY[i]);
      RecFlightTime[i] = T - RecPropLength[i] / TMath::C() * Ng * 1e6;
    }
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
  double HyFlightTime = FL / TMath::C() / beta * 1e6;
  //cout << "[+] | - Start Reconstruction " << endl;
  //cout<< "beta= "<<beta<<endl;
  RecHitPos();
  RecTOF(beta);
  for (int i = 1; i <= MaxReflection; i++)
  {
    for (int j = 1; j <= 3; j++)
    {
      for (int s = 0; s < i; s++)
      {

        MirrorTransformation((enum lightpath)j, i, s);
        RecTOF(beta);
        SetTrackHit();
      }
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
    int ctr = 0;
    memset(FT, 0, sizeof(FT));
    memset(TL, 0, sizeof(TL));
    memset(HitXpos, 0, sizeof(HitXpos));
    memset(HitYpos, 0, sizeof(HitYpos));
    SimuFT = FTOFdett - T0dett;
    SetTrackHit();
    for (int i = 0; i < 16 * sensorN; i++)
    {

      if (FTOFeletot[i] < 0)
        continue;
      SetPhotonHit(i);
      Reconstruction(Mass_Mu / 1e3);
      FT[i] = BestFlightTime;
      TL[i] = BestPropLength;
      TOP[i] = BestPropLength / TMath::C() * Ng * 1e6;
      ctr++;
    };
    if (ctr > 0)
    {
      MeanFT = Mean(16 * sensorN, FT) * 16 * sensorN / ctr;
      MeanTL = Mean(16 * sensorN, TL) * 16 * sensorN / ctr;
      TT.push_back(MeanFT);
      AA.push_back(MeanTL);
      ot->Fill();
      exctr++;
    }
  }
  fout->WriteTObject(ot);
}

#endif
