#include "DircRecSvc.h"

using namespace TMath;
using namespace std;

DircRecSvc::DircRecSvc(string in, string out, TTree* tree):DircData(in,out,tree),
  LUT_n("LUT/n.txt",400,2)
{
  fRandom = new TRandom3(time(0));
  NpMax = 1.48779, Np = 1.46979, Ng = 1.51404;
  MaxReflection = 2;

  SectorNu = 12;
  interval=10, PhotonDetW=27.5, CathodeW=5.5 ; // mm
  double angle = Pi()/2.-Pi()/SectorNu;
  QuartzR1 = 560./Sin(angle) - interval/Cos(angle)*Tan(angle/2);
  QuartzR2 = 1050. - interval/Cos(angle)/Tan(angle/2);
  R = QuartzR2*Sin(angle)-PhotonDetW/2.;
  PhotonDetNu = int((R/Tan(angle)*2.-PhotonDetW/Tan(angle))/PhotonDetW);
  Gap = (R/Tan(angle)*2.-PhotonDetW/Tan(angle)-PhotonDetW*PhotonDetNu)/PhotonDetNu;
}

DircRecSvc::~DircRecSvc()
{
}

DircRecSvc::lightpath DircRecSvc::Mirror(DircRecSvc::lightpath a)
{
  if(a == left) return right;
  if(a == right) return left;
  if(a == direct) return direct;
}

void DircRecSvc::MirrorTransformation(DircRecSvc::lightpath a, int i)
{
  if(a == direct) return;
  if(a == left)
  {
    double angle,angleD, x1,y1,x2,y2, tmpX,tmpY;
    angle = Pi()/2.-Pi()/SectorNu;
    x1 = (QuartzR1*Sin(angle)+interval/Cos(angle))*Cos(Pi()/4.)-QuartzR1*Cos(angle)*Sin(Pi()/4.);
    y1 = (QuartzR1*Sin(angle)+interval/Cos(angle))*Sin(Pi()/4.)+QuartzR1*Cos(angle)*Cos(Pi()/4.);
    x2 = (QuartzR2*Sin(angle)+interval/Cos(angle))*Cos(Pi()/4.)-QuartzR2*Cos(angle)*Sin(Pi()/4.);
    y2 = (QuartzR2*Sin(angle)+interval/Cos(angle))*Sin(Pi()/4.)+QuartzR2*Cos(angle)*Cos(Pi()/4.);
    A = y1-y2;
    B = x2-x1;
    C = x1*y2-x2*y1;
    tmpX = ((B*B-A*A)*Px-2.*A*B*Py-2.*A*C)/(A*A+B*B);
    tmpY = ((A*A-B*B)*Py-2.*A*B*Px-2.*B*C)/(A*A+B*B);
    Px = tmpX;
    Py = tmpY;

    angle = Pi()/4.+Pi()/SectorNu;
    if(Dx>0) angleD = ATan(Dy/Dx);
    if(Dx<0) angleD = Pi()+ATan(Dy/Dx);
    tmpX = Sqrt(Dx*Dx+Dy*Dy)*Cos(2.*angle-angleD);
    tmpY = Sqrt(Dx*Dx+Dy*Dy)*Sin(2.*angle-angleD);
    Dx = tmpX;
    Dy = tmpY;
  }
  if(a == right)
  {    
    double angle,angleD, x1,y1,x2,y2, tmpX,tmpY;
    angle = Pi()/2.-Pi()/SectorNu;
    x1 = (QuartzR1*Sin(angle)+interval/Cos(angle))*Cos(Pi()/4.)+QuartzR1*Cos(angle)*Sin(Pi()/4.);
    y1 = (QuartzR1*Sin(angle)+interval/Cos(angle))*Sin(Pi()/4.)-QuartzR1*Cos(angle)*Cos(Pi()/4.);
    x2 = (QuartzR2*Sin(angle)+interval/Cos(angle))*Cos(Pi()/4.)+QuartzR2*Cos(angle)*Sin(Pi()/4.);
    y2 = (QuartzR2*Sin(angle)+interval/Cos(angle))*Sin(Pi()/4.)-QuartzR2*Cos(angle)*Cos(Pi()/4.);
    A = y1-y2;
    B = x2-x1;
    C = x1*y2-x2*y1;
    tmpX = ((B*B-A*A)*Px-2.*A*B*Py-2.*A*C)/(A*A+B*B);
    tmpY = ((A*A-B*B)*Py-2.*A*B*Px-2.*B*C)/(A*A+B*B);
    Px = tmpX;
    Py = tmpY;

    angle = Pi()/4.-Pi()/SectorNu;
    if(Dx>0) angleD = ATan(Dy/Dx);
    if(Dx<0) angleD = Pi()+ATan(Dy/Dx);
    tmpX = Sqrt(Dx*Dx+Dy*Dy)*Cos(2.*angle-angleD);
    tmpY = Sqrt(Dx*Dx+Dy*Dy)*Sin(2.*angle-angleD);
    Dx = tmpX;
    Dy = tmpY;
  }
  if(i>1)
  {
    MirrorTransformation(Mirror(a),i-1);
  }
}

void DircRecSvc::Clear()
{
  RecDeltaX.clear();
  RecDeltaY.clear();
  RecDircX.clear();
  RecDircY.clear();
  RecDircZ.clear();
  RecPropLength.clear();
  RecFlightTime.clear();
}

void DircRecSvc::SetTrackHit()
{
  T0 = fRandom->Gaus(0,0.04);
  Dx = TrackDx->at(0);
  Dy = TrackDy->at(0);
  Dz = TrackDz->at(0);
  Px = TrackPx->at(0);
  Py = TrackPy->at(0);
  FL = FlightLength->at(0);
  double DieL = Sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
  Dx = Dx/DieL;
  Dy = Dy/DieL;
  Dz = Dz/DieL;
}

void DircRecSvc::SetPhotonHit(int i)
{
  TTS = fRandom->Gaus(0,0.07);
  ChX = ChannelX->at(i);
  ChY = ChannelY->at(i);
  T = GlobalTime->at(i)+TTS-T0;
}

void DircRecSvc::RecHitPos()
{
  int PhotonDetID = ChX/4;
  double posX0 = (R + interval/Cos(Pi()/2.-Pi()/SectorNu))*Cos(Pi()/4.);
  double posY0 = (R + interval/Cos(Pi()/2.-Pi()/SectorNu))*Sin(Pi()/4.);
  double PDW = -(PhotonDetW+Gap)*(PhotonDetNu-1)/2.+(PhotonDetW+Gap)*PhotonDetID;
  HitX = posX0 + (PDW+CathodeW*(ChX%4-1.5))*Sin(Pi()/4.) + CathodeW*(ChY-1.5)*Cos(Pi()/4.);
  HitY = posY0 - (PDW+CathodeW*(ChX%4-1.5))*Cos(Pi()/4.) + CathodeW*(ChY-1.5)*Sin(Pi()/4.);
}
void DircRecSvc::RecTOF(lightpath a, double beta, int num=1)
{
  MirrorTransformation(a, num);

  int i = RecDeltaX.size();
  double CosThetaC = 1./Np/beta;
  RecDeltaX.push_back(HitX-Px);
  RecDeltaY.push_back(HitY-Py);
  RecDircX.push_back(RecDeltaX[i]/Abs(RecDeltaY[i]));
  RecDircY.push_back(RecDeltaY[i]/Abs(RecDeltaY[i]));
  RecDircZ.push_back(0);
  RecPropLength.push_back(0);
  RecFlightTime.push_back(0);

  A = Dz*Dz-CosThetaC*CosThetaC;
  B = 2.*Dz*(Dx*RecDircX[i]+Dy*RecDircY[i]);
  C = Power(Dx*RecDircX[i]+Dy*RecDircY[i],2)-CosThetaC*CosThetaC*(Power(RecDircX[i],2)+Power(RecDircY[i],2));

  double Delta = B*B-4.*A*C;
  if(Delta >= 0)
  {
    double Z1 = (-B+Sqrt(Delta))/2./A, Z2 = (-B-Sqrt(Delta))/2./A;
    int tag1=0, tag2=0;
    if(Dx*RecDircX[i]+Dy*RecDircY[i]+Dz*Z1>0&&
        (Power(RecDircX[i],2)+Power(RecDircY[i],2))/(Power(RecDircX[i],2)+Power(RecDircY[i],2)+Z1*Z1)>=0.95/NpMax/NpMax) tag1=1;
    if(Dx*RecDircX[i]+Dy*RecDircY[i]+Dz*Z2>0&&
        (Power(RecDircX[i],2)+Power(RecDircY[i],2))/(Power(RecDircX[i],2)+Power(RecDircY[i],2)+Z2*Z2)>=0.95/NpMax/NpMax) tag2=1;
    if(tag1 || tag2)
    {
      if(tag1 && tag2)  RecDircZ[i] = Abs(Z1)<Abs(Z2)? Z1:Z2;
      if(tag1 && !tag2) RecDircZ[i] = Z1;
      if(tag2 && !tag1) RecDircZ[i] = Z2;
      RecPropLength[i] = Abs(RecDeltaY[i]*Sqrt(Power(RecDircX[i],2)+Power(RecDircY[i],2)+Power(RecDircZ[i],2))/RecDircY[i]);
      RecFlightTime[i] = T - RecPropLength[i]/TMath::C()*Ng*1e6;
    }
  }

  if(num%2 == 1) MirrorTransformation(a,num);
  else MirrorTransformation(Mirror(a),num);
}

void DircRecSvc::Reconstruction(double mass)
{
  Clear();
  double beta = 1./Sqrt(1+Power(mass/PrimaryMomentum,2));
  double HyFlightTime = FL/TMath::C()/beta*1e6;

  RecHitPos();
  RecTOF(direct,beta);
  for(int i=1;i<=MaxReflection;i++)
  {
    RecTOF(left,beta,i);
    RecTOF(right,beta,i);
  }

  int tag = 0;
  for(int id = 0;id<(int)RecFlightTime.size(); id++)
  {
    if(Abs(RecFlightTime[id]-HyFlightTime) < Abs(RecFlightTime[tag]-HyFlightTime)) tag = id;
  }
  BestPropLength = RecPropLength[tag];
  BestFlightTime = RecFlightTime[tag];
}
