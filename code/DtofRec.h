#ifndef DtofRec_h
#define DtofRec_h

#include "CRdata.h"
#include "LUT.h"
#include "TMath.h"
#include "TRandom3.h"
class DtofRec:public CRdata{
  public :
    DtofRec(string in="data.root",
        string out="Rec.root",
        TTree *tree=0);
    ~DtofRec();

    enum lightpath {direct, left, right};
    lightpath Mirror(lightpath);

    /// ReConstruction
    void MirrorTransformation(lightpath, int);
    void RecHitPos(); 
    void RecTOF(lightpath,double,int); 
    void Reconstruction(double);

    /// Set Hit 
    void Clear();
    void SetTrackHit();
    void SetPhotonHit(int i);

    void Loop();
    /// LUT
    LUT LUT_n;

    /// Hit
    double Dx;
    double Dy;
    double Dz;
    double Px;
    double Py;
    double FL;
    int    ChX;
    int    ChY;
    double T;
    double TTS,T0;

    /// result of Reconstruction
    double HitX,HitY;
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


    /// parameters
    double A,B,C;
    double NpMax,Np,Ng;
    double interval,PhotonDetW,CathodeW; // mm
    double R,Gap,QuartzR1,QuartzR2; // mm
    int    SectorNu, PhotonDetNu;
    int MaxReflection;
    double PrimaryMomentum;

    /// random number
    TRandom3* fRandom;
};

#endif
