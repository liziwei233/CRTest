#ifndef DircRecAlg_h
#define DircRecAlg_h

#include "DircData.h"
#include "LUT.h"
#include "DircRecSvc.h"

class DircRecAlg:public DircRecSvc{
  public :
    DircRecAlg(string in="data/DIRC.root",
        string out="result/DircRec.root",
        TTree *tree=0);
    ~DircRecAlg();
    virtual void Loop();
    void book();

    double TruThetaC,  TruPropLength;
    double RecPropLengthHyPi,  RecFlightTimeHyPi;
    double RecPropLengthHyK,  RecFlightTimeHyK;
    double ProbL;

    std::map<string,TH1F*> hist1;
    std::map<string,TH2F*> hist2;

    std::map<string,double> mass;
    std::map<string,TF1*> PDF;
    
    std::vector<double> TOFHyPi, TOFHyK;

    TTree* fNtuple;
};

#endif
