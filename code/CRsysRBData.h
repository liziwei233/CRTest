#ifndef CRsysRBData_h
#define CRsysRBData_h 1
#include <iostream>
#include "TROOT.h"
#include "CRsysData.h"
;
class CRsysRBData : public CRsysData
{
public:
    CRsysRBData()
    {
        RBInitial();
    }
    ~CRsysRBData() {}

    //CRsysData *p;
    // *** T0 hitted by muon
    double T0detRBx;
    double T0detRBy;
    double T0detRBz;

    // *** T0 rebuild signal
    vector<int> T0eleid; // mV
    vector<double> T0elethrd; // mV
    vector<double> T0eleU;    //amplitude
    vector<double> T0eletot;
    vector<double> T0elethtime[2]; // 0-of leading edge, 1-of falling edge
    vector<double> T0elefittot;
    vector<double> T0elefittime[2];

    // *** FTOF hitted by muon
    double FTOFdetRBx;
    double FTOFdetRBy;
    double FTOFdetRBz;

    // *** FTOF rebuild signal
    vector<int> FTOFeleid; // mV
    vector<double> FTOFelethrd; // mV
    vector<double> FTOFeleU;    //amplitude
    vector<double> FTOFeletot;
    vector<double> FTOFelethtime[2]; // 0-of leading edge, 1-of falling edge
    vector<double> FTOFelefittot;
    vector<double> FTOFelefittime[2];

    //for CR

    double CRRBpx;
    double CRRBpy;
    double CRRBpz;
    double CRRBtheta;

    void RBInitial()
    {
        //p->Initial();
        this->Initial();

        T0detRBx = -999;
        T0detRBy = -999;
        T0detRBz = -999;

        // *** T0 rebuild signal
    vector<int>().swap(T0eleid); // mV
    vector<double>().swap(T0elethrd); // mV
    vector<double>().swap(T0eleU);    //amplitude
    vector<double>().swap(T0eletot);
    vector<double>().swap(T0elethtime[0]); // 0-of leading edge, 1-of falling edge
    vector<double>().swap(T0elethtime[1]); // 0-of leading edge, 1-of falling edge
    vector<double>().swap(T0elefittot);
    vector<double>().swap(T0elefittime[0]);
    vector<double>().swap(T0elefittime[1]);

        

        FTOFdetRBx = -999;
        FTOFdetRBy = -999;
        FTOFdetRBz = -999;

             // *** FTOF rebuild signal
    vector<int>().swap(FTOFeleid); // mV
    vector<double>().swap(FTOFelethrd); // mV
    vector<double>().swap(FTOFeleU);    //amplitude
    vector<double>().swap(FTOFeletot);
    vector<double>().swap(FTOFelethtime[0]); // 0-of leading edge, 1-of falling edge
    vector<double>().swap(FTOFelethtime[1]); // 0-of leading edge, 1-of falling edge
    vector<double>().swap(FTOFelefittot);
    vector<double>().swap(FTOFelefittime[0]);
    vector<double>().swap(FTOFelefittime[1]);

        //for CR

        CRRBpx = -999;
        CRRBpy = -999;
        CRRBpz = -999;
        CRRBtheta = -999;
    }
};

#endif