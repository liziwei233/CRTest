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
    int T0eleid[4];      // mV
    double T0elethrd[4]; // mV
    double T0eleU[4];    //amplitude
    double T0eletot[4];
    double T0elethtime1[4]; // 0-of leading edge, 1-of falling edge
    double T0elethtime2[4]; // 0-of leading edge, 1-of falling edge
    double T0elefittot[4];
    double T0elefittime1[4];
    double T0elefittime2[4];

    // *** FTOF hitted by muon
    double FTOFdetRBx;
    double FTOFdetRBy;
    double FTOFdetRBz;

    // *** FTOF rebuild signal
    int FTOFeleid[128];      // mV
    double FTOFelethrd[128]; // mV
    double FTOFeleU[128];    //amplitude
    double FTOFeletot[128];
    double FTOFelethtime1[128]; // 0-of leading edge, 1-of falling edge
    double FTOFelethtime2[128]; // 0-of leading edge, 1-of falling edge
    double FTOFelefittot[128];
    double FTOFelefittime1[128];
    double FTOFelefittime2[128];

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
        for (int i = 0; i < 4; i++)
        {
            T0eleid[i] = -999;   // mV
            T0elethrd[i] = -999; // mV
            T0eleU[i] = -999;    //amplitude
            T0eletot[i] = -999;
            T0elethtime1[i] = -999; // 0-of leading edge, 1-of falling edge
            T0elethtime2[i] = -999; // 0-of leading edge, 1-of falling edge
            T0elefittot[i] = -999;
            T0elefittime1[i] = -999;
            T0elefittime2[i] = -999;
        }

        for (int i = 0; i < 128; i++)
        {

            FTOFeleid[i] = -999;   // mV
            FTOFelethrd[i] = -999; // mV
            FTOFeleU[i] = -999;    //amplitude
            FTOFeletot[i] = -999;
            FTOFelethtime1[i] = -999; // 0-of leading edge, 1-of falling edge
            FTOFelethtime2[i] = -999; // 0-of leading edge, 1-of falling edge
            FTOFelefittot[i] = -999;
            FTOFelefittime1[i] = -999;
            FTOFelefittime2[i] = -999;
        }

        FTOFdetRBx = -999;
        FTOFdetRBy = -999;
        FTOFdetRBz = -999;

        //for CR

        CRRBpx = -999;
        CRRBpy = -999;
        CRRBpz = -999;
        CRRBtheta = -999;
    }
};

#endif