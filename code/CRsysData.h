#ifndef CRsysData_h
#define CRsysData_h 1
#include <iostream>
#include "TROOT.h"
;
class CRsysData
{
public:
    CRsysData()
    {
        Initial();
    }
    ~CRsysData() {}

    // for T0
    // *** T0 photon collection info
    vector<int> T0photonid; // channel data
    vector<double> T0photonE;
    vector<double> T0photonTOP;
    vector<double> T0photonx;
    vector<double> T0photony;
    vector<double> T0photonz;
    vector<double> T0photont;
    

    // *** T0 hitted by muon
    int T0detid; // detector id
    double T0detx;
    double T0dety;
    double T0detz;
    double T0dett;

    // for FTOF
    // *** FTOF photon collection info
    vector<int> FTOFphotonid; // channel data
    vector<double> FTOFphotonE;
    vector<double> FTOFphotonTOP;
    vector<double> FTOFphotonx;
    vector<double> FTOFphotony;
    vector<double> FTOFphotonz;
    vector<double> FTOFphotont;

    // *** FTOF hitted by muon
    int FTOFdetid; // detector id
    double FTOFdetx;
    double FTOFdety;
    double FTOFdetz;
    double FTOFdett;

    // for MM
    double Trackerdetx[4];
    double Trackerdety[4];
    double Trackerdetz[4];
    double Trackerdett[4];

    //for CR
    double CRE;
    double CRpx;
    double CRpy;
    double CRpz;
    double CRtheta;
    double CRphi;

    void Initial()
    {
        T0detid = -999; // detector id
        T0detx = -999;
        T0dety = -999;
        T0detz = -999;
        T0dett = -999;
        //for CR
        CRE = -999;
        CRpx = -999;
        CRpy = -999;
        CRpz = -999;
        CRtheta = -999;
        FTOFdetid = -999; // detector id
        FTOFdetx = -999;
        FTOFdety = -999;
        FTOFdetz = -999;
        FTOFdett = -999;
        for (int i = 0; i < 4; i++)
        {

            Trackerdetx[i] = -999;
            Trackerdety[i] = -999;
            Trackerdetz[i] = -999;
            Trackerdett[i] = -999;
        }
        vector<int>().swap(FTOFphotonid); // channel data
        vector<int>().swap(T0photonid);   // channel data
        vector<double>().swap(FTOFphotonE);
        vector<double>().swap(FTOFphotonTOP);
        vector<double>().swap(FTOFphotonx);
        vector<double>().swap(FTOFphotony);
        vector<double>().swap(FTOFphotonz);
        vector<double>().swap(FTOFphotont);
        vector<double>().swap(T0photonE);
        vector<double>().swap(T0photonTOP);
        vector<double>().swap(T0photonx);
        vector<double>().swap(T0photony);
        vector<double>().swap(T0photonz);
        vector<double>().swap(T0photont);
    }
};

#endif