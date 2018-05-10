/*
*   File : CRTest/inc/OpRecorder.hh
*   Brief: Record optical process information
*   Description:
*       Record optical photon : count, energy, direction
*       Scintillation, ScintToFiber, WLS, PMT
*/

#ifndef CRTest_OpRecorder_h
#define CRTest_OpRecorder_h

#include "VirtualRecorder.hh"

#include "G4RootAnalysisManager.hh"
#include "G4Track.hh"
#include "globals.hh"
#include<vector>

enum OpPhotonType{
	Nothing = 0,
	Quartz2Air,
	Quartz2GlueL,
	Quartz2GlueR,
	Glue2PmtL,
	Glue2PmtR,
	GlueRef,
	Detected
};

class OpRecorder : public VirtualRecorder{
public:
    OpRecorder();
    ~OpRecorder();
    static OpRecorder* Instance();
    void Reset();
    void Print();
	void SetBoundaryName(G4String);

public: // for class Analysis
	virtual void CreateEntry(
		G4int ntupleID, G4RootAnalysisManager*);
	virtual void FillEntry(
		G4int ntupleID, G4RootAnalysisManager*);

	G4bool Record(const G4Track*);

// TODO : Convert to 'private' and add Increment method
public:
    G4int nCerenkov;
    G4int nScintTotal;
    G4int nQuartz2Air;
	G4int nQuartz2GlueL;
	G4int nQuartz2GlueR;
    G4int nWlsEmit;
    G4int nGlue2PMTL;
	G4int nGlue2PMTR;
	G4int nDetection; // Detected by PmtSD
// For deatil probe
    
    G4int nBoundaryRefraction;
	G4int nBoundaryAbsorption;
    G4int nBoundaryTransmission;
    G4int nBoundaryUndefined;
    G4int nBoundaryWARNNING;
	G4int nFresnelReflection;
	G4int nTotalInternalReflection;
	G4int nLambertianReflection;
	G4int nLobeReflection;
	G4int nSpikeReflection;
	G4int nBackScattering;

    G4int nDebug; // use for debug

	std::vector<int>* fID;
	std::vector<double>* fL;
	std::vector<int>* fBounce;
	std::vector<double>* fWaveL;

private:
    static OpRecorder* fgInstance;
	G4String boundaryName;

	std::vector<double>* fCount;
	std::vector<double>* fEk;
	std::vector<double>* fTime;
	std::vector<double>* fX;
	std::vector<double>* fY;
	std::vector<double>* fZ;
	std::vector<double>* fPX;
	std::vector<double>* fPY;
	std::vector<double>* fPZ;

	
};


#endif // CRTest_OpRecorder_h