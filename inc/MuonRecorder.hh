/*
*   File : CRTest/inc/MuonRecorder.hh
*   Brief: Record Muon tracking information
*   Description:
*       Record primary muon
*		energy, time, position, direction
*/

#ifndef CRTest_MuonRecorder_h
#define CRTest_MuonRecorder_h

#include "VirtualRecorder.hh"

#include "G4RootAnalysisManager.hh"

#include "G4Track.hh"

#include "globals.hh"
#include <string>
#include <map>
using namespace std;
class MuonRecorder : public VirtualRecorder {

public:
	MuonRecorder();
	virtual ~MuonRecorder();

	static MuonRecorder* Instance();
	void Reset();

public:
	virtual void CreateEntry(
		G4int ntupleID, G4RootAnalysisManager*);
	virtual void FillEntry(
		G4int ntupleID, G4RootAnalysisManager*);
	
	G4bool Record(const G4Track*);
	std::vector<int>* fDetID;
	std::map<int, bool> flag; 
	std::vector<string>* fName;
	std::vector<int>* fID;

private:
	static MuonRecorder* fgInstance;

private: // Muon Track Information
	std::vector<int>* fCount;
	std::vector<double>* fEk;
	std::vector<double>* fTime;
	std::vector<double>* fX;
	std::vector<double>* fY;
	std::vector<double>* fZ;
	std::vector<double>* fPX;
	std::vector<double>* fPY;
	std::vector<double>* fPZ;
	
};

#endif