/*
*   File : CRTest/inc/CryPositionSD.hh
*
*   Brief: Postion-sensitive Detector for Cosmic-ray or others
*
*   Description:
*       Derive from G4LogicalVolume
*       Data/Physical Member // G4LogicalVolume(s)
*           TotalBox - Hold all components
*           Components - Detector detail structure
*       Inner Class?
*           CryHit, CrySD
*/

#ifndef CRTest_CryPositionSD_h
#define CRTest_CryPositionSD_h

#include "VirtualSD.hh"

#include "CryHit.hh"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"

#include "G4String.hh"

class CryPositionSD : public VirtualSD {

public:
    CryPositionSD(G4String& name);
    virtual ~CryPositionSD();
public:
    virtual void Initialize(G4HCofThisEvent* hce);
    virtual G4bool ProcessHits(
		G4Step* aStep, G4TouchableHistory* roHist);
	virtual void EndOfEvent(G4HCofThisEvent* hce);

	virtual void CreateEntry(
		G4int ntupleID, G4RootAnalysisManager* rootData);
	virtual void FillEntry(
		G4int ntupleID, G4RootAnalysisManager* rootData);

private:
    G4int fHCID;
    CryHC* fHC;

	G4int fFirstColID;
	std::vector<double>* fEdep;
	std::vector<int>* fCounter;

	G4String fname;
	G4String buff;
};

#endif /*CRTest_CryPositionSD_h*/