/**
*   FILE : CRTest/src/CryPositionSD.cc
*   Brief: Implementation of class CryPositionSD
*/

#include "CryPositionSD.hh"

#include "G4SDManager.hh"

#include "G4String.hh"
#include "G4ios.hh"

CryPositionSD::CryPositionSD(G4String &name)
	: VirtualSD(name), fHCID(0), fFirstColID(-1),
	fEdep(NULL),fCounter(NULL)
{
    collectionName.insert("CryHC");

	fEdep = new std::vector<double>;
	fCounter = new std::vector<int>;
}

CryPositionSD::~CryPositionSD()
{
	delete fEdep;
	delete fCounter;
}

void CryPositionSD::Initialize(G4HCofThisEvent *hce)
{
    G4cout << "[-] INFO - CryPositionSD Initialized."
           << " - by CryPositionSD" << G4endl;
    fHC = new CryHC(SensitiveDetectorName, collectionName[0]);

    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    hce->AddHitsCollection(fHCID, fHC);

	if(fEdep) 
		fEdep->clear();
	else
		fEdep = new std::vector<double>;
	
	std::vector<int>().swap(*fCounter);
	std::vector<double>().swap(*fHitEk);
	std::vector<double>().swap(*fHitTime);
	std::vector<double>().swap(*fHitX);
	std::vector<double>().swap(*fHitY);
	std::vector<double>().swap(*fHitZ);
	std::vector<double>().swap(*fHitPX);
	std::vector<double>().swap(*fHitPY);
	std::vector<double>().swap(*fHitPZ);
	std::vector<int>().swap(*fHitID);
}

G4bool CryPositionSD::ProcessHits(G4Step *theStep, G4TouchableHistory *)
{
	// Sensitive only for Primary track
	if(theStep->GetTrack()->GetParentID() != 0)
		return false;
	
	G4Track *theTrack = theStep->GetTrack();	
    G4double edep = theStep->GetTotalEnergyDeposit();
    if(edep <= 0) return false;
    
	G4StepPoint* theParticle = theStep->GetPostStepPoint();
	if(!fNphysvol)
		CalculateNoPhysvols(theStep->GetPreStepPoint());
    
    CryHit* newHit = new CryHit();
    newHit->SetEdep(edep);
	newHit->SetDetectorID(CalculateCopyNo(theStep->GetPreStepPoint()));

    fHC->insert(newHit);
	
	fHitEk->push_back(theParticle->GetKineticEnergy());
	fHitTime->push_back(theParticle->GetGlobalTime());

	fHitID->push_back(theTrack->GetTrackID());

	G4ThreeVector pos = theParticle->GetPosition();
	fHitX->push_back(pos.x());
	fHitY->push_back(pos.y());
	fHitZ->push_back(pos.z());

	G4ThreeVector dir = theParticle->GetMomentumDirection();
	fHitPX->push_back(dir.x());
	fHitPY->push_back(dir.y());
	fHitPZ->push_back(dir.z());

    return true;
}

void CryPositionSD::EndOfEvent(G4HCofThisEvent*){

   	for(int i = 0 ; i < fNvolume ; i++)
		{fEdep->push_back(0.);
		fCounter->push_back(0);}

    for (int i = 0; i < fHC->entries(); i++)
    {
		{(*fEdep)[(*fHC)[i]->GetDetectorID()]
			+= (*fHC)[i]->GetEdep();
		(*fCounter)[(*fHC)[i]->GetDetectorID()]++;}
    }

}

void CryPositionSD::CreateEntry(
	G4int ntupleID, G4RootAnalysisManager* rootData)
{
	fFirstColID =
		rootData->CreateNtupleDColumn(ntupleID, "int.Edep", *fEdep);
	rootData->CreateNtupleIColumn(ntupleID, "int.count", *fCounter);
	rootData->CreateNtupleDColumn(ntupleID, "int.E", *fHitEk);
	rootData->CreateNtupleDColumn(ntupleID, "int.t", *fHitTime);
	rootData->CreateNtupleDColumn(ntupleID, "int.x", *fHitX);
	rootData->CreateNtupleDColumn(ntupleID, "int.y", *fHitY);
	rootData->CreateNtupleDColumn(ntupleID, "int.z", *fHitZ);
	rootData->CreateNtupleDColumn(ntupleID, "int.px", *fHitPX);
	rootData->CreateNtupleDColumn(ntupleID, "int.py", *fHitPY);
	rootData->CreateNtupleDColumn(ntupleID, "int.pz", *fHitPZ);
	rootData->CreateNtupleIColumn(ntupleID, "int.trackID", *fHitID);
}

void CryPositionSD::FillEntry(
	G4int, G4RootAnalysisManager*)
{}