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
	fname = name.erase(name.length()-2);
	buff = fname+"HC";
    collectionName.insert(buff);

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
	std::vector<int>().swap(*fHitCopyNo);
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
	G4VPhysicalVolume *thePostPV = theParticle->GetPhysicalVolume();
	G4cout<< "[-] thePostPV = "<<thePostPV->GetName()<<G4endl;
	buff = fname+"_PV";
	G4cout<< "[-] fname = "<<buff<<G4endl;
	if (thePostPV->GetName() != buff) return false;
	if(!fNphysvol)
		CalculateNoPhysvols(theStep->GetPreStepPoint());
    
    CryHit* newHit = new CryHit();
	G4int copyNo = CalculateCopyNo(theStep->GetPreStepPoint());
    newHit->SetEdep(edep);
	newHit->SetDetectorID(copyNo);

    fHC->insert(newHit);
	fHitCopyNo->push_back(copyNo);
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
	buff = fname+".Edep";
	fFirstColID =
	rootData->CreateNtupleDColumn(ntupleID, buff, *fEdep);

	buff = fname+".count";
	rootData->CreateNtupleIColumn(ntupleID, buff, *fCounter);
	
	buff = fname+".id";
	rootData->CreateNtupleIColumn(ntupleID, buff, *fHitCopyNo);

	buff = fname+".E";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitEk);

	buff = fname+".t";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitTime);

	buff = fname+".x";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitX);

	buff = fname+".y";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitY);

	buff = fname+".z";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitZ);

	buff = fname+".px";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitPX);

	buff = fname+".py";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitPY);

	buff = fname+".pz";
	rootData->CreateNtupleDColumn(ntupleID, buff, *fHitPZ);

	buff = fname+".trackID";
	rootData->CreateNtupleIColumn(ntupleID, buff, *fHitID);
}

void CryPositionSD::FillEntry(
	G4int, G4RootAnalysisManager*)
{}