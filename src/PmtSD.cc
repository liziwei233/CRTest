/**
*   FILE : CRTest/src/PmtSD.cc
*   Brief: Implementation of class PmtSD
*/

#include "PmtSD.hh"

#include "G4SDManager.hh"

#include "G4String.hh"
#include "G4ios.hh"

PmtSD::PmtSD(G4String &name)
	: VirtualSD(name), fHCID(0), fHC(NULL), fFirstColID(-1),
	fCounter(NULL)
{
	fname=name;
    //collectionName.insert("PmtHC");
    collectionName.insert(fname.replace(4,10,"HC"));
	fCounter = new std::vector<int>;

}

PmtSD::~PmtSD()
{
	delete fCounter;

}

void PmtSD::Initialize(G4HCofThisEvent *hce)
{
    G4cout << "[-] INFO - PmtSD Initialized."
           << " - by PmtSD" << G4endl;
    fHC = new PmtHC(SensitiveDetectorName, collectionName[0]);

    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    hce->AddHitsCollection(fHCID, fHC);

	std::vector<int>().swap(*fCounter);
	std::vector<int>().swap(*fHitCopyNo);
	std::vector<double>().swap(*fHitEk);
	std::vector<double>().swap(*fHitTime);
	std::vector<double>().swap(*fFlyTime);
	std::vector<double>().swap(*fHitX);
	std::vector<double>().swap(*fHitY);
	std::vector<double>().swap(*fHitZ);
	std::vector<double>().swap(*fHitPX);
	std::vector<double>().swap(*fHitPY);
	std::vector<double>().swap(*fHitPZ);
	std::vector<int>().swap(*fHitID);
	


}

G4bool PmtSD::ProcessHits(G4Step *theStep, G4TouchableHistory *)
{

	G4Track *theTrack = theStep->GetTrack();	
    G4double edep = theStep->GetTotalEnergyDeposit();
    if(edep <= 0) return false;
    
	G4StepPoint* theParticle = theStep->GetPostStepPoint();

	if(!fNphysvol)
		CalculateNoPhysvols(theParticle);

	PmtHit* newHit = new PmtHit;

	G4int copyNo = CalculateCopyNo(theParticle);
	newHit->SetPmtID(copyNo);

    fHC->insert(newHit);
	
	// TODO : #ifdef CRTest_SD_MORE
	// TODO : move to DumpHit and call in EndOfEvent
	fHitCopyNo->push_back(copyNo);
	fHitEk->push_back(theParticle->GetKineticEnergy());
	fHitTime->push_back(theParticle->GetGlobalTime());
	fFlyTime->push_back(theParticle->GetLocalTime());

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

void PmtSD::EndOfEvent(G4HCofThisEvent*){
	G4cout << "fNvolume= "<< fNvolume << G4endl;

	for(int i = 0 ; i < fNvolume ; i++)
		fCounter->push_back(0);
    
	for (int i = 0; i < fHC->entries(); i++)
    {
        (*fCounter)[(*fHC)[i]->GetPmtID()]++;
    }

}

void PmtSD::CreateEntry(
	G4int ntupleID, G4RootAnalysisManager* rootData)
{
	fFirstColID = 
	rootData->CreateNtupleIColumn(ntupleID, fname.replace(4,10,".count"), *fCounter);
	// TODO : #ifdef CRTest_SD_MORE
	rootData->CreateNtupleIColumn(ntupleID, fname.replace(4,10,".id"), *fHitCopyNo);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".E"), *fHitEk);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".t"), *fHitTime);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".flyt"), *fFlyTime);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".x"), *fHitX);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".y"), *fHitY);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".z"), *fHitZ);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".px"), *fHitPX);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".py"), *fHitPY);
	rootData->CreateNtupleDColumn(ntupleID, fname.replace(4,10,".pz"), *fHitPZ);
	rootData->CreateNtupleIColumn(ntupleID, fname.replace(4,10,".trackID"), *fHitID);
}

void PmtSD::FillEntry(
	G4int, G4RootAnalysisManager*)
{
}