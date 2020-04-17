/**
*   FILE : CRTest/src/MuonRecorder.cc
*   Brief: Implementation of class MuonRecorder
*/

#include "MuonRecorder.hh"

#include "G4SystemOfUnits.hh"

MuonRecorder* MuonRecorder::fgInstance = NULL;

MuonRecorder::MuonRecorder()
	: VirtualRecorder()
{
	fCount = new std::vector<int>;
	fID = new std::vector<int>;
	fEk = new std::vector<double>;
	fTime = new std::vector<double>;
	fX = new std::vector<double>;
	fY = new std::vector<double>;
	fZ = new std::vector<double>;
	fPX = new std::vector<double>;
	fPY = new std::vector<double>;
	fPZ = new std::vector<double>;

	fDetID = new std::vector<int>;
}

MuonRecorder::~MuonRecorder(){
	fCount->clear();delete fCount;
	fID->clear();delete fID;
	fEk->clear();delete fEk;
	fTime->clear();delete fTime;
	fX->clear();delete fX;
	fY->clear();delete fY;
	fZ->clear();delete fZ;
	fPX->clear();delete fPX;
	fPY->clear();delete fPY;
	fPZ->clear();delete fPZ;

	fDetID->clear();delete fDetID;
}

MuonRecorder* MuonRecorder::Instance(){
	if(!fgInstance)
		fgInstance = new MuonRecorder;
	return fgInstance;
}

void MuonRecorder::Reset()
{
	std::vector<int>().swap(*fCount);
	std::vector<int>().swap(*fID);
	std::vector<double>().swap(*fEk);
	std::vector<double>().swap(*fTime);
	std::vector<double>().swap(*fX );
	std::vector<double>().swap(*fY );
	std::vector<double>().swap(*fZ );
	std::vector<double>().swap(*fPX);
	std::vector<double>().swap(*fPY);
	std::vector<double>().swap(*fPZ);

	std::vector<int>().swap(*fDetID);
}
void MuonRecorder::CreateEntry(
	G4int ntupleID, G4RootAnalysisManager* rootData)
{
	fFirstColID = 
		rootData->CreateNtupleIColumn(ntupleID, "mu.Count", *fCount);
	rootData->CreateNtupleIColumn(ntupleID, "mu.ID", *fID);
	rootData->CreateNtupleDColumn(ntupleID, "mu.E", *fEk);
	rootData->CreateNtupleDColumn(ntupleID, "mu.t", *fTime);
	rootData->CreateNtupleDColumn(ntupleID, "mu.x", *fX);
	rootData->CreateNtupleDColumn(ntupleID, "mu.y", *fY);
	rootData->CreateNtupleDColumn(ntupleID, "mu.z", *fZ);
	rootData->CreateNtupleDColumn(ntupleID, "mu.px", *fPX);
	rootData->CreateNtupleDColumn(ntupleID, "mu.py", *fPY);
	rootData->CreateNtupleDColumn(ntupleID, "mu.pz", *fPZ);

	rootData->CreateNtupleIColumn(ntupleID, "mu.DetID", *fDetID);
}

void MuonRecorder::FillEntry(G4int,G4RootAnalysisManager*)
{}

G4bool MuonRecorder::Record(const G4Track* theMuon){
	if(theMuon->GetParentID() != 0)
		return false;
	//if(fCount->size()>=1) return false;
	fCount->push_back(0);
	fID->push_back(theMuon->GetTrackID() );
	fEk->push_back(theMuon->GetKineticEnergy() / GeV );
	fTime->push_back(theMuon->GetGlobalTime() / ns );

	G4ThreeVector pos = theMuon->GetPosition();
	fX->push_back(pos.x() / mm );
	fY->push_back(pos.y() / mm );
	fZ->push_back(pos.z() / mm );

	G4ThreeVector pmu = theMuon->GetMomentumDirection();
	fPX->push_back(pmu.x());
	fPY->push_back(pmu.y());
	fPZ->push_back(pmu.z());
}