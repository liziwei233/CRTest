/**
*   FILE : CRTest/src/OpRecorder.cc
*   Brief: Implementation of class OpRecorder
*/

#include "OpRecorder.hh"
#include "G4SystemOfUnits.hh"
#include<vector>

OpRecorder *OpRecorder::fgInstance = 0;

OpRecorder::OpRecorder()
    : VirtualRecorder(),
	nCerenkov(0),nScintTotal(0),nQuartz2Air(0),nQuartz2GlueL(0),
	nQuartz2GlueR(0),nWlsEmit(0),nGlue2PMTL(0),nGlue2PMTR(0),
	nDetectionL(0),nDetectionR(0),nCathodL(0),nCathodR(0),
      nBoundaryAbsorption(0), nBoundaryTransmission(0),
      nFresnelReflection(0),nTotalInternalReflection(0),nLambertianReflection(0),
      nLobeReflection(0),nSpikeReflection(0),nBackScattering(0),nBoundaryRefraction(0),
      nDebug(0),fID(NULL),fL(NULL),fBounce(NULL),fWaveL(NULL)
{
	boundaryName = "NULL";
    fCount = new std::vector<double>;
	fEk = new std::vector<double>;
	fTime = new std::vector<double>;
    fX = new std::vector<double>;
	fY = new std::vector<double>;
	fZ = new std::vector<double>;
	fPX = new std::vector<double>;
	fPY = new std::vector<double>;
	fPZ = new std::vector<double>;

	fID = new std::vector<int>;
	fL = new std::vector<double>;
	fWaveL = new std::vector<double>;
	fBounce = new std::vector<int>;
}

OpRecorder::~OpRecorder() {
    fCount->clear();delete fCount;
	fEk->clear();delete fEk;
	fTime->clear();delete fTime;
    fX->clear();delete fX;
	fY->clear();delete fY;
	fZ->clear();delete fZ;
	fPX->clear();delete fPX;
	fPY->clear();delete fPY;
	fPZ->clear();delete fPZ;

	fID->clear();delete fID;
	fL->clear();delete fL;
	fWaveL->clear();delete fWaveL;
	fBounce->clear();delete fBounce;
}

OpRecorder *OpRecorder::Instance()
{
    if (fgInstance == NULL)
        fgInstance = new OpRecorder();
    return fgInstance;
}

void OpRecorder::Reset()
{
    nCerenkov = 0;
    nScintTotal = 0;
    nQuartz2Air = 0;
    nQuartz2GlueL = 0;
    nQuartz2GlueR = 0;
	nWlsEmit = 0;
	nGlue2PMTL = 0;
    nGlue2PMTR = 0;
	nCathodL = 0;
	nCathodR = 0;
	nDetectionL = 0;
	nDetectionR = 0;
    
    nBoundaryAbsorption = 0;
    nBoundaryTransmission = 0;
    nBoundaryUndefined = 0;
    nBoundaryWARNNING = 0;

    nFresnelReflection = 0;
	nTotalInternalReflection = 0;
	nLambertianReflection = 0;
	nLobeReflection = 0;
	nSpikeReflection = 0;
	nBackScattering = 0;
	nBoundaryRefraction = 0;

    nDebug = 0;


    std::vector<double>().swap(*fCount);
    std::vector<double>().swap(*fEk);
    std::vector<double>().swap(*fTime);
    std::vector<double>().swap(*fX);
    std::vector<double>().swap(*fY);
    std::vector<double>().swap(*fZ);
    std::vector<double>().swap(*fPX);
    std::vector<double>().swap(*fPY);
    std::vector<double>().swap(*fPZ);

	std::vector<int>().swap(*fID);
	std::vector<double>().swap(*fL);
	std::vector<double>().swap(*fWaveL);
	std::vector<int>().swap(*fBounce);

}

void OpRecorder::SetBoundaryName(G4String name){
	boundaryName = name;
}

void OpRecorder::Print()
{
    G4cout << " | + Scintillation Total Count\t: " << nScintTotal << G4endl
           << " | + Cerenkov Total Count\t: " << nCerenkov << G4endl
           << " | + Quartz. to air Boundary\t: " << nQuartz2Air << G4endl
		   << " | + Quartz. to silicone Oil (LEFT)\t\t: " << nQuartz2GlueL << G4endl
           << " | + Quartz. to silicone Oil (Right)\t\t: " << nQuartz2GlueR << G4endl
		   << " | + Oil to Window (LEFT)\t\t: " << nGlue2PMTL << G4endl
           << " | + Oil to Window (Right)\t\t: " << nGlue2PMTR << G4endl
		   << " | + PMT (Left) Hits\t\t: " << nCathodL << G4endl
           << " | + PMT (Right) Hits\t\t: " <<  nCathodR << G4endl
           << " | + Detected by PMT (Left)\t\t: " << nDetectionL << G4endl
		   << " | + Detected by PMT (Right)\t\t: " << nDetectionR << G4endl
		   << " | + Boundary Details for " << boundaryName <<G4endl
		   << " | + + Boundary Transmission\t: " << nBoundaryTransmission << G4endl
		   << " | + + Boundary FresnelRefraction\t: " << nBoundaryRefraction << G4endl
           << " | + + Boundary FresnelReflection\t: " << nFresnelReflection << G4endl
           << " | + + Boundary TotalInternalReflection\t: " << nTotalInternalReflection << G4endl
           << " | + + Boundary LambertianReflection\t: " << nLambertianReflection << G4endl
           << " | + + Boundary LobeReflection\t: " << nLobeReflection << G4endl
           << " | + + Boundary SpikeReflection\t: " << nSpikeReflection << G4endl
           << " | + + Boundary BackScattering\t: " << nBackScattering << G4endl
           << " | + + Boundary Absorption\t: " << nBoundaryAbsorption << G4endl
           << " | + + Boundary Transmission\t: " << nBoundaryTransmission << G4endl
           << " | + + Boundary Undefined\t: " << nBoundaryUndefined << G4endl
           << " | + + Boundary WARNNING\t: " << nBoundaryWARNNING << G4endl
           << " | + X Count for Debug\t\t: " << nDebug << G4endl;
}

void OpRecorder::CreateEntry(G4int ntupleID, G4RootAnalysisManager* rootData)
{
	fFirstColID = 
		rootData->CreateNtupleIColumn(ntupleID, "op.crkov");
	rootData->CreateNtupleIColumn(ntupleID, "op.q2a");
	rootData->CreateNtupleIColumn(ntupleID, "op.q2gL");
	rootData->CreateNtupleIColumn(ntupleID, "op.q2gR");
	rootData->CreateNtupleIColumn(ntupleID, "op.wls");
	rootData->CreateNtupleIColumn(ntupleID, "op.g2pL");
    rootData->CreateNtupleIColumn(ntupleID, "op.g2pR");
	rootData->CreateNtupleIColumn(ntupleID, "op.det");

	rootData->CreateNtupleIColumn(ntupleID, "op.ID",*fID);
	rootData->CreateNtupleDColumn(ntupleID, "op.L",*fL);
	rootData->CreateNtupleDColumn(ntupleID, "op.WaveL",*fWaveL);
	rootData->CreateNtupleIColumn(ntupleID, "op.Bounce",*fBounce);

    rootData->CreateNtupleDColumn(ntupleID, "ph.Count", *fCount);
	rootData->CreateNtupleDColumn(ntupleID, "ph.E", *fEk);
	rootData->CreateNtupleDColumn(ntupleID, "ph.t", *fTime);
    rootData->CreateNtupleDColumn(ntupleID, "ph.x", *fX);
	rootData->CreateNtupleDColumn(ntupleID, "ph.y", *fY);
	rootData->CreateNtupleDColumn(ntupleID, "ph.z", *fZ);
	rootData->CreateNtupleDColumn(ntupleID, "ph.px", *fPX);
	rootData->CreateNtupleDColumn(ntupleID, "ph.py", *fPY);
	rootData->CreateNtupleDColumn(ntupleID, "ph.pz", *fPZ);
}

void OpRecorder::FillEntry(G4int ntupleID, G4RootAnalysisManager* rootData)
{
	rootData->FillNtupleIColumn(ntupleID, fFirstColID, nCerenkov);
	rootData->FillNtupleIColumn(ntupleID, fFirstColID+1, nQuartz2Air);
	rootData->FillNtupleIColumn(ntupleID, fFirstColID+2, nQuartz2GlueL);
	rootData->FillNtupleIColumn(ntupleID, fFirstColID+3, nQuartz2GlueR);
	rootData->FillNtupleIColumn(ntupleID, fFirstColID+4, nWlsEmit);
	rootData->FillNtupleIColumn(ntupleID, fFirstColID+5, nGlue2PMTL);
    rootData->FillNtupleIColumn(ntupleID, fFirstColID+6, nGlue2PMTR);
	rootData->FillNtupleIColumn(ntupleID, fFirstColID+7, nDetectionL+nDetectionR);
}

G4bool OpRecorder::Record(const G4Track* thePhoton)
{
    //more details refers to StepAction.C
    
    if(thePhoton->GetParentID() != 1)
		return false;
	fCount->push_back(1);
	fEk->push_back(thePhoton->GetKineticEnergy() );
	fTime->push_back(thePhoton->GetGlobalTime()  );

	G4ThreeVector pos = thePhoton->GetPosition();
	fX->push_back(pos.x() );
	fY->push_back(pos.y() );
	fZ->push_back(pos.z() );

	G4ThreeVector pph = thePhoton->GetMomentumDirection();
	fPX->push_back(pph.x());
	fPY->push_back(pph.y());
	fPZ->push_back(pph.z());
}
