
#include "GPSgenerator.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

GPSgenerator::GPSgenerator()
    : G4VUserPrimaryGeneratorAction(),
    fGeneralParticleSource(0)
{
    fGeneralParticleSource = new G4GeneralParticleSource();
}

GPSgenerator::~GPSgenerator()
{
    delete fGeneralParticleSource;
	G4cout << "[-] INFO - GPSgenerator deleted. " << G4endl;
}

void GPSgenerator::GeneratePrimaries(G4Event *anEvent)
{
/*
    fGeneralParticleSource->AddaSource(1);
    G4ParticleDefinition *particleDefinition =
    G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
    fGeneralParticleSource->SetParticleDefinition(particleDefinition);
    SetParticleTime();
    
    G4double theta,phi;
    G4int n=4;

    
    G4ThreeVector position = ();
    G4ThreeVector direction = ();

    G4SPSPosDistribution* posGenerator;
	G4SPSAngDistribution* angGenerator;
	G4SPSEneDistribution* eneGenerator;
	G4SPSRandomGenerator* biasRndm;   

    posGenerator->SetPosDisType(Point); 
    posGenerator->SetCentreCoords(position);
    posGenerator->SetPosRot1(0,0,1);
    posGenerator->SetPosRot2(1,0,0);

    angGenerator->SetAngDistType(cos);
    angGenerator->SetMinTheta(0. *deg);
    angGenerator->SetMaxTheta(28.57 *deg);

    eneGenerator->SetEnergyDisType()
*/
    fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
}
