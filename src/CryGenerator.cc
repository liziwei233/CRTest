/**
*   FILE : CRTest/src/CryDetector.cc
*   Brief: Implementation of class CryDetector
*   Reference : $CRYHOME/geant4/src/CRYPrimaryGeneratorAction.cc
*/

#include "CryGenerator.hh"

#include "RNGWrapper.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include<fstream>

CryGenerator::CryGenerator(G4String setupFile)
{
    fParticleGun = new G4ParticleGun();

    CRYFromFile(setupFile);

    fParticles = new std::vector<CRYParticle*>;

    fParticleTable = G4ParticleTable::GetParticleTable();
}

CryGenerator::~CryGenerator(){
    delete fParticleGun;
    delete fCryGenerator;
    delete fParticles;
}

void CryGenerator::CRYFromFile(G4String setupFile)
{
    // Read CRY setup file
    std::ifstream input;
    input.open(setupFile, std::ios::in);
    
    
    if(input.fail()){
        G4cout << "[X] ERROR - Failed to open CRY setup file : "
            << setupFile << " - by CRYGenerator."
            << G4endl;
        fInputState = -1; // Not setup.
        return;
    }// File ERROR
    
    G4String setupString("");
    char buffer[1024];
    while(!input.getline(buffer,1024).eof()){
        setupString += buffer;
        setupString += " ";
    }// Read settings
    
    // CRY Initialize
    UpdateCRY(setupString);
}

void CryGenerator::InputCRY(){
    fInputState = 1;
}

void CryGenerator::UpdateCRY(G4String input){
    // TODO : Check Memory Leak
    CRYSetup *setup = new CRYSetup(input,getenv("CRYDATAPATH"));
    fCryGenerator = new CRYGenerator(setup);
    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(
        CLHEP::HepRandom::getTheEngine(),
        &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    fInputState=0;
}

void CryGenerator::GeneratePrimaries(G4Event* anEvent)
{
    // Check CRYGenerator
    if(fInputState != 0){
        G4cout << "[X] ERROR - CryGenerator was not successfully initialized"
            << "- by CryGenerator."
            << G4endl;
        fParticleGun->GeneratePrimaryVertex(anEvent);
        return;
    }

    // Generate Secondary Particle(s)
    fParticles->clear();
    fCryGenerator->genEvent(fParticles);

    for(int i = 0 ; i < fParticles->size() ; i++){
        CRYParticle* particle = (*fParticles)[i];
        G4String pName = CRYUtils::partName(particle->id());
        fParticleGun->SetParticleDefinition(
            fParticleTable->FindParticle(particle->PDGid()));
        fParticleGun->SetParticleEnergy(particle->ke()*MeV);
        
        //G4ThreeVector position = GetWorldBoundary();

/*
//
//** put genaretor vetex upside y axis -> - y axis
//       
        position.setX(position.x()*(1-2*G4UniformRand()));
        position.setY(position.y());
        position.setZ(position.z()*(1-2*G4UniformRand()));
        fParticleGun->SetParticlePosition(position);
        fParticleGun->SetParticleMomentumDirection(
            G4ThreeVector(particle->u(),particle->w(),-particle->v()));
  */      
/*
//
//** put generator vetex on side plane   z axis -> - z axis
 
        position.setX(position.x()*(1-2*G4UniformRand()));
        position.setY(position.y()*(1-2*G4UniformRand()));
        position.setZ(position.z());
        fParticleGun->SetParticlePosition(position);
        fParticleGun->SetParticleMomentumDirection(
            G4ThreeVector(particle->u(),particle->v(),particle->w()));
*/

//
//** put generator vetex on side plane   x axis -> - x axis
        G4ThreeVector position(680.001,100.001,100.001);
        position.setX(position.x());
        //position.setY(0);
        //position.setZ(0);

        position.setY(position.y()*(1-2*G4UniformRand()));
        position.setZ(position.z()*(1-2*G4UniformRand()));
        fParticleGun->SetParticlePosition(position);
        fParticleGun->SetParticleMomentumDirection(
            G4ThreeVector(particle->w(),particle->v(),particle->u()));

        fParticleGun->GeneratePrimaryVertex(anEvent);

        delete particle;
    }// foreach in fParticles


}
