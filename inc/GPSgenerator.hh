/*
*   File : CRTest/inc/Generator.hh
*
*   Brief: Generate primary event vertex and define messenger commands
*
*/

#ifndef GPSgenerator_h
#define GPSgenerator_h

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"

class GPSgenerator : public G4VUserPrimaryGeneratorAction {

public:
    GPSgenerator();
    virtual ~GPSgenerator();
public:
    virtual void GeneratePrimaries(G4Event* anEvent);
private:
    G4GeneralParticleSource*  fGeneralParticleSource;
};

#endif /*CRTest_GPSgenerator_h*/