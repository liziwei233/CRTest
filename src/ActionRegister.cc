/**
*   FILE : CRTest/src/ActionRegister.cc
*   Brief: Implementation of class ActionRegister
*/

#include "ActionRegister.hh"

#include "Generator.hh"
#include "CryGenerator.hh"
#include "PduGenerator.hh"
#include "GPSgenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackAction.hh"
#include "TrackingAction.hh"
#include "StepAction.hh"

ActionRegister::ActionRegister()
    : G4VUserActionInitialization()
{}

ActionRegister::~ActionRegister()
{}

void ActionRegister::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

void ActionRegister::Build() const
{
    SetUserAction(new GPSgenerator);
    SetUserAction(new Generator);
    //SetUserAction(new CryGenerator("./mac/setup.file"));
    SetUserAction(new PduGenerator);
    //theRunManager->SetUserAction(new PrimaryGeneratorAction(""));
    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
    SetUserAction(new StackAction);
	SetUserAction(new TrackingAction);
	SetUserAction(new StepAction);
}