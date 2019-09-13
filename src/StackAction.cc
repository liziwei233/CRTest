/**
*   FILE : CRTest/src/StackAction.cc
*   Brief: Implementation of class StackAction
*/

#include "StackAction.hh"

#include "Analysis.hh"
#include "OpRecorder.hh"
#include "MuonRecorder.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"

#include "G4OpticalPhoton.hh"

StackAction::StackAction()
    : G4UserStackingAction()
{
}

StackAction::~StackAction()
{
	G4cout << "[-] INFO - StackAction deleted. " << G4endl;
}

G4ClassificationOfNewTrack
StackAction::ClassifyNewTrack(const G4Track *theTrack)
{


    //theTrack->GetUserInformation()->Print();
	//G4cout<<"particle TrackID is : "<<theTrack->GetTrackID()<<G4endl;
	
	// Record muon
	if(theTrack->GetParentID() == 0){
	//G4cout<<"Mark:particle process is : "<<theTrack->GetCreatorProcess()->GetProcessName() <<G4endl;
	G4cout<<"particle name is : "<<theTrack->GetParticleDefinition()->GetParticleName()<<G4endl;
	//G4cout<<"Mark:particle ParentID is : "<<theTrack->GetParentID()<<G4endl;
	//G4cout<<"Mark:particle TrackID is : "<<theTrack->GetTrackID()<<G4endl;
		MuonRecorder::Instance()->Record(theTrack);
		return fUrgent;
	}

    OpRecorder *Recorder = OpRecorder::Instance();
	

	//Count what process generated the optical photons
	if (theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
	{
		
		// particle is secondary
		if (theTrack->GetCreatorProcess()->GetProcessName() 
			== "Cerenkov")
		{
			//Analysis::Instance()->FillOpPhotonTrackForEvent(
			//	theTrack, OpPhotonType::Scintillation);
			Recorder->nCerenkov++;
			//Recorder->Record(theTrack);
		}
		else if (theTrack->GetCreatorProcess()->GetProcessName() 
			== "Scintillation")
		{
			//Analysis::Instance()->FillOpPhotonTrackForEvent(
			//	theTrack, OpPhotonType::Scintillation);
			Recorder->nScintTotal++;
		}
		else if (theTrack->GetCreatorProcess()->GetProcessName() 
			== "OpWLS")
		{
			//Analysis::Instance()->FillOpPhotonTrackForEvent(
			//	theTrack, OpPhotonType::OpWLS);
			Recorder->nWlsEmit++;
		}
	return fWaiting;

	}

    return fKill;
}

void StackAction::NewStage() {}

void StackAction::PrepareNewEvent() {}
