/**
*   FILE : CRTest/src/StepAction.cc
*   Brief: Implementation of class StepAction
*/

#include "StepAction.hh"

#include "Analysis.hh"
#include "OpRecorder.hh"
#include "MuonRecorder.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4SystemOfUnits.hh"

StepAction::StepAction()
{
}

StepAction::~StepAction()
{
    G4cout << "[-] INFO - StepAction deleted. " << G4endl;
}

void StepAction::UserSteppingAction(const G4Step *aStep)
{

    OpRecorder *Recorder = OpRecorder::Instance();

    G4Track *theTrack = aStep->GetTrack();
    G4StepPoint *thePrePoint = aStep->GetPreStepPoint();
    G4StepPoint *thePostPoint = aStep->GetPostStepPoint();
    G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
    G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();

    /*
    Recorder->fID->push_back(theTrack->GetTrackID());
    Recorder->fL->push_back(aStep->GetStepLength()/cm);
    Recorder->fWaveL->push_back(1240/(theTrack->GetKineticEnergy()/eV));
    */
    const G4VProcess *theProcess = fpSteppingManager->GetfCurrentProcess();

    // for Muon (primary track)
    if (theTrack->GetParentID() == 0)
    {

        if (thePostPoint->GetStepStatus() == fGeomBoundary)
        {

            if (thePrePV->GetName() == "World_PV" && thePostPV->GetName() == "Detector_PV")
            {
                MuonRecorder::Instance()->Record(theTrack);
                MuonRecorder::Instance()->fDetID->push_back(thePostPV->GetCopyNo());
            }
        }
        return;
    }

    //  for Optical
    if (theTrack->GetParticleDefinition() !=
        G4OpticalPhoton::OpticalPhotonDefinition())
        return;

    //
    // Boundary Check
    //
    OpPhotonType type = OpPhotonType::Nothing;
    if (thePostPoint->GetStepStatus() == fGeomBoundary)
    {
        assert(theProcess->GetProcessName() == "OpBoundary");
        G4OpBoundaryProcess *boundary = (G4OpBoundaryProcess *)theProcess;
        G4OpBoundaryProcessStatus status = boundary->GetStatus();
        G4bool gotThrough =
            (status == Transmission || status == FresnelRefraction);
        if (thePrePV->GetName() == "medium_PV" &&
                 thePostPV->GetName() == "Detector_PV")
        {
            type = Quartz2Air;
            Recorder->nQuartz2Air++;
            Recorder->fBounce->push_back(theTrack->GetTrackID());

            Recorder->SetBoundaryName("Quartz2Air");
            BoundaryStats(boundary);
        }
        else if (gotThrough)
        {
            // OpPthoton got through boundary
            if (thePrePV->GetName() == "medium_PV" &&
                thePostPV->GetName() == "SO_PV")
            {
                type = Quartz2Glue;
                Recorder->nQuartz2Glue += 1;
                //Recorder->SetBoundaryName("Quartz2Glue");
                //BoundaryStats(boundary);
            }
            
            else if (thePrePV->GetName() == "SO_PV" &&
                     thePostPV->GetName() == "Window_PV")
            {
                type = Glue2Pmt;
                Recorder->nGlue2PMT += 1;

                //Recorder->SetBoundaryName("Glue2PmtR");
                //BoundaryStats(boundary);
            }
            
        }
        //else if (thePrePV->GetName() == "lightguide_left_PV" &&
        //         thePostPV->GetName() == "PMT_left_PV")
        else if (thePrePV->GetName() == "Window_PV" &&
                 thePostPV->GetName() == "PMT_PV")
        {
            // OpPhoton hit PMT photocathode
            type = Cathod;
            Recorder->nCathod += 1;
            if (status == Detection)
            {
                type = Detected;
                Recorder->nDetection += 1;

                //return;
            }
            //Recorder->SetBoundaryName("CathodL");
            //BoundaryStats(boundary);
        }
        
        // For Debug boundary details
        /*
        else if (thePrePV->GetName() == "medium_PV" &&
                 thePostPV->GetName() == "Detector_PV")
        {
            type = Quartz2Air;
            Recorder->nQuartz2Air++;
            Recorder->fBounce->push_back(theTrack->GetTrackID());

            Recorder->SetBoundaryName("Quartz2Air");
            BoundaryStats(boundary);
        }
        */
    }
    //Analysis::Instance()->FillOpPhotonTrackForEvent(theTrack, type);
}

G4bool StepAction::BoundaryStats(G4OpBoundaryProcess *boundary)
{
    OpRecorder *Recorder = OpRecorder::Instance();
    switch (boundary->GetStatus())
    {
    case FresnelRefraction:
        Recorder->nBoundaryRefraction++;
        break;
    case Transmission:
        Recorder->nBoundaryTransmission++;
        break;
    case Absorption:;
    case Detection:
        Recorder->nBoundaryAbsorption++;
        break;
    case FresnelReflection:
        Recorder->nFresnelReflection++;
        break;
    case TotalInternalReflection:
        Recorder->nTotalInternalReflection++;
        break;
    case LambertianReflection:
        Recorder->nLambertianReflection++;
        break;
    case LobeReflection:
        Recorder->nLobeReflection++;
        break;
    case SpikeReflection:
        Recorder->nSpikeReflection++;
        break;
    case BackScattering:
        Recorder->nBackScattering++;
        break;
    case Undefined:
        Recorder->nBoundaryUndefined++;
        break;
    case NotAtBoundary:;
    case SameMaterial:;
    case StepTooSmall:;
    case NoRINDEX:;
        Recorder->nBoundaryWARNNING++;
        break;
    default:
        return false;
    }
    return true;
}