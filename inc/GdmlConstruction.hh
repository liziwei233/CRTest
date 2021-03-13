/*
*   File : CRTest/inc/GdmlConstruciton.hh
*
*   Brief: Define of Detector System by GDML
*
*/

#ifndef CRTest_GdmlConstruction_h
#define CRTest_GdmlConstruction_h

#include "SysConstruction.hh"

#include "G4GDMLParser.hh"
#include "G4GDMLAuxStructType.hh"

#include "G4VPhysicalVolume.hh"

class G4UserLimits;

class GdmlConstruction : public SysConstruction
{
public:
    GdmlConstruction(G4GDMLParser* gdml);
	GdmlConstruction(G4String gdmlFileName);
    virtual ~GdmlConstruction();

    virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();
	
	void DumpStructure();

private:
	void Init();
	void DumpVolume(G4VPhysicalVolume* physvol, 
		G4String prefix, G4bool expanded=true) const;

	void ReadAuxiliary();
	void PrintAuxiliary(const G4GDMLAuxListType*,G4String);
	void ReadProperty(const G4GDMLAuxListType*,G4String);
	G4bool ReadBorderProperty(const G4GDMLAuxListType*,G4String);
	G4bool ReadSkinProperty(const G4GDMLAuxListType*,G4String);

private:
    G4VPhysicalVolume* fWorldPV;
	
	G4LogicalVolume* fTrigger;
	G4LogicalVolume* fTracker;
	G4LogicalVolume* fFTOF_medium;
	G4LogicalVolume* fR10754_window;
	G4LogicalVolume* fR10754_Sioil;
	G4LogicalVolume* fR10754_anode;
	G4LogicalVolume* fR10754_sensor;
	G4LogicalVolume* fFTOF_detector;
	G4LogicalVolume* fT0_medium;
	G4LogicalVolume* fT0_Alwrap;
	G4LogicalVolume* fR3809_anode;
	G4LogicalVolume* fR3809_window;
	G4LogicalVolume* fT0_Sioil;
	G4LogicalVolume* fT0_lightguide;
	G4LogicalVolume* fR3809_sensor;
	G4LogicalVolume* fT0_detector;
	
	G4GDMLParser* fGdml;
	G4UserLimits* fStepLimit;  
};

#endif // CRTest_GdmlConstruction_h