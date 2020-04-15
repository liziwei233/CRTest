/*
*   File : CRTest/src/GdmlConstruction.cc
*   Brief: Implementation of class GdmlConstruction
*/

#include "GdmlConstruction.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4GlobalMagFieldMessenger.hh"

#include "Analysis.hh"
#include "CryPositionSD.hh"
#include "PmtSD.hh"

#include "G4GDMLParser.hh"
#include "G4GDMLAuxStructType.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RotationMatrix.hh"
#include "G4UserLimits.hh"

#include "G4LogicalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

//#include "cadmeshConstruction.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4SystemOfUnits.hh"
#include <string>

#include <vector>

GdmlConstruction::GdmlConstruction(G4GDMLParser *gdml)
	: SysConstruction(), fWorldPV(NULL), fGdml(gdml), fStepLimit(NULL),
	  fPmtL(NULL), fPmtR(NULL), fRadianer(NULL), fPDetector(NULL), fQuartz(NULL)
{
	Init();
}

GdmlConstruction::GdmlConstruction(G4String gdmlFileName)
	: SysConstruction(), fWorldPV(NULL), fGdml(NULL), fStepLimit(NULL),
	  fPmtL(NULL), fPmtR(NULL), fRadianer(NULL), fPDetector(NULL), fQuartz(NULL)
{
	fGdml = new G4GDMLParser;
	fGdml->Read(gdmlFileName, false);

	this->Init();
}

GdmlConstruction::~GdmlConstruction()
{
	G4cout << "[-] INFO - GdmlConstruction deleted. " << G4endl;
	delete fGdml;
	fGdml = NULL;
	delete fStepLimit;
}

void GdmlConstruction::Init()
{
	assert(fGdml != NULL);
	fWorldPV = fGdml->GetWorldVolume();
	if (!fWorldPV)
	{
		G4cout << "[#] ERROR - CAN NOT FOUND WORLD SETUP - "
			   << __FILE__ << " " << __func__ << G4endl;
		exit(1);
	}
	G4LogicalVolumeStore *lvStore = G4LogicalVolumeStore::GetInstance();
	fWorld = lvStore->GetVolume("World", false);
	if (!fWorld)
		fWorld = fWorldPV->GetLogicalVolume();
	//fTarget = lvStore->GetVolume("Target", false);
	fDetector = lvStore->GetVolume("Detector",false);
	fRadianer = lvStore->GetVolume("medium", false);
	fPDetector = lvStore->GetVolume("PMT", false);

	fPmtL = lvStore->GetVolume("PMT_left", false);
	fPmtR = lvStore->GetVolume("PMT_right", false);

	fStepLimit = new G4UserLimits();
	fStepLimit->SetUserMaxTime(200 * ns);
	fRadianer->SetUserLimits(fStepLimit);

	/*
  G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  flgPV = pvStore->GetVolume("lightguide_right_PV",false);
  fLightguide = lvStore->GetVolume("lightguide_right",false);
  cadmeshConstruction cad("../models/lightguidepos.STL");
  G4cout << "[?]lightguide.STL has been read succesfully! - ";
  fLightguide->SetSolid(cad.GetCADSolid());
  G4RotationMatrix* rotm  = new G4RotationMatrix();
  rotm->rotateZ(180*deg);
  rotm->rotateY(180*deg);
  G4ThreeVector pos1 = G4ThreeVector(-7.5*mm,7.5*mm,-160*mm);
  flgPV->SetRotation(rotm);
  flgPV->SetTranslation(pos1);

  flgPV = pvStore->GetVolume("lightguide_left_PV",false);
  fLightguide = lvStore->GetVolume("lightguide_left",false);
  //cadmeshConstruction cad("../models/lightguidepos.STL");
  //G4cout << "[?]lightguide.STL has been read succesfully! - ";
  fLightguide->SetSolid(cad.GetCADSolid());
  //G4RotationMatrix* rotm  = new G4RotationMatrix();
  //rotm->rotateZ(180*deg);
  //rotm->rotateY(180*deg);
  G4ThreeVector pos2 = G4ThreeVector(-7.5*mm,-7.5*mm,160*mm);
  //flgPV->SetRotation(rotm);
  flgPV->SetTranslation(pos2);
*/
	DumpStructure();

	ReadAuxiliary();
}

G4VPhysicalVolume *GdmlConstruction::Construct()
{
	fWorld->SetVisAttributes(G4VisAttributes::Invisible);
	fDetector->SetVisAttributes(G4VisAttributes::Invisible);

	G4VisAttributes *visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0,0.3)); //white
	 visAttributes = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.5)); //green
	 visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0, 1.0)); //megenta
	 
	 visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 1.0)); //blue
	 //visAttributes->SetVisibility(false);
	 
	fRadianer->SetVisAttributes(visAttributes);
	 visAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 1.0)); //red
	fPDetector->SetVisAttributes(visAttributes);
	/*
	G4bool checkOverlaps = true;

	// ------------ Air ------------
	G4Material *air =
		G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

	// ------------- Al ----------------
	G4Material *Al =
		G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	// ------------ SiO2 ---------------
	G4Element *Si = new G4Element("Silicon", "Si", 14., 28.1 * g / mole);
	G4Element *O = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);

	G4Material *SiO2 = new G4Material("quartz", 2.2 * g / cm3, 2);
	SiO2->AddElement(Si, 1);
	SiO2->AddElement(O, 2);
	G4Material *PhotonDet = new G4Material("PhotonDet", 2.2 * g / cm3, 1);
	PhotonDet->AddElement(Si, 1);
	G4Material *Cathode = new G4Material("Cathode", 2.2 * g / cm3, 1);
	Cathode->AddElement(Si, 1);

	// ----------------- quartz ------------------
	//

	//
	// -----------------optical photon energy-----------------
	//
	G4MaterialPropertiesTable *MPT_Quartz = new G4MaterialPropertiesTable();

	G4double ErefractiveIndex_quartz[] =
		{
			1.909 * eV, 1.939 * eV, 1.969 * eV, 2.001 * eV, 2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
			2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV, 2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV,
			2.757 * eV, 2.820 * eV, 2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV, 3.353 * eV, 3.446 * eV,
			3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV, 4.002 * eV, 4.136 * eV}; //
	G4int nEntries = sizeof(ErefractiveIndex_quartz) / sizeof(G4double);
	G4double refractiveIndex_quartz[] =
		{
			1.4565, 1.4568, 1.4571, 1.4574, 1.4577, 1.4580, 1.4584, 1.4587, 1.4591, 1.4595,
			1.4599, 1.4603, 1.4608, 1.4613, 1.4618, 1.4623, 1.4629, 1.4635, 1.4641, 1.4648,
			1.4656, 1.4663, 1.4672, 1.4681, 1.4691, 1.4701, 1.4713, 1.4725, 1.4738, 1.4753,
			1.4769, 1.4787, 1.4806, 1.4827, 1.4851, 1.4878};
	assert(sizeof(refractiveIndex_quartz) == sizeof(ErefractiveIndex_quartz));
	MPT_Quartz->AddProperty("RINDEX", ErefractiveIndex_quartz, refractiveIndex_quartz, nEntries)->SetSpline(true);

	G4double Eabsorption[] =
		{
			1.907 * eV, 1.921 * eV, 1.946 * eV, 1.978 * eV, 2.019 * eV, 2.071 * eV, 2.117 * eV, 2.162 * eV, 2.215 * eV, 2.261 * eV,
			2.311 * eV, 2.375 * eV, 2.442 * eV, 2.498 * eV, 2.566 * eV, 2.634 * eV, 2.706 * eV, 2.772 * eV, 2.852 * eV, 2.937 * eV,
			3.011 * eV, 3.088 * eV, 3.174 * eV, 3.255 * eV, 3.335 * eV, 3.435 * eV, 3.524 * eV, 3.630 * eV, 3.737 * eV, 3.863 * eV,
			3.991 * eV, 4.120 * eV, 4.216 * eV};
	nEntries = sizeof(Eabsorption) / sizeof(G4double);
	G4double absorption[] =
		{
			2349 * m, 2283 * m, 2163 * m, 2026 * m, 1866 * m, 1689 * m, 1543 * m, 1423 * m, 1290 * m, 1188 * m,
			1086 * m, 975 * m, 873 * m, 798 * m, 714 * m, 647 * m, 581 * m, 527 * m, 465 * m, 417 * m,
			377 * m, 341 * m, 306 * m, 275 * m, 253 * m, 222 * m, 199 * m, 182 * m, 160 * m, 137 * m,
			124 * m, 106 * m, 102 * m};
	assert(sizeof(absorption) == sizeof(Eabsorption));
	MPT_Quartz->AddProperty("ABSLENGTH", Eabsorption, absorption, nEntries)->SetSpline("true");
	SiO2->SetMaterialPropertiesTable(MPT_Quartz);

	G4Box *solidQuartz = new G4Box("Quartz", 5, 5, 20);
	fQuartz = new G4LogicalVolume(solidQuartz, SiO2, "Quartz");
	G4VPhysicalVolume *physiQuartz = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), fQuartz, "Quartz", fWorld, false, 0, checkOverlaps);

	G4Box *solidPhotonDet = new G4Box("PhotonDet", 5, 5, 1);
	fPDetector = new G4LogicalVolume(solidPhotonDet, Al, "PhotonDet");
	G4VPhysicalVolume *physiPhotonDet = new G4PVPlacement(0, G4ThreeVector(0, 0, 21), fPDetector, "PhotonDet", fWorld, false, 0, checkOverlaps);

	// Quartz-Air Surface
	G4OpticalSurface *opQuartzAirSurface = new G4OpticalSurface("QuartzAirSurface");
	opQuartzAirSurface->SetType(dielectric_dielectric);
	opQuartzAirSurface->SetFinish(ground);
	opQuartzAirSurface->SetModel(unified);
	opQuartzAirSurface->SetSigmaAlpha(0.1 * deg);
	//new G4LogicalBorderSurface("QuartzAirSurface", physiQuartz, fWorldPV, opQuartzAirSurface);

	// Mirror Surface
	G4OpticalSurface *opMirrorSurface = new G4OpticalSurface("MirrorSurface");
	opMirrorSurface->SetType(dielectric_metal);
	opMirrorSurface->SetFinish(polished);
	opMirrorSurface->SetModel(unified);

	new G4LogicalBorderSurface("SideCoatingQuartzSurface", physiQuartz, fWorldPV, opMirrorSurface);

	// Quartz-Cathode Surface
	G4OpticalSurface *opQuartzCathodeSurface = new G4OpticalSurface("QuartzCathodeSurface");
	opQuartzCathodeSurface->SetType(dielectric_metal);
	opQuartzCathodeSurface->SetFinish(polished);
	opQuartzCathodeSurface->SetModel(unified);
	new G4LogicalBorderSurface("CathodeSurface", physiQuartz, physiPhotonDet, opQuartzCathodeSurface);

	//
	// ************************ Generate & Add Material Properties Table attached to the optical surfaces ******************************
	//
	G4double ephoton[2] = {0.4 * eV, 6.02 * eV};

	// ------------------- Optical Quartz Air Surface ---------------------
	G4double EreflecQuartzAir[] =
		{
			1.904 * eV, 2.055 * eV, 2.264 * eV, 2.584 * eV, 3.034 * eV, 3.163 * eV, 3.237 * eV, 3.314 * eV, 3.395 * eV, 3.474 * eV,
			3.581 * eV, 3.663 * eV, 3.737 * eV, 3.821 * eV, 3.881 * eV, 3.943 * eV, 3.993 * eV, 4.044 * eV, 4.074 * eV, 4.127 * eV,
			4.150 * eV, 4.181 * eV};
	nEntries = sizeof(EreflecQuartzAir) / sizeof(G4double);
	G4double reflecQuartzAir[] =
		{
			0.982, 0.981, 0.980, 0.979, 0.978, 0.977, 0.976, 0.976, 0.974, 0.973,
			0.971, 0.969, 0.967, 0.965, 0.963, 0.961, 0.960, 0.958, 0.957, 0.955,
			0.954, 0.953};
	assert(sizeof(EreflecQuartzAir) == sizeof(reflecQuartzAir));

	G4double specularlobeQuartzAir[2] = {0.0, 0.0};
	G4double specularspikeQuartzAir[2] = {1.0, 1.0};
	G4double backscatterQuartzAir[2] = {0, 0};

	G4MaterialPropertiesTable *SMPT_QuartzAir = new G4MaterialPropertiesTable();

	SMPT_QuartzAir->AddProperty("REFLECTIVITY", EreflecQuartzAir, reflecQuartzAir, nEntries);
	SMPT_QuartzAir->AddProperty("SPECULARLOBECONSTANT", ephoton, specularlobeQuartzAir, 2);
	SMPT_QuartzAir->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularspikeQuartzAir, 2);
	SMPT_QuartzAir->AddProperty("BACKSCATTERCONSTANT", ephoton, backscatterQuartzAir, 2);

	opQuartzAirSurface->SetMaterialPropertiesTable(SMPT_QuartzAir);

	// ------------------- Optical Mirror Surface ---------------------
	G4double EreflecMirror[] =
		{
			1.906 * eV, 1.914 * eV, 1.922 * eV, 1.931 * eV, 1.941 * eV, 1.958 * eV, 1.977 * eV, 1.997 * eV, 2.024 * eV, 2.055 * eV,
			2.089 * eV, 2.124 * eV, 2.150 * eV, 2.176 * eV, 2.200 * eV, 2.217 * eV, 2.235 * eV, 2.256 * eV, 2.279 * eV, 2.300 * eV,
			2.339 * eV, 2.361 * eV, 2.381 * eV, 2.400 * eV, 2.420 * eV, 2.439 * eV, 2.450 * eV, 2.466 * eV, 2.477 * eV, 2.522 * eV,
			2.554 * eV, 2.643 * eV, 2.671 * eV, 2.701 * eV, 2.821 * eV, 2.902 * eV, 3.016 * eV, 3.171 * eV, 3.259 * eV, 3.343 * eV,
			3.399 * eV, 3.431 * eV, 3.506 * eV, 3.546 * eV, 3.586 * eV, 3.668 * eV, 3.705 * eV, 3.749 * eV, 3.846 * eV, 3.893 * eV,
			3.935 * eV, 3.991 * eV, 4.094 * eV, 4.163 * eV};
	nEntries = sizeof(EreflecMirror) / sizeof(G4double);
	G4double reflecMirror[] =
		{
			0.869, 0.872, 0.875, 0.878, 0.881, 0.883, 0.886, 0.890, 0.894, 0.898,
			0.903, 0.908, 0.912, 0.915, 0.918, 0.920, 0.921, 0.923, 0.924, 0.925,
			0.926, 0.927, 0.928, 0.929, 0.930, 0.931, 0.932, 0.933, 0.935, 0.936,
			0.937, 0.938, 0.939, 0.940, 0.939, 0.938, 0.937, 0.936, 0.935, 0.934,
			0.933, 0.932, 0.931, 0.930, 0.929, 0.928, 0.927, 0.926, 0.925, 0.924,
			0.923, 0.922, 0.921, 0.920};
	assert(sizeof(EreflecMirror) == sizeof(reflecMirror));
	G4MaterialPropertiesTable *SMPT_Mirror = new G4MaterialPropertiesTable();
	SMPT_Mirror->AddProperty("REFLECTIVITY", EreflecMirror, reflecMirror, nEntries);

	opMirrorSurface->SetMaterialPropertiesTable(SMPT_Mirror);

	// ------------------- Optical Quartz PhotonDet-Cathode Surface ---------------------
	G4double reflecQuartzPhotonDet[2] = {0, 0};
	G4double EeffiQuartzPhotonDet[] =
		{
			1.900 * eV, 1.928 * eV, 1.961 * eV, 1.990 * eV, 2.048 * eV, 2.118 * eV, 2.178 * eV, 2.231 * eV, 2.314 * eV, 2.386 * eV,
			2.436 * eV, 2.509 * eV, 2.565 * eV, 2.646 * eV, 2.724 * eV, 2.789 * eV, 2.867 * eV, 2.931 * eV, 3.027 * eV, 3.108 * eV,
			3.239 * eV, 3.333 * eV, 3.458 * eV, 3.592 * eV, 3.753 * eV, 3.896 * eV, 4.122 * eV};
	nEntries = sizeof(EeffiQuartzPhotonDet) / sizeof(G4double);
	G4double effiQuartzPhotonDet[] =
		{
			0.0386, 0.0414, 0.0456, 0.0490, 0.0572, 0.0660, 0.0743, 0.0837, 0.0966, 0.1088,
			0.1197, 0.1317, 0.1466, 0.1632, 0.1795, 0.2023, 0.2173, 0.2362, 0.2537, 0.2725,
			0.2928, 0.2928, 0.2894, 0.2793, 0.2632, 0.2571, 0.2452};
	assert(sizeof(EeffiQuartzPhotonDet) == sizeof(effiQuartzPhotonDet));
	G4MaterialPropertiesTable *SMPT_Quartz_Cathode = new G4MaterialPropertiesTable();
	SMPT_Quartz_Cathode->AddProperty("REFLECTIVITY", ephoton, reflecQuartzPhotonDet, 2);
	SMPT_Quartz_Cathode->AddProperty("EFFICIENCY", EeffiQuartzPhotonDet, effiQuartzPhotonDet, nEntries);
	opQuartzCathodeSurface->SetMaterialPropertiesTable(SMPT_Quartz_Cathode);
*/
	return fWorldPV;
}

void GdmlConstruction::ConstructSDandField()
{
	G4SDManager* SDManager = G4SDManager::GetSDMpointer();
	
	//PmtSD *pmtSD = new PmtSD(sdName);
	//fPDetector->SetSensitiveDetector(pmtSD);
	//auto hodoscope1 = new B5HodoscopeSD("/hodoscope1");

	//fPDetector->SetSensitiveDetector(hodoscope1);

	//Analysis::Instance()->RegisterSD(pmtSD);
	if(fPDetector){
		// Create, Set & Register PmtSD
		G4String sdName = "PmtSD";
		PmtSD* pmtSD = new PmtSD(sdName);
        SDManager->AddNewDetector(pmtSD);  
		fPDetector->SetSensitiveDetector(pmtSD);

		Analysis::Instance()->RegisterSD(pmtSD);
		G4cout << "[-] INFO - pmtRSD has been registed succesfully!" << G4endl;
	}

	
	if(fRadianer){
		G4String sdName = "CryPostionSD";
		CryPositionSD* crySD = new CryPositionSD(sdName);
        G4cout << "[-] INFO - crySD has been set succesfully!" << G4endl;
		//SetSensitiveDetector(fRadianer, crySD);

		//Analysis::Instance()->RegisterSD(crySD);
		G4cout << "[-] INFO - crySD has been registed succesfully!" << G4endl;
	}
	if(fPmtL){
		// Create, Set & Register PmtSD
		G4String sdName = "PmtLSD";
		PmtSD* pmtLSD = new PmtSD(sdName);
        SDManager->AddNewDetector(pmtLSD);  
		fPmtL->SetSensitiveDetector(pmtLSD);
		

		Analysis::Instance()->RegisterSD(pmtLSD);
		G4cout << "[-] INFO - pmtLSD has been registed succesfully!" << G4endl;
	}
	if(fPmtR){
		// Create, Set & Register PmtSD
		G4String sdName = "PmtRSD";
		PmtSD* pmtRSD = new PmtSD(sdName);
        SDManager->AddNewDetector(pmtRSD);  
		fPmtR->SetSensitiveDetector(pmtRSD);

		Analysis::Instance()->RegisterSD(pmtRSD);
		G4cout << "[-] INFO - pmtRSD has been registed succesfully!" << G4endl;
	}
	
	return;
}

void GdmlConstruction::ReadAuxiliary()
{
	// Volume Auxiliary
	const G4LogicalVolumeStore *LVstore = G4LogicalVolumeStore::GetInstance();
	G4cout << "[-] INFO - Auxiliary Info. for Logical Volumes" << G4endl;
	std::vector<G4LogicalVolume *>::const_iterator LViter;
	for (LViter = LVstore->begin(); LViter != LVstore->end(); LViter++)
	{
		G4GDMLAuxListType auxList = fGdml->GetVolumeAuxiliaryInformation(*LViter);
		if (auxList.size() == 0)
			continue;
		G4cout << " - Logical Volume : "
			   << (*LViter)->GetName() << G4endl;
		PrintAuxiliary(&auxList, " | ");
	}
	// Userinfo Auxiliary
	G4cout << "[-] INFO - Auxiliary Info. for Global/Userinfo " << G4endl;
	PrintAuxiliary(fGdml->GetAuxList(), " | ");

	return;
}

void GdmlConstruction::PrintAuxiliary(
	const G4GDMLAuxListType *auxList, G4String prefix)
{
	for (std::vector<G4GDMLAuxStructType>::const_iterator auxIter = auxList->begin();
		 auxIter != auxList->end(); auxIter++)
	{
		G4cout << prefix << auxIter->type << " : "
			   << auxIter->value << " " << auxIter->unit
			   << G4endl;
		if (!auxIter->auxList)
			continue;

		if (auxIter->type == "Property")
			ReadProperty(auxIter->auxList, prefix + " + ");
		else
			PrintAuxiliary(auxIter->auxList, prefix + " | ");
	} // print contend in <aux></aux>
	return;
}

void GdmlConstruction::ReadProperty(
	const G4GDMLAuxListType *auxList, G4String prefix)
{

	G4cout << prefix << " Property for OpticalSurface" << G4endl;

	for (std::vector<G4GDMLAuxStructType>::const_iterator auxIter = auxList->begin();
		 auxIter != auxList->end(); auxIter++)
	{
		if (auxIter->type == "Skin")
		{
			G4cout << prefix << "Skin Surface Property : "
				   << auxIter->value << G4endl;
			if (auxIter->auxList)
				ReadSkinProperty(auxIter->auxList, prefix + " + ");
		}
		else if (auxIter->type == "Border")
		{
			G4cout << prefix << "Border Surface Property : "
				   << auxIter->value << G4endl;
			if (auxIter->auxList)
				ReadBorderProperty(auxIter->auxList, prefix + " + ");
		}
		else
		{
			G4cout << "[#] - ERROR - WRONG AUXTYPE for Property - "
				   << " Support 'Skin' and 'Border' ONLY " << G4endl;
		}
	}
}

G4bool GdmlConstruction::ReadSkinProperty(
	const G4GDMLAuxListType *auxList, G4String prefix)
{

	G4String SurfaceName;
	G4LogicalVolume *lvptr = NULL;
	G4Material *matptr = NULL;

	G4cout << prefix;
	for (std::vector<G4GDMLAuxStructType>::const_iterator auxIter = auxList->begin();
		 auxIter != auxList->end(); auxIter++)
	{

		if (auxIter->type == "SurfaceName")
		{
			SurfaceName = auxIter->value;
			G4cout << "Surface-" << auxIter->value << " | ";
		}
		else if (auxIter->type == "LVname")
		{
			G4cout << "LogVol-" << auxIter->value << " | ";
			lvptr =
				G4LogicalVolumeStore::GetInstance()->GetVolume(auxIter->value);
			if (!lvptr)
			{
				G4cerr << "[#] ERROR - Logical Volume NOT FOUND" << G4endl;
				return false;
			}
		}
		else if (auxIter->type == "Material")
		{
			G4cout << "Material-" << auxIter->value << " | ";
			matptr =
				G4Material::GetMaterial(auxIter->value);
			if (!matptr || !matptr->GetMaterialPropertiesTable())
			{
				G4cerr << "[#] ERROR - G4Material NOT FOUND" << G4endl;
				return false;
			}
		}
		else
		{
			G4cout << "[#] - ERROR - WRONG AUXTYPE for Skin - "
				   << " NEEDED 'LVname' and 'Material' " << G4endl;
		}
	}
	G4cout << G4endl;

	G4LogicalSkinSurface *Surface = NULL;
	G4OpticalSurface *OpSurf = NULL;
	Surface = G4LogicalSkinSurface::GetSurface(lvptr);
	if (Surface)
		OpSurf =
			dynamic_cast<G4OpticalSurface *>(Surface->GetSurfaceProperty());
	if (OpSurf && !OpSurf->GetMaterialPropertiesTable())
	{
		OpSurf->SetMaterialPropertiesTable(
			matptr->GetMaterialPropertiesTable());
		assert(SurfaceName == OpSurf->GetName());
	}

	return true;
}

G4bool GdmlConstruction::ReadBorderProperty(
	const G4GDMLAuxListType *auxList, G4String prefix)
{

	G4String SurfaceName;
	G4VPhysicalVolume *thePrePV = NULL;
	G4VPhysicalVolume *thePostPV = NULL;
	G4Material *matptr = NULL;

	G4cout << prefix;
	for (std::vector<G4GDMLAuxStructType>::const_iterator auxIter = auxList->begin();
		 auxIter != auxList->end(); auxIter++)
	{

		if (auxIter->type == "SurfaceName")
		{
			SurfaceName = auxIter->value;
			G4cout << "Surface-" << SurfaceName << " | ";
		}
		else if (auxIter->type == "PVname")
		{
			G4cout << "PhysVol-" << auxIter->value << " | ";
			G4VPhysicalVolume *physvol =
				G4PhysicalVolumeStore::GetInstance()->GetVolume(auxIter->value);
			if (!physvol)
			{
				G4cerr << "[#] ERROR - Physical Volume NOT FOUND" << G4endl;
				return false;
			}
			else if (!thePrePV)
				thePrePV = physvol;
			else if (!thePostPV)
				thePostPV = physvol;
			else
			{
				G4cerr << "[#] ERROR - TOO many Physical Volume FOUND" << G4endl;
				return false;
			}
		}
		else if (auxIter->type == "Material")
		{
			G4cout << "Material-" << auxIter->value << " | ";
			matptr =
				G4Material::GetMaterial(auxIter->value);
			if (!matptr || !matptr->GetMaterialPropertiesTable())
			{
				G4cerr << "[#] ERROR - G4Material NOT FOUND" << G4endl;
				return false;
			}
		}
		else
		{
			G4cout << "[-] - ERROR - WRONG AUXTYPE for Skin - "
				   << " NEEDED 'PVname' and 'Material' " << G4endl;
		}
	}
	G4cout << G4endl;

	G4LogicalBorderSurface *Surface = NULL;
	G4OpticalSurface *OpSurf = NULL;
	Surface = G4LogicalBorderSurface::GetSurface(thePrePV, thePostPV);
	if (Surface)
		OpSurf =
			dynamic_cast<G4OpticalSurface *>(Surface->GetSurfaceProperty());
	if (OpSurf && !OpSurf->GetMaterialPropertiesTable())
	{
		OpSurf->SetMaterialPropertiesTable(
			matptr->GetMaterialPropertiesTable());
		assert(SurfaceName == OpSurf->GetName());
	}

	return true;
}

void GdmlConstruction::DumpStructure()
{
	G4cout << "[-] INFO - CRTest Geometry Structure : " << G4endl;

	DumpVolume(fWorldPV, " | ");
}

void GdmlConstruction::DumpVolume(
	G4VPhysicalVolume *physvol,
	G4String prefix, G4bool expanded) const
{
	G4ThreeVector pos = physvol->GetTranslation();

	G4cout << prefix << physvol->GetName()
		   << "[" << physvol->GetCopyNo() << "] : "
		   << "Position(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")"
		   << G4endl;

	if (!expanded)
		return;

	G4LogicalVolume *lvptr = physvol->GetLogicalVolume();
	G4cout << prefix << " + " << lvptr->GetName() << " : ";
	lvptr->GetSolid()->DumpInfo();

	G4String lastPVName = "";
	for (int i = 0; i < lvptr->GetNoDaughters(); i++)
	{
		G4VPhysicalVolume *thePV = lvptr->GetDaughter(i);
		if (thePV->GetName() != lastPVName)
		{
			expanded = true;
			lastPVName = thePV->GetName();
		}
		else
			expanded = false;
		DumpVolume(thePV, prefix + " | ", expanded);
	}
}
