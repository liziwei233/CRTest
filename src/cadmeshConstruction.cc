/* ************************************************
 * GEANT4 VCGLIB/CAD INTERFACE - basic example
 *
 * File:      DetectorConstruction.cc
 *
 * Author:    Christopher M Poole,
 * Email:     mail@christopherpoole.net
 *
 * Date:      20th March, 2011
 **************************************************/

// USER //
#include "cadmeshConstruction.hh"

// CADMESH //
#include "CADMesh.hh"

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include <string>

cadmeshConstruction::cadmeshConstruction(const std::string & name)
{
    std::cout<<"read cadmesh file :"<<name<<std::endl;
    filename = name;


    G4NistManager * nist_manager = G4NistManager::Instance();
    G4Material * air = nist_manager->FindOrBuildMaterial("G4_AIR");
    G4Material * water = nist_manager->FindOrBuildMaterial("G4_WATER");
    // CAD model rotation. 
    //G4RotationMatrix * rot = new G4RotationMatrix();
    //rot->rotateZ(90*deg);
    
    // Load CAD file as tessellated solid //
    offset = G4ThreeVector(0, 0, 0);
    
    // Note that offset is applied to the points in mesh directly before placement.
    CADMesh * mesh = new CADMesh((char*) filename.c_str());
    mesh->SetScale(mm);
    mesh->SetOffset(offset);
    mesh->SetReverse(false);
    
    std::cout<<"cadmesh has been built !"<<std::endl;

    cad_solid = mesh->TessellatedMesh();
    //cad_logical = new G4LogicalVolume(cad_solid, water, "cad_logical", 0, 0, 0);
    //cad_physical = new G4PVPlacement(rot, G4ThreeVector(), cad_logical,
    //                                 "cad_physical", world_logical, false, 0);
    //cad_logical->SetVisAttributes(G4Color(0.5, 0.3, 1, 1));
    
    std::cout<<"cadmesh has been translated to TessellatedMesh !!"<<std::endl;

}

cadmeshConstruction::~cadmeshConstruction()
{
}



