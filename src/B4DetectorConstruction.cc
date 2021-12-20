//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
    G4GlobalMagFieldMessenger *B4DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B4DetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{
    // Lead material defined using NIST Manager
    auto nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_Fe");
    G4Element *C = nistManager->FindOrBuildElement("C");
    G4Element *H = nistManager->FindOrBuildElement("H");

    G4Material *sci = new G4Material("scintillator", 0.909 * g / cm3, 2);
    sci->AddElement(C, 8);
    sci->AddElement(H, 8);

    G4Material *plastic = new G4Material("plastic", 0.909 * g / cm3, 2);
    plastic->AddElement(C, 8);
    plastic->AddElement(H, 8);

    // Vacuum
    new G4Material("Galactic", 1., 1.01 * g / mole, universe_mean_density,
                   kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B4DetectorConstruction::DefineVolumes()
{
    auto calorimeter_size_x = cell_width_on_x * cell_num_on_x_per_layer;
    auto calorimeter_size_y = cell_width_on_y * cell_num_on_y_per_layer;
    auto calorimeter_size_z = scintillator_layer_num * scintillator_thickness + (scintillator_layer_num - 1) * gap_thickness + 2 * (Fe_thickness + plastic_thickness);
    auto world_size_x = 1.2 * calorimeter_size_x;
    auto world_size_y = 1.2 * calorimeter_size_y;
    auto world_size_z = 1.2 * (calorimeter_size_z + 2 * 5 * m);
    auto first_position = -calorimeter_size_z / 2;

    // Get materials
    auto default_material = G4Material::GetMaterial("Galactic");
    auto Fe_material = G4Material::GetMaterial("G4_Fe");
    auto plastic_material = G4Material::GetMaterial("plastic");
    auto scintillator_material = G4Material::GetMaterial("scintillator");
    auto gap_material = G4Material::GetMaterial("G4_Pb");

    if (!default_material || !scintillator_material || !gap_material)
    {
        G4ExceptionDescription msg;
        msg << "Cannot retrieve materials already defined.";
        G4Exception("B4DetectorConstruction::DefineVolumes()",
                    "MyCode0001", FatalException, msg);
    }

    //
    // World
    //
    auto worldS = new G4Box("World",                                               // its name
                            world_size_x / 2, world_size_y / 2, world_size_z / 2); // its size

    auto worldLV = new G4LogicalVolume(
        worldS,           // its solid
        default_material, // its material
        "World");         // its name

    auto worldPV = new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        worldLV,         // its logical volume
        "World",         // its name
        0,               // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    //
    // Calorimeter
    //
    auto calorimeterS = new G4Box("Calorimeter",                                                           // its name
                                  calorimeter_size_x / 2, calorimeter_size_y / 2, calorimeter_size_z / 2); // its size

    auto calorLV = new G4LogicalVolume(
        calorimeterS,     // its solid
        default_material, // its material
        "Calorimeter");   // its name

    new G4PVPlacement(
        0,               // no rotation
        G4ThreeVector(), // at (0,0,0)
        calorLV,         // its logical volume
        "Calorimeter",   // its name
        worldLV,         // its mother  volume
        false,           // no boolean operation
        0,               // copy number
        fCheckOverlaps); // checking overlaps

    //
    // Fe in the front
    //
    auto front_Fe_S = new G4Box("front_Fe",                                                        // its name
                                calorimeter_size_x / 2, calorimeter_size_y / 2, Fe_thickness / 2); // its size

    auto front_Fe_LV = new G4LogicalVolume(
        front_Fe_S,  // its solid
        Fe_material, // its material
        "front_Fe"); // its name

    new G4PVPlacement(
        0,                                                        // no rotation
        G4ThreeVector(0., 0., first_position + Fe_thickness / 2), // its position
        front_Fe_LV,                                              // its logical volume
        "front_Fe",                                               // its name
        calorLV,                                                  // its mother  volume
        false,                                                    // no boolean operation
        0,                                                        // copy number
        fCheckOverlaps);                                          // checking overlaps

    //
    // plastic in the front
    //
    auto front_plastic_S = new G4Box("front_plastic",                                                        // its name
                                     calorimeter_size_x / 2, calorimeter_size_y / 2, plastic_thickness / 2); // its size

    auto front_plastic_LV = new G4LogicalVolume(
        front_plastic_S,  // its solid
        plastic_material, // its material
        "front_plastic"); // its name

    new G4PVPlacement(
        0,                                                                            // no rotation
        G4ThreeVector(0., 0., first_position + Fe_thickness + plastic_thickness / 2), // its position
        front_plastic_LV,                                                             // its logical volume
        "front_plastic",                                                              // its name
        calorLV,                                                                      // its mother  volume
        false,                                                                        // no boolean operation
        0,                                                                            // copy number
        fCheckOverlaps);                                                              // checking overlaps

    //
    // cell
    //
    G4int copy_number = 0;
    for (G4int iz = 0; iz < scintillator_layer_num; iz++)
    {
        for (G4int iy = 0; iy < cell_num_on_y_per_layer; iy++)
        {
            for (G4int ix = 0; ix < cell_num_on_x_per_layer; ix++)
            {
                G4double x_position = -calorimeter_size_x / 2 + ix * cell_width_on_x + cell_width_on_x / 2;
                G4double y_position = -calorimeter_size_y / 2 + iy * cell_width_on_y + cell_width_on_y / 2;
                G4double z_position = first_position + Fe_thickness + plastic_thickness + iz * (scintillator_thickness + gap_thickness) + scintillator_thickness / 2;

                auto cell_S = new G4Box("cell",                                                                // its name
                                        cell_width_on_x / 2, cell_width_on_y / 2, scintillator_thickness / 2); // its size

                auto cell_LV = new G4LogicalVolume(
                    cell_S,                // its solid
                    scintillator_material, // its material
                    "cell");               // its name

                new G4PVPlacement(
                    0,                                                 // no rotation
                    G4ThreeVector(x_position, y_position, z_position), // its position
                    cell_LV,                                           // its logical volume
                    "cell",                                            // its name
                    calorLV,                                           // its mother  volume
                    false,                                             // no boolean operation
                    copy_number,                                       // copy number
                    fCheckOverlaps);                                   // checking overlaps

                copy_number++;
            }
        }
    }

    //
    // Gap
    //
    for (G4int i = 0; i < scintillator_layer_num - 1; i++)
    {
        G4double z_position = first_position + Fe_thickness + plastic_thickness + i * (scintillator_thickness + gap_thickness) + scintillator_thickness + gap_thickness / 2;

        auto gap_S = new G4Box("gap",                                                              // its name
                               calorimeter_size_x / 2, calorimeter_size_y / 2, gap_thickness / 2); // its size

        auto gap_LV = new G4LogicalVolume(
            gap_S,        // its solid
            gap_material, // its material
            "gap");       // its name

        new G4PVPlacement(
            0,                               // no rotation
            G4ThreeVector(0, 0, z_position), // its position
            gap_LV,                          // its logical volume
            "gap",                           // its name
            calorLV,                         // its mother  volume
            false,                           // no boolean operation
            0,                               // copy number
            fCheckOverlaps);                 // checking overlaps
    }

    //
    // Fe in the back
    //
    auto back_Fe_S = new G4Box("back_Fe",                                                         // its name
                               calorimeter_size_x / 2, calorimeter_size_y / 2, Fe_thickness / 2); // its size

    auto back_Fe_LV = new G4LogicalVolume(
        back_Fe_S,   // its solid
        Fe_material, // its material
        "back_Fe");  // its name

    new G4PVPlacement(
        0,                                                                                                                                                                                              // no rotation
        G4ThreeVector(0., 0., first_position + (Fe_thickness + plastic_thickness) + scintillator_layer_num * scintillator_thickness + (scintillator_layer_num - 1) * gap_thickness + Fe_thickness / 2), // its position
        back_Fe_LV,                                                                                                                                                                                     // its logical volume
        "back_Fe",                                                                                                                                                                                      // its name
        calorLV,                                                                                                                                                                                        // its mother  volume
        false,                                                                                                                                                                                          // no boolean operation
        0,                                                                                                                                                                                              // copy number
        fCheckOverlaps);                                                                                                                                                                                // checking overlaps

    //
    // plastic in the back
    //
    auto back_plastic_S = new G4Box("back_plastic",                                                         // its name
                                    calorimeter_size_x / 2, calorimeter_size_y / 2, plastic_thickness / 2); // its size

    auto back_plastic_LV = new G4LogicalVolume(
        back_plastic_S,   // its solid
        plastic_material, // its material
        "back_plastic");  // its name

    new G4PVPlacement(
        0,                                                                                                                                                                                                                  // no rotation
        G4ThreeVector(0., 0., first_position + (Fe_thickness + plastic_thickness) + scintillator_layer_num * scintillator_thickness + (scintillator_layer_num - 1) * gap_thickness + Fe_thickness + plastic_thickness / 2), // its position
        back_plastic_LV,                                                                                                                                                                                                    // its logical volume
        "back_plastic",                                                                                                                                                                                                     // its name
        calorLV,                                                                                                                                                                                                            // its mother  volume
        false,                                                                                                                                                                                                              // no boolean operation
        0,                                                                                                                                                                                                                  // copy number
        fCheckOverlaps);                                                                                                                                                                                                    // checking overlaps

    //
    // print parameters
    //
    G4cout
        << G4endl
        << "------------------------------------------------------------" << G4endl
        << "---> The calorimeter is " << scintillator_layer_num << " layers of: [ "
        << scintillator_thickness / mm << "mm of " << scintillator_material->GetName()
        << " + "
        << gap_thickness / mm << "mm of " << gap_material->GetName() << " ] " << G4endl
        << "------------------------------------------------------------" << G4endl;

    //
    // Visualization attributes
    //
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    simpleBoxVisAtt->SetVisibility(true);
    calorLV->SetVisAttributes(simpleBoxVisAtt);

    //
    // Always return the physical World
    //
    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
