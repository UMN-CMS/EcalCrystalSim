//
// ********************************************************************
// * License and Disclaimer *
// * *
// * The Geant4 software is copyright of the Copyright Holders of *
// * the Geant4 Collaboration. It is provided under the terms and *
// * conditions of the Geant4 Software License, included in the file *
// * LICENSE and available at http://cern.ch/geant4/license . These *
// * include a list of copyright holders. *
// * *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work make any representation or warranty, express or implied, *
// * regarding this software system or assume any liability for its *
// * use. Please see the license in the file LICENSE and URL above *
// * for the full disclaimer and the limitation of liability. *
// * *
// * This code implementation is the result of the scientific and *
// * technical work of the GEANT4 collaboration. *
// * By using, copying, modifying or distributing the software (or *
// * any work based on the software) you agree to acknowledge its *
// * use in resulting scientific publications, and indicate your *
// * acceptance of all terms of the Geant4 Software license. *
// ********************************************************************
//
/// \file hadronic/Hadr01/src/G4EmUserPhysics.cc
/// \brief Implementation of the G4EmUserPhysics class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName: G4EmUserPhysics
//
// Author: V.Ivanchenko 11.07.2012
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "globals.hh"
#include "G4EmUserPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VProcess.hh"
#include "G4EmProcessOptions.hh"
#include "G4AntiProton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4Transportation.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4EmSaturation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmUserPhysics::G4EmUserPhysics(G4int ver)
  : G4VPhysicsConstructor("User Optical Options"), fVerbose(ver)
{
  G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmUserPhysics::~G4EmUserPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmUserPhysics::ConstructParticle()
{
  //G4AntiProton::AntiProton();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmUserPhysics::ConstructProcess()
{
  G4cout << "G4EmUserPhysics::ConstructProcess()\n";
  //G4Transportation * theTransportationProcess = new G4Transportation(1);
  G4Cerenkov * theCerenkovProcess = new G4Cerenkov("Cerenkov");
  G4Scintillation * theScintillationProcess = new G4Scintillation("Scintillation");
  G4OpAbsorption * theAbsorptionProcess = new G4OpAbsorption();
  G4OpRayleigh * theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpMieHG * theMieHGScatteringProcess = new G4OpMieHG();
  G4OpBoundaryProcess * theBoundaryProcess = new G4OpBoundaryProcess();

  //theCerenkovProcess->DumpPhysicsTable();
  //theScintillationProcess->DumpPhysicsTable();
  //theRayleighScatteringProcess->DumpPhysicsTable();
  
  theCerenkovProcess->SetMaxNumPhotonsPerStep(20);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);
  
  theScintillationProcess->SetScintillationYieldFactor(1.);
  theScintillationProcess->SetTrackSecondariesFirst(true);
  
  // Use Birks Correction in the Scintillation process
  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  theScintillationProcess->AddSaturation(emSaturation);

  //auto aParticleIterator=GetParticleIterator();
  aParticleIterator->reset();
  while( (*aParticleIterator)() )
  {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
     if (theCerenkovProcess->IsApplicable(*particle))
     {
       pmanager->AddProcess(theCerenkovProcess);
       pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
     }
    
    if (theScintillationProcess->IsApplicable(*particle))
    {
      pmanager->AddProcess(theScintillationProcess);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
    }
    
    if (particleName == "opticalphoton")
    {
      
      //pmanager->AddProcess(theTransportationProcess);
      //pmanager->SetProcessOrderingToFirst(theTransportationProcess,idxAlongStep);
      //pmanager->SetProcessOrderingToFirst(theTransportationProcess,idxPostStep);
      pmanager->AddDiscreteProcess(theAbsorptionProcess);
      pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(theMieHGScatteringProcess);
      pmanager->AddDiscreteProcess(theBoundaryProcess);
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
