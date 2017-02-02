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
// $Id: B2TrackerSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "B2TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TH1.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4UnitsTable.hh"

#include <fstream>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  std::vector<int> collectedPhoton;
  std::vector<int> totalPhoton;
  std::vector<int> trackIDs;
  TH1I* fractionalCollection = new TH1I("collection", "Optical Photon Fractional Collection", 100, 0, 20);
  TCanvas *c1 = new TCanvas("c1", "Total Energy Deposition in 25 Crystals", 200,10,700,900);
  TPad *pad1 = new TPad("pad1", "The pad with the 3D Hist",0.03,0.62,0.50,0.92,21);
 
  ofstream out;
B2TrackerSD::B2TrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
  collectedPhoton.push_back(0);
  totalPhoton.push_back(0);

  
  /*
  //make some tree things
    outtreeTrack = new TTree("track1GeV", "1GeV");
    outtreeTrack->Branch("eventID", &meventID, "eventID/I");
    outtreeTrack->Branch("particle", &mpdgID, "particle/I");
    outtreeTrack->Branch("trackID", &mtrackID, "trackID/I");
    outtreeTrack->Branch("parentID", &mparentID, "parentID/I");
    outtreeTrack->Branch("stepNum", &mstepNum, "stepNum/I");
    outtreeTrack->Branch("crystNum", &mcrystNum, "crystNum/I");
    outtreeTrack->Branch("xpos", &mxpos, "xpos/D");
    outtreeTrack->Branch("ypos", &mypos, "ypos/D");
    outtreeTrack->Branch("zpos", &mzpos, "zpos/D");
    outtreeTrack->Branch("kinE", &mkinE, "kinE/D");
    outtreeTrack->Branch("deltE", &mdeltE, "deltE/D");
    outtreeTrack->Branch("stepLen", &mstepLen, "stepLen/D");
    outtreeTrack->Branch("trackLen", &mtrackLen, "trackLen/D");    */  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//some root stuffs
	//std::vector<int> events;
	//TH1I* h1 = new TH1I("h1", "Step Length Per Event", 100, 0, 10000);
	//TCanvas *c1 = new TCanvas("c1", "Histogram in 1D", 200,10,700,900);
	//TPad *pad1 = new TPad("pad1","The pad with the 3D Hist",0.03,0.62,0.50,0.92,21);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::~B2TrackerSD() 
{ 
	
	collectedPhoton.clear();
	totalPhoton.clear();
	out.close();
	pad1->Draw();
	
	
	delete pad1;
	delete c1;
	delete fractionalCollection;
	
	//TFile* out = new TFile("/home/crystSim/TreeDataStorage/60-1GeV_tracks.root", "RECREATE");
	//outtreeTrack->Write();
	//out->Close();
	//delete outtreeTrack;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new B2TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // Collecting and counting optical photons
  G4String particle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4String next = aStep->GetTrack()->GetNextVolume()->GetName();
  G4int partID = aStep->GetTrack()->GetTrackID();
  bool found = false;
  for( int i = 0; i < trackIDs.size(); i++){
  	if( partID == trackIDs[i]){
  		found = true;
  		break;
  	}//end if 
  }//end for
  
  if( particle == "opticalphoton" && found == false){
  	trackIDs.push_back(partID);
  	int eventNum = totalPhoton.size() - 1;
  	totalPhoton[eventNum] += 1;
  	
  }//end if
  
  if(next == "Photodetector" || next == "World" && particle == "opticalphoton"){
  	int eventNum = totalPhoton.size() - 1;
  	collectedPhoton[eventNum] += 1;
  	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }//end if
  
  // Kill the photon if it starts wandering around too much
  G4double lenth = aStep->GetTrack()->GetTrackLength();
  if(lenth > 660){ //4x the crystal length
  	aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  	if(edep != 0){
		// writing it all to the TTree:
		//   but first lets define some values
		/*meventID = events.size();
		mpdgID = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
		mtrackID = aStep->GetTrack()->GetTrackID();
		mparentID = aStep->GetTrack()->GetParentID();
		mstepNum = aStep->GetTrack()->GetCurrentStepNumber();
		G4ThreeVector temp = aStep->GetTrack()->GetPosition();
		mcrystNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
		mxpos = temp.x();
		mypos = temp.y();
		mzpos = temp.z();
		mkinE = aStep->GetTrack()->GetKineticEnergy();
		mdeltE = aStep->GetTotalEnergyDeposit();
		mstepLen = aStep->GetStepLength();
		mtrackLen = aStep->GetTrack()->GetTrackLength();
		
		// now we can actually write it to the tree
		outtreeTrack->Fill(); */
	}
	
  if (edep==0.) return false;
  return true;
  B2TrackerHit* newHit = new B2TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert( newHit );
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent*)
{	
	int events = collectedPhoton.size();
	totalPhoton.push_back(0);
  	collectedPhoton.push_back(0);
  	trackIDs.clear();
  	double average = 0;
  	double total = 0;
  	for( int i = 0; i < collectedPhoton.size(); i++){
		total += collectedPhoton[i];
	}
	average = total/collectedPhoton.size();
  	double fracCollection = 100 * double(collectedPhoton[events - 1])/double(totalPhoton[events - 1]);
  	G4cout << "In event " << events << " there were a total of " << G4endl;
  	G4cout << totalPhoton[events - 1] << " optical photons generated." << G4endl;
  	G4cout << collectedPhoton[events - 1] << " of these were collected. The collection rate is" << G4endl;
  	G4cout << fracCollection << "% \n"; 
  	G4cout << average << " is the average number of collected photons, so far." << G4endl << G4endl;
  	
  	pad1->Draw();
  	fractionalCollection->Fill(collectedPhoton[events - 1]);
  	fractionalCollection->Draw();
  	
  	
  	 
  	out.open("/home/crystSim/photonStats001GEV.txt", std::ofstream::app);
  	out << totalPhoton[events - 1] << " " << collectedPhoton[events - 1] << std::endl;
  	out.close();
  if ( verboseLevel>1 ) { 
  	
     G4int nofHits = fHitsCollection->entries();
     //G4cout << G4endl
           // << "-------->Hits Collection: in this event they are " << nofHits 
           // << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
