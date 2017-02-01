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
// $Id: B2EventAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file B2EventAction.cc
/// \brief Implementation of the B2EventAction class

#include "B2EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "TH1.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


B2EventAction::B2EventAction()
: G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2EventAction::~B2EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories

  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  G4int eventID = event->GetEventID();
	/*
	//Only do this once, in a run of 1000 events.
	if(eventID % 999 == 0 && eventID != 0){	;
	//Declare the histogram
		TH1I* h1 = new TH1I("h1", "Step Length Per Event", 1000, 1000, 1000000);
		TH1F* h2 = new TH1F("h2", "Step Lengths Per Event", 1000, 0, 100);
	//Declare a read-in variable
		int hits;
		double lengths;
	//Filling the histogram
  		std::ifstream in;
 		in.open("/home/crystSim/trackInfo.txt");
 		while(in >> hits){
			h1->Fill(hits);
		}//end while
 		in.close();
		
		
		in.open("/home/crystSim/trackLenInfo.txt");
		while(in >> lengths){
			h2->Fill(lengths);
		}//end while
		in.close();

	//Get some information from the histogram,
	//  then write it to file.
		std::ofstream out;
		out.open("/home/crystSim/meanStepsEnergy.txt",std::ofstream::app);
		out << h1->GetMean(1) << " " << h1->GetMeanError(1) << std::endl;
		out.close();
		delete h1;
		
		out.open("/home/crystSim/meanTrackLen.txt", std::ofstream::app);
		out << h2->GetMean(1) << " " << h2->GetMeanError(1) << std::endl;
		out.close();
		delete h2;
		
		out.open("/home/crystSim/EnergyData.txt", std::ofstream::app);
		out << -999 << " " << -999 << " " << -999 << " x " << -999 << " " << -999 << " " << -999 << " x" << std::endl;
		out.close(); 

	//If we want to run at multiple energies
	//  we must delete all of the information in 
	//  out read-from file "/home/crystSim/trackInfo.txt"
	//  so that when we write to it in the next event
	//  the file is clean.

		std::ofstream delInfo;
		delInfo.open("/home/crystSim/trackInfo.txt", std::ofstream::out | std::ofstream::trunc);
		delInfo.close(); 
		
		delInfo.open("/home/crystSim/trackLenInfo.txt", std::ofstream::out | std::ofstream::trunc);
		delInfo.close();

	}//endif
	
	std::ofstream out;
	out.open("/home/crystSim/EnergyData.txt", std::ofstream::app);
		out << -999 << " " << -999 << " " << -999 << " x " << -999 << " " << -999 << " " << -999 << " x" << std::endl;
		out.close(); 
  //if ( eventID < 100 || eventID % 100 == 0) {
    //G4cout << ">>> Event: " << eventID  << G4endl;
   // if ( trajectoryContainer ) {
     // G4cout << "    " << n_trajectories
            // << " trajectories stored in this event." << G4endl;
    //}
    //G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
   // G4cout << "    "  
         //  << hc->GetSize() << " hits stored in this event" << G4endl;
  //}

	//want to do some root things right here.*/
	
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
