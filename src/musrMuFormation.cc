/***************************************************************************
 *  musrSim - the program for the simulation of (mainly) muSR instruments. *
 *          More info on http://lmu.web.psi.ch/simulation/index.html .     *
 *          musrSim is based od Geant4 (http://geant4.web.cern.ch/geant4/) *
 *                                                                         *
 *  Copyright (C) 2009 by Paul Scherrer Institut, 5232 Villigen PSI,       *
 *                                                       Switzerland       *
 *                                                                         *
 *  This program is free software; you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation; either version 2 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program; if not, write to the Free Software            *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              *
 ***************************************************************************/

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muonium Formation according to yield.cc function (through GetYields method).
//  Id    : musrMuFormation.cc, v 1.4
//  Author: Taofiq PARAISO, T. Shiroka
//  Date  : 2007-12
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "musrMuFormation.hh"

using namespace std;

musrMuFormation::musrMuFormation(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType)
{
                   
                     random = new TRandom();
}

musrMuFormation::~musrMuFormation(){}

double musrMuFormation::landauMPV   = 0.1;
double musrMuFormation::landauSigma = 0.1;
double musrMuFormation::landauLoss = 0.0;
double musrMuFormation::MuFormation = 1.0;
double musrMuFormation::Mu2S_probability = 0.1;



G4VParticleChange* musrMuFormation::PostStepDoIt(const G4Track& trackData,
                                                 const G4Step&  aStep)
{ // Initialize ParticleChange  (by setting all its members equal to
  //                             the corresponding members in G4Track)
  fParticleChange.Initialize(trackData);

  G4Track theNewTrack;
  if(CheckCondition(aStep))
    {
      GetDatas(&aStep);
      G4Step theStep;
      PrepareSecondary( trackData);
      
      fParticleChange.AddSecondary(aSecondary);
      fParticleChange.ProposeTrackStatus(fStopAndKill) ;
    }  
  else
    {
      fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
    }
  return &fParticleChange;
}


G4bool musrMuFormation::CheckCondition(const G4Step& aStep)
{ // Decide when to call the MuFormation process - i.e. for muons going through the C foil.
  G4bool condition=false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name  
  //if(p_name == "mu+"&&aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName()=="log_CFoil") 
  std::string logVolName = aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName();
  if(p_name == "mu+" && ((logVolName=="log_coulombCFoil")||(logVolName=="log_CFoil")))
    {
      condition=true;
    }
  return condition;
}


G4double musrMuFormation::GetMeanFreePath(const G4Track&, 
                                                G4double,
                                                G4ForceCondition* condition)
{
  *condition = Forced;
   return DBL_MAX;
}


void musrMuFormation::GetDatas(const G4Step* aStep)
{    // Particle generation according to yield table
     particleTable=G4ParticleTable::GetParticleTable();
     rnd=G4UniformRand();
     
     G4double Eloss = 0;
     //cout << "Landau Loss is set to " << landauLoss << endl;
     //cout << "Landau MPV is set to " << landauMPV << endl;

     if(landauLoss == 1.0){
        do{
            Eloss = random->Landau(landauMPV,landauSigma);
        if (Eloss > 0.) break;
        } while (1); 
     }
     
     G4double E = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV-Eloss;
     if (E < 0. ) E = 0.; 
     //cout << E << " and " << Eloss << endl;
     Gonin.GetYields(E,105.658369*1000,yvector); // Energy [keV], muon mass [keV/c2], yield table
     G4String p_new = "Mu";
         
     // Positive muon
     if(p_name=="mu+")
       {
       if(MuFormation >= 1.0){
           if(MuFormation == 2.0){yvector[1] = 1.0;}
	 if(rnd>yvector[1]) 
	   {
	     particle = particleTable->FindParticle(p_name) ;
	   }
	 else 
	   {
            //cout << Mu2S_probability << endl;
            if(G4UniformRand() < Mu2S_probability){ //form Mu2S
                particle = particleTable->FindParticle("Mu2S");

            }
            else{ //form Mu
                particle = particleTable->FindParticle(p_new);
            }  
	   }
       }
       else{
	     particle = particleTable->FindParticle(p_name);
	   }

     // Set the new dynamic particle DP
        DP = new G4DynamicParticle(particle, 
	         aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(), E*CLHEP::keV);
	 
	 // IMPORTANT : COPY THOSE DATA TO GET THE SAME PARTICLE PROPERTIES!!!
	 // SHOULD BE KEPT WHEN BUILDING A PARTICLE CHANGE  
	 DP->SetProperTime(  aStep->GetTrack()->GetDynamicParticle()->GetProperTime());
	 DP->SetPolarization(aStep->GetTrack()->GetDynamicParticle()->GetPolarization().x(),
	                     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().y(),
			     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().z());
	 DP->SetPreAssignedDecayProperTime(aStep->GetTrack()->GetDynamicParticle()->GetPreAssignedDecayProperTime());
       }
}


void musrMuFormation::PrepareSecondary(const G4Track& track)
{
  if(p_name=="mu+")
    {
      aSecondary = new G4Track(DP,track.GetGlobalTime(),track.GetPosition());
    }
}
