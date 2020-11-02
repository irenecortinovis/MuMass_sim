//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muonium Formation
//  Id    : mumassMuFormation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-05
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "mumassMu2SRelaxation.hh"
#include "F04GlobalField.hh"

using namespace std;

mumassMu2SRelaxation::mumassMu2SRelaxation(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType){}

mumassMu2SRelaxation::~mumassMu2SRelaxation(){}


G4VParticleChange* mumassMu2SRelaxation::PostStepDoIt(const G4Track& trackData,
                                                 const G4Step&  aStep)
{ // Initialize ParticleChange  (by setting all its members equal to
  //                             the corresponding members in G4Track)
  fParticleChange.Initialize(trackData);

  G4Track theNewTrack;
  if(CheckCondition(aStep))
    {
      G4bool decayed = GetDatas(&aStep);
      G4Step theStep;
      PrepareSecondary( trackData);
      if(decayed){
         fParticleChange.AddSecondary(aSecondary);
         fParticleChange.AddSecondary(aGamma);
         fParticleChange.ProposeTrackStatus(fStopAndKill) ;
        //cout << "Mu(2S) relaxed to the groundstate and emitted a Ly-a photon " << endl;

      }
      else{fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;}	
    }  
  else
    {
      fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;

    }
  return &fParticleChange;
}


G4bool mumassMu2SRelaxation::CheckCondition(const G4Step& aStep)
{ // Decide when to call the Mu(2S) decay process
  G4bool condition=false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name  
  if(p_name == "Mu2S")
    {
      condition=true;
    }
  return condition;
}


G4double mumassMu2SRelaxation::GetMeanFreePath(const G4Track&, 
                                                G4double,
                                                G4ForceCondition* condition)
{
  *condition = Forced;
   return DBL_MAX;
}


G4bool mumassMu2SRelaxation::GetDatas(const G4Step* aStep)
{    // Particle generation according to yield table
     particleTable=G4ParticleTable::GetParticleTable();
     rnd=G4UniformRand();
     //G4double E = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;

     //define here probability of Mu(2S) decaying
     G4double point[4];
     point[0] = aStep->GetPostStepPoint()->GetPosition().x();
     point[1] = aStep->GetPostStepPoint()->GetPosition().y();
     point[2] = aStep->GetPostStepPoint()->GetPosition().z();
     point[3] = 0;
     G4double Bfi[6]={0,0,0,0,0,0};
     F04GlobalField::getObject()->GetFieldValue(point,Bfi);
     G4ThreeVector E_field = G4ThreeVector(Bfi[3]/(CLHEP::volt/CLHEP::m),Bfi[4]/(CLHEP::volt/CLHEP::m),Bfi[5]/(CLHEP::volt/CLHEP::m)); //only taking E-Field
     
     G4double delta_time = aStep->GetDeltaTime()/CLHEP::s; //in s
     mu2s_decay_prob = Prob2SRelaxation.GetProbability(E_field, delta_time); // E_field [V/m], yield table
     G4bool decayed = false;
     G4String p_new = "Mu";
     //Excited Mu
     if(p_name=="Mu2S")
       {
	 if(rnd<mu2s_decay_prob) 
	   {
           
            G4LorentzVector mu2s_4vector = aStep->GetTrack()->GetDynamicParticle()->Get4Momentum();
            
            
            //define here momentum direction
            G4double random_theta = G4UniformRand()*TMath::Pi();
            G4double random_phi = G4UniformRand()*TMath::Pi()*2;
            G4double r = CLHEP::h_Planck * CLHEP::c_light / (122*CLHEP::nm);           
            gamma_direction = G4ThreeVector(r*cos(random_theta)*sin(random_phi),r*cos(random_phi)*sin(random_theta),r*cos(random_theta));
            DP2 = new G4DynamicParticle(particleTable->FindParticle("gamma"), gamma_direction);
            DP2->SetProperTime(aStep->GetTrack()->GetDynamicParticle()->GetProperTime());
            
            //for now, recoil caused by photon emission is neglected
            mu_momentum_direction = aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection();
            mu_kinetic_energy = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
            // we neglect the energy aspect because we haven't increased Mu(2S) energy compared to Mu(1S) either
	         particle = particleTable->FindParticle(p_new) ;
             DP = new G4DynamicParticle(particle, mu_momentum_direction, mu_kinetic_energy);
             decayed = true;
	   }
	 else 
	   {
             particle = particleTable->FindParticle(p_name);
             DP = new G4DynamicParticle(particle, 
	         aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),
	         aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
             DP2 = new G4DynamicParticle(particle, //just a dummy to prevent Seg Faults, will not be used afterwards
	         aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),
	         aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
	   }


     // Set the new dynamic particle DP

	 
	 // IMPORTANT : COPY THOSE DATA TO GET THE SAME PARTICLE PROPERTIES!!!
	 // SHOULD BE KEPT WHEN BUILDING A PARTICLE CHANGE  
	 DP->SetProperTime(  aStep->GetTrack()->GetDynamicParticle()->GetProperTime());
	 DP->SetPolarization(aStep->GetTrack()->GetDynamicParticle()->GetPolarization().x(),
	                     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().y(),
			     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().z());
	 DP->SetPreAssignedDecayProperTime(aStep->GetTrack()->GetDynamicParticle()->GetPreAssignedDecayProperTime());
       }

  return decayed;
}


void mumassMu2SRelaxation::PrepareSecondary(const G4Track& track)
{
  if(p_name=="Mu2S")
    {
      aSecondary = new G4Track(DP,track.GetGlobalTime(),track.GetPosition());
      aGamma = new G4Track(DP2, track.GetGlobalTime(), track.GetPosition());
 
    }
}


