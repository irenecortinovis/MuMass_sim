//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Laser Excitation for Mu(1S) -> Mu(2S)
//  Id    : mumassLaserExcitation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-06
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "mumassLaserExcitation.hh"

using namespace std;

mumassLaserExcitation::mumassLaserExcitation (const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType){}

mumassLaserExcitation::~mumassLaserExcitation(){}

G4VParticleChange* mumassLaserExcitation::PostStepDoIt(const G4Track& trackData,
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


G4bool mumassLaserExcitation::CheckCondition(const G4Step& aStep)
{ // Decide when to call the Mu2S-formation process - i.e. for Mu in "laser volume"
  G4bool condition=false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name  
  std::string logVolName = aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName();
  if(p_name == "Mu" && logVolName=="log_saveLaser")
    {
      condition=true;
    }
  return condition;
}


G4double mumassLaserExcitation::GetMeanFreePath(const G4Track&, 
                                                G4double,
                                                G4ForceCondition* condition)
{
  *condition = Forced;
   return DBL_MAX;
}


void mumassLaserExcitation::GetDatas(const G4Step* aStep)
{    // Particle generation according to yield table
     particleTable=G4ParticleTable::GetParticleTable();
     rnd=G4UniformRand();
     
     Mu_Prob_LaserExc.Mu_speed = aStep->GetTrack()->GetVelocity();
     Mu_Prob_LaserExc.Mu_location = aStep->GetPreStepPoint()->GetPosition();
     Mu_Prob_LaserExc.Mu_momentum_direction = aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection();
     Mu_Prob_LaserExc.delta_time = aStep->GetDeltaTime();
     
     G4double delta_time = aStep->GetDeltaTime();
     mu2s_formation_prob = Mu_Prob_LaserExc.GetTransitionProbability(delta_time); // delta time
     //cout << "The probability is " << mu2s_formation_prob <<endl;
     //cout << "The random nr is " << std::setprecision(15) << rnd << endl;
     
     //define here momentum direction
     //mu_momentum_direction = G4ThreeVector(1,0,0);

     //define here kinetic energy, comes from macro input 
     //mu_kinetic_energy = 30*CLHEP::eV;

     G4String p_new = "Mu2S";
    
     if(p_name=="Mu")

       {
     //mu2s_formation_prob *= 10000;
     //std::cout << mu2s_formation_prob << std::endl;
	 if(rnd<(mu2s_formation_prob)) 
	   {
	     particle = particleTable->FindParticle(p_new) ;
             DP = new G4DynamicParticle(particle, aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),
             aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
             cout << "Mu2S was formed"<< endl;
             std::cout << mu2s_formation_prob << std::endl;
	   }
	 else 
	   {
	     particle = particleTable->FindParticle(p_name);
             DP = new G4DynamicParticle(particle, 
	         aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),
	         aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
             	//cout << "No Mu2S formation<< endl;
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
}


void mumassLaserExcitation::PrepareSecondary(const G4Track& track)
{
  if(p_name=="Mu")
    {
      aSecondary = new G4Track(DP,track.GetGlobalTime(),track.GetPosition());
    }
}
