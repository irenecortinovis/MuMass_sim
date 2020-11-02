//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muonium Formation in a silica target
//  Id    : mumassMuFormation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-05
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#include "mumassMuFormation.hh"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

mumassMuFormation::mumassMuFormation(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType){}

mumassMuFormation::~mumassMuFormation(){}

double mumassMuFormation::mu_temperature = 250;

G4VParticleChange* mumassMuFormation::PostStepDoIt(const G4Track& trackData,
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


G4bool mumassMuFormation::CheckCondition(const G4Step& aStep)
{ // Decide when to call the MuFormation process - i.e. for muons implanted into silica target
  G4bool condition=false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name  
  std::string logVolName = aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName();
  if(p_name == "mu+" && logVolName=="log_SilicaTarget")
    {
      //cout << "Condition true" << endl;
      condition=true;
    }
  return condition;
}


G4double mumassMuFormation::GetMeanFreePath(const G4Track&, 
                                                G4double,
                                                G4ForceCondition* condition)
{
  *condition = Forced;
   return DBL_MAX;
}


void mumassMuFormation::GetDatas(const G4Step* aStep)
{    // Particle generation according to yield table
     particleTable=G4ParticleTable::GetParticleTable();
     rnd=G4UniformRand();
     //G4double E = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;

     //define here probability of Mu formation -->  comes from probability_MuFormation.cc
     G4double implantation_energy = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
     
     mu_formation_prob = SilicaProbability.GetProbability(implantation_energy, mu_temperature); 
     //cout << "The efficiency is " << mu_formation_prob << " at an implant energy of " << implantation_energy/CLHEP::keV << endl;
      //mu_formation_prob = 0;

     G4String p_new = "Mu";
      
     // Positive muon
     if(p_name=="mu+")

       {
	 if(rnd<mu_formation_prob) 
	   {
          //define here momentum direction (cosTheta distribution)
          G4double cosTheta = 2.*G4UniformRand()-1.;
          G4double sinTheta = sqrt(1.-cosTheta*cosTheta);
          G4double phi     = 2*TMath::Pi() * G4UniformRand();
          if(cosTheta>0.){ cosTheta*=-1.;}
          mu_momentum_direction = G4ThreeVector(sinTheta*cos(phi),sinTheta*sin(phi),cosTheta);
          //mu_momentum_direction = G4ThreeVector(0,0,-1);

          //define here kinetic energy, according to Maxwell-Boltzmann, temperature comes from macro input 
          G4double mu_kinetic_energy = MuFormationEKin(mu_temperature);
          particle = particleTable->FindParticle(p_new) ;
          DP = new G4DynamicParticle(particle, mu_momentum_direction, mu_kinetic_energy);
          
          //define here diffusion time
          G4double implant_depth = SilicaImpDepth.GetDepth(implantation_energy)*CLHEP::angstrom; //in angstrom
          
          
          //diffusion constants and formula for diffusion time taken from:
          //DOI: 10.1103/PhysRevLett.108.143401, Muonium Emission into Vacuum from Mesoporous Thin Films at Cryogenic Temperatures
          if(mu_temperature == 100){diffusion_constant = 4.2*pow(10,-5)*pow(CLHEP::cm,2)/CLHEP::s;} //G4units in mm^2 / ns, in paper given als cm^2/s
          else if(mu_temperature == 250){diffusion_constant = 1.6*pow(10,-4)*pow(CLHEP::cm,2)/CLHEP::s;}
          else{diffusion_constant = 10000000;} //no diffusion defined, making time negligible here          
          diffusion_time = pow(implant_depth,2)/diffusion_constant;
          //diffusion_time = 0;
          
          //cout << "Mu was formed with energy of " << mu_kinetic_energy /CLHEP::eV *1000.0 << " meV" << endl;
          //cout << "The diffusion took " << diffusion_time << " ns" << endl;

	   }
	 else 
	   {
	     particle = particleTable->FindParticle(p_name);
         diffusion_time = 0;
             DP = new G4DynamicParticle(particle, 
	         aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),0); //aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
              DP->SetProperTime(  aStep->GetTrack()->GetDynamicParticle()->GetProperTime());
              

             	//cout << "No Mu formation due to probability of " << mu_formation_prob << endl;
	   }


     // Set the new dynamic particle DP

	 
	 // IMPORTANT : COPY THOSE DATA TO GET THE SAME PARTICLE PROPERTIES!!!
	 // SHOULD BE KEPT WHEN BUILDING A PARTICLE CHANGE  
    DP->SetProperTime(aStep->GetTrack()->GetDynamicParticle()->GetProperTime()+diffusion_time);
	 DP->SetPolarization(aStep->GetTrack()->GetDynamicParticle()->GetPolarization().x(),
	                     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().y(),
			     aStep->GetTrack()->GetDynamicParticle()->GetPolarization().z());
	 DP->SetPreAssignedDecayProperTime(aStep->GetTrack()->GetDynamicParticle()->GetPreAssignedDecayProperTime());
       }
}


void mumassMuFormation::PrepareSecondary(const G4Track& track)
{
  if(p_name=="mu+")
    {
      aSecondary = new G4Track(DP,track.GetGlobalTime()+diffusion_time,track.GetPosition());
    }
}

G4double mw_boltzmann(Double_t *x, Double_t *par) {
    
    //par[0] = m
    //par[1] = kb
    //par[2] = temperature
    return (4*TMath::Pi()*pow(par[0]/(2*TMath::Pi()*par[1]*par[2]),3/2)*pow(x[0],2)*TMath::Exp(-par[0]*pow(x[0],2)/(2*par[1]*par[2])));
}
    
G4double mumassMuFormation::MuFormationEKin(G4double temperature)
{
  //Calculate the energy distribution of the Mu assuming a Maxwell-Boltzmann profile for a given temperature

  G4double k_b = 1.380658e-23; //Boltzmann constant J/K
  G4double m_e = 9.1093897e-31; //mass of electron kg
  G4double m_u = 1.8835316e-28; //mass of muon kg
  
  G4double m_mu = m_e + m_u; //mass of muonium kg
  
  TF1 *maxwell_boltzmann = new TF1("maxwell_boltzmann", mw_boltzmann, 0, 20000, 3);
  maxwell_boltzmann->FixParameter(0, m_mu);
  maxwell_boltzmann->FixParameter(1, k_b);
  maxwell_boltzmann->FixParameter(2, temperature);
  maxwell_boltzmann->SetNpx(100);
  
  G4double rnd_velocity = maxwell_boltzmann->GetRandom(); // m/s
  G4double rnd_energy = 0.5*pow(rnd_velocity,2)*m_mu*(6.242*pow(10,18))*CLHEP::eV; //we only give it here a G4unit!
  delete maxwell_boltzmann;
  
  return rnd_energy;
}
