//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muonium Formation in a silica target
//  Id    : mumassMuFormation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-05
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                              
#ifndef   mumassMuFormation_h
#define   mumassMuFormation_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleTable.hh"
#include "probability_MuFormation.hh"
#include "implantation_depth.hh"


/*! mumassMuFormation class defines the muonium formation process in a silica target.

 * The process is executed at the END of a step, i.e. the muon is converted into 
 * Muonium AFTER flying through the target */

class mumassMuFormation : public G4VDiscreteProcess
{
 public:
   
   mumassMuFormation(const G4String& name = "MuFormation_silica", // process description
		   G4ProcessType aType = fElectromagnetic);

  ~mumassMuFormation();

  
  //! - Main method. Muonium formation process is executed at the END of a step. */
  G4VParticleChange* PostStepDoIt(
			     const G4Track&,
			     const G4Step&);
  
  G4double GetMeanFreePath(const G4Track& aTrack,
			   G4double previousStepSize,
			   G4ForceCondition* condition);


  //! Condition for process application (step Object).
  G4bool CheckCondition(const G4Step& aStep);
  
  //! Condition for process application (step Pointer).
  G4bool CheckCondition(const G4Step* aStep);
  
  
  G4String  p_name;
  G4bool condition;
  
  G4double mu_formation_prob;
  SilicaMuProbability SilicaProbability; 
  void GetDatas( const G4Step* aStep);
  // model parameters
  G4ParticleTable* particleTable; 
  G4ParticleDefinition* particle;
  G4ThreeVector mu_momentum_direction;
  G4double mu_kinetic_energy;
  static double mu_temperature;
  G4double diffusion_time;

  
  G4double rnd;
  G4DynamicParticle *DP;
  
  //! The particle change object.
  G4VParticleChange fParticleChange; 
 
  void  PrepareSecondary(const G4Track&);
  G4Track* aSecondary;

  void InitializeSecondaries(const G4Track&);
  
  G4double MuFormationEKin(G4double temperature);
  G4double diffusion_constant;
  SilicaImplantDepth SilicaImpDepth;

};



#endif
