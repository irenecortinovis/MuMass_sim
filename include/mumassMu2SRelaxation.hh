//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muonium Formation in a silica target
//  Id    : mumassMu2SRelaxation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-05
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                              
#ifndef   mumassMu2SRelaxation_h
#define   mumassMu2SRelaxation_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleTable.hh"
#include "probability_Mu2SRelaxation.hh"



/*! mumassMu2SRelaxation class defines the Mu(2S) decay process into a Mu(1S) + gamma.
*/

class mumassMu2SRelaxation : public G4VDiscreteProcess
{
 public:
   
   mumassMu2SRelaxation(const G4String& name = "Mu2SRelaxation", // process description
		   G4ProcessType aType = fElectromagnetic);

  ~mumassMu2SRelaxation();

  
  //! - Main method. Decay is exexcuted at the END of a step. */
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
  G4bool decayed;
  
  
  G4bool GetDatas( const G4Step* aStep);
  // model parameters
  G4ParticleTable* particleTable; 
  G4ParticleDefinition* particle;
  G4ThreeVector mu_momentum_direction;
  G4ThreeVector gamma_direction;

  Mu2SProbability Prob2SRelaxation; //reading from proability_Mu2SRelaxation.cc the probability
  G4double mu_kinetic_energy;
  G4double mu2s_decay_prob;
  G4double gamma_energy;

  G4double rnd;
  G4DynamicParticle *DP;
  G4DynamicParticle *DP2;
  
  //! The particle change object.
  G4VParticleChange fParticleChange; 
 
  void  PrepareSecondary(const G4Track&);
  G4Track* aSecondary;
  G4Track* aGamma;

  void InitializeSecondaries(const G4Track&);
};

#endif
