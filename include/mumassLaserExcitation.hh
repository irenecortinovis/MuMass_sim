//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Laser Excitation for Mu(1S) -> Mu(2S)
//  Id    : mumassLaserExcitation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-06
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                              
#ifndef   mumassLaserExcitation_h
#define   mumassLaserExcitation_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleTable.hh"
#include "probability_laserexcitation.hh"



/*! mumassLaserExcitation class defines the muonium(2s) formation process through laser excitation.
 */

class mumassLaserExcitation : public G4VDiscreteProcess
{
 public:
   
   mumassLaserExcitation(const G4String& name = "Mu_LaserExcitation", // process description
		   G4ProcessType aType = fElectromagnetic);

  ~mumassLaserExcitation();

  
  //! - Main method. Muonium(2S) formation process is executed at the END of a step. */
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
  
  Prob_LaserExc Mu_Prob_LaserExc; //reading from proability_laserexcitation the probability

  
  void GetDatas( const G4Step* aStep);
  // model parameters
  G4ParticleTable* particleTable; 
  G4ParticleDefinition* particle;
  G4ThreeVector mu_momentum_direction;
  G4double mu2s_formation_prob;
  G4double rnd;
  G4DynamicParticle *DP;
  
  //! The particle change object.
  G4VParticleChange fParticleChange; 
 
  void  PrepareSecondary(const G4Track&);
  G4Track* aSecondary;

  void InitializeSecondaries(const G4Track&);
};

#endif
