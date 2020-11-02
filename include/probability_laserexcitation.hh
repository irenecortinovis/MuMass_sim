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
//  Muonium2S as a function of laser power.
//          The method GetYields is used by LaserExcitation.
//  Id    : probability_laserexcitation.cc, v 1.0
//  Author: G. janka
//  Date  : 2020/06
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifndef Prob_LaserExc_h
#define Prob_LaserExc_h 1

#include "globals.hh"
#include "G4ios.hh" 
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh" 
#include "G4Step.hh" 


class Prob_LaserExc
{

public:
  Prob_LaserExc(); // Class constructor
 ~Prob_LaserExc(); // Class destructor

  G4double GetTransitionProbability(double TOF);

  virtual
  G4double LaserIntensity(G4ThreeVector MuLaserIntersection, G4ThreeVector MuInLaserPropagation, G4double time, G4double MuVelocity, G4bool twophoton);	
  virtual		  
   void Rkqs(G4double y[], G4double dydx[], G4int n, G4double x, G4double htry, G4double eps, G4double yscal[], G4double *hdid, G4double *hnext);
  // void Rkqs(G4double *y, G4double *dydx, G4int n, G4double x, G4double htry, G4double eps,
  //G4double *yscal, G4double *hdid, G4double *hnext);
  virtual
    void Rkck(G4double y[], G4double G4dydx[], G4int n, G4double x, G4double h,G4double yout[], G4double yerr[]);  
  //void Rkck(G4double *y, G4double *dydx, G4int n, G4double x, G4double h, G4double *yout,
  //G4double *yerr);  
  virtual
  void cderiv(G4double time, G4double c[], G4double dc[]);

  
  //from macro
  static double betarabi, betaioni, betastark;
  static double gdecay, edecay;
  static double detuning, LaserTimeOffset, DistanceSourceLaserX;
  static double LaserFWHM, lambda, w0, Epr;
  
  static double LaserPos_X, LaserPos_Y, LaserPos_Z;
  
  G4ThreeVector Mu_momentum_direction, Mu_location;
  G4double Mu_speed, delta_time;
  
private:
  // Some internal variables

  G4int nvar;
  G4double eps = 10e-12;
  G4double dtdid, dtnext, dt;
  
  G4double Pt;
  G4ThreeVector MuLaserIntersection1, MuInLaserDirection;
  G4double MuVelocity;
  
  G4double c[5];
  G4double dc[5];
  G4double  ProbMu2S;
  G4double  ProbMu1S;

};

#endif
