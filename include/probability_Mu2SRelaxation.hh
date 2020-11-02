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
//  Muonium yields as a function of initial mu+ energies.
//          The method GetYields is used by MuFormation.
//  Id    : yields.cc, v 1.1
//  Author: Taofiq PARAISO, T. Shiroka
//  Date  : 2007-12
//  Notes : First implemented in Fortran by A. Hofer
//          C++ conversion by T.K. Paraiso 04-2005
//          Slight modifications by T. Shiroka
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifndef Mu2SProbability_h
#define Mu2SProbability_h 1

#include "globals.hh"
#include "TMath.h"
#include "G4ThreeVector.hh"

/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 The Muonium Yield function as well as the parameters are taken from:
 M. Gonin, R. Kallenbach, P. Bochsler: "Charge exchange of hydrogen atoms
 in carbon foils at 0.4 - 120 keV", Rev.Sci.Instrum. 65 (3), March 1994
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

class Mu2SProbability
{

public:
  Mu2SProbability(); // Class constructor
 ~Mu2SProbability(); // Class destructor

  G4double GetProbability(G4ThreeVector E, G4double delta_time);

private:   // Some internal variables
  double e_charge;
  double bohr_radius;
  double energy_diff;
  double gamma_2p;
  double hbar, hbar_ev;
  double e_field;

  double decay_rate_2s;
  double lifetime;
  double probability;
};

#endif
