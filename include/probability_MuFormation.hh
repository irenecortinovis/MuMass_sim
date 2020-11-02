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
//  Muonium formation on a SiO2 target and emission into vaccum
//  Id    : probability_Mu2SFormation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-06
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifndef SilicaMuProbability_h
#define SilicaMuProbability_h 1

#include "globals.hh"
#include "TMath.h"
#include "G4ThreeVector.hh"
/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 The efficiencies were extracted from DOI: 10.1103/PhysRevLett.108.143401
 * Unfortunately, the fit couldn't be reproduced (P(x,E) missing)
 * Right now, some values off the fit were extracted with MatLab & grabit
 * A spline through those datapoints describes the fit, returning the probability
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

class SilicaMuProbability
{

public:
  SilicaMuProbability(); // Class constructor
 ~SilicaMuProbability(); // Class destructor

  G4double GetProbability(G4double E, G4int Temperature);

private:   // Some internal variables

    G4double probability;
};

#endif
