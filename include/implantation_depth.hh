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
//  Implantation Depth of Muonium formation on a SiO2 target and emission into vaccum
//  Id    : implantation_depth.hh, v 1.0
//  Author: G. Janka
//  Date  : 2020-09
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifndef SilicaImplantDepth_h
#define SilicaImplantDepth_h 1

#include "globals.hh"
#include "TMath.h"
#include "G4ThreeVector.hh"


class SilicaImplantDepth
{

public:
  SilicaImplantDepth(); // Class constructor
 ~SilicaImplantDepth(); // Class destructor

  G4double GetConstant(G4double E);
  G4double GetMean(G4double E);
  G4double GetAlpha(G4double E);
  G4double GetSigma(G4double E);
  G4double GetDepth(G4double E);

private:   // Some internal variables

    G4double linear_constant;
    G4double linear_slope;
    
};

#endif
