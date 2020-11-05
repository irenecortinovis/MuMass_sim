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

#include "probability_MuFormation.hh"
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4RunManager.hh"
#include "musrErrorMessage.hh"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdlib.h>

/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 The efficiencies were extracted from DOI: 10.1103/PhysRevLett.108.143401
 * Unfortunately, the fit couldn't be reproduced (P(x,E) missing)
 * Right now, some values off the fit were extracted with MatLab & grabit
 * A spline through those datapoints describes the fit, returning the probability
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

SilicaMuProbability:: SilicaMuProbability(){;}
SilicaMuProbability::~SilicaMuProbability(){;}

G4double SilicaMuProbability::GetProbability(
       G4double E,         // implantation energy
       G4int Temperature)    // temperature of target
{

    probability = 0.0;
    
    if(Temperature == 100)
    {
        //extracted data points, energy defined in keV
        double_t energy_100k[16] = {0.404, 0.531,0.759,1.240,1.971,2.955,3.888,5.702,7.366,9.154,11.143,13.258,15.071,16.481,18.117,18.923};
        double_t efficiency_100k[16] = {0.609, 0.576, 0.540, 0.478, 0.423, 0.366, 0.321, 0.252, 0.196, 0.155, 0.117, 0.088, 0.068, 0.053, 0.044, 0.038};
        TGraph *gr = new TGraph(16, energy_100k, efficiency_100k);
        TSpline3 *s = new TSpline3("grs",gr);
        probability = s->Eval(E/CLHEP::keV);

    }
    else if(Temperature == 250)
    {
        //energy defined in keV
        //double_t energy_250k[14] = {0.480, 1.085, 1.767, 2.623, 4.084, 5.596, 7.032, 9.022, 10.836, 12.850, 14.890, 16.753, 18.163, 18.918};
        //double_t efficiency_250k[14] = {0.608, 0.566, 0.530, 0.496, 0.446, 0.399, 0.362, 0.318, 0.277, 0.237, 0.202, 0.176, 0.158, 0.149};
        //TGraph *gr = new TGraph(14, energy_250k, efficiency_250k);
        //TSpline3 *s = new TSpline3("grs",gr);
        //probability = s->Eval(E/CLHEP::keV) ; 
        probability = 0;
    }
    
    else
    { //temperature not defined
        musrErrorMessage::GetInstance()->musrError(INFO,
			        "For this temperature, no Mu-formation efficiency is known",true);
        G4RunManager* fRunManager = G4RunManager::GetRunManager();
		fRunManager->AbortRun();
    }
  
  return probability;

}
