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
//  Muonium(2S) relaxation depending on the experienced electrical field
//  Id    : probability_Mu2SRelaxation.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-06
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "probability_Mu2SRelaxation.hh"
#include "TF1.h"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdlib.h>

/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 The Muonium Yield function as well as the parameters are taken from:
 M. Gonin, R. Kallenbach, P. Bochsler: "Charge exchange of hydrogen atoms
 in carbon foils at 0.4 - 120 keV", Rev.Sci.Instrum. 65 (3), March 1994
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

Mu2SProbability:: Mu2SProbability(){;}
Mu2SProbability::~Mu2SProbability(){;}

G4double Mu2SProbability::GetProbability(
       G4ThreeVector E,          // electrical field vector
       G4double delta_time)     // time step size
{
  // Parameter NAMES  for the decay rate
  double e_charge;
  double bohr_radius;
  double energy_diff;
  double gamma_2p;
  double hbar, hbar_ev;
  double e_field;

  double decay_rate_2s;
  double lifetime;
  double probability;
  

  // VALUES for the decay rate
  e_charge = 1.60217662*pow(10,-19); //electrical charge, coulomb
  bohr_radius = 0.532*pow(10,-10); //bohr radius for muonium in m
  hbar = 1.054571817*pow(10,-34); //plank constant, J/s
  energy_diff = 6.9377*pow(10,-25);// energy difference 2S-2P (hbar*w) in J
  gamma_2p = 1/(1.6*pow(10,-9)); // in s
  e_field = sqrt(pow(E.x(),2)+pow(E.y(),2)+pow(E.z(),2));//field strengh in V/m
  decay_rate_2s = 9*pow(e_charge,2)*pow(bohr_radius,2)*pow(e_field,2)*gamma_2p*(1/(pow(energy_diff,2)+pow(hbar_ev*gamma_2p,2)/4));
  lifetime = pow(decay_rate_2s,-1);
  
  probability = 1-TMath::Exp(-delta_time*decay_rate_2s); 
  /*if(probability > 0.05) {
      std::cout << "Probability: " << probability << std::endl;
      std::cout << "E Field: " << e_field/100 << " V/cm" << std::endl;
      std::cout << "Time difference: " << delta_time << std::endl;}*/
  return probability;

}
