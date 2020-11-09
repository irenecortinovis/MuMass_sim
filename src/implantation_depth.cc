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
//  Id    : implantation_depth.cc, v 1.0
//  Author: G. Janka
//  Date  : 2020-09
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "implantation_depth.hh"
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
 The parameters were extracted from TrimSP simulations of Thomas Prokscha
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

SilicaImplantDepth:: SilicaImplantDepth(){;}
SilicaImplantDepth::~SilicaImplantDepth(){;}

G4double SilicaImplantDepth::GetConstant(
       G4double E)        // implantation energy
{

  double_t graph_energy[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  double_t graph_constant[20] = {0.121814,0.102205,0.0927978,0.0870097,0.0823915,0.0788111,0.0764954,0.0734189,0.0711445,0.0692638,0.0678075,0.0654901,0.0638953,0.0623306,0.0610625,0.0598554,0.0585537,0.0575061,0.0555962,0.0547255};
  TGraph *gr = new TGraph(20, graph_energy, graph_constant);
  TSpline3 *s = new TSpline3("grs",gr);
  G4double constant = s->Eval(E/CLHEP::keV);

  return constant;
}

G4double SilicaImplantDepth::GetAlpha(
       G4double E)        // implantation energy
{

  double_t graph_energy[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  double_t graph_alpha[20] = {1.33882,0.912562,0.735104,0.676613,0.61516,0.587658,0.547345,0.547428,0.525994,0.544047,0.527182,0.518806,0.51848,0.498583,0.523297,0.500482,0.495397,0.488719,0.479994,0.477145};
  TGraph *gr = new TGraph(20, graph_energy, graph_alpha);
  TSpline3 *s = new TSpline3("grs",gr);
  G4double alpha = s->Eval(E/CLHEP::keV);

  return alpha;
}


G4double SilicaImplantDepth::GetMean(
       G4double E)        // implantation energy
{
  
  linear_constant = 76.527;
  linear_slope = 153.666;
  G4double mean = linear_constant + linear_slope * E/CLHEP::keV;
  return mean;
}

G4double SilicaImplantDepth::GetSigma(
       G4double E)        // implantation energy
{
  
  linear_constant = 139.91;
  linear_slope = 3.44;
  G4double sigma= linear_constant + linear_slope * E/CLHEP::keV;
  return sigma;
}

G4double SilicaImplantDepth::GetDepth(
       G4double E)        // implantation energy
{
  
  auto f1  = new TF1("f1","crystalball",0,4000);
  f1->FixParameter(0, GetConstant(E));
  f1->FixParameter(1, GetMean(E));
  f1->FixParameter(2, GetSigma(E));
  f1->FixParameter(3, GetAlpha(E));
  f1->FixParameter(4, 10000);
  f1->SetNpx(100);
    
  std::cout << GetConstant(E) << ", " << GetMean(E) << ", " << GetSigma(E) << ", " << GetAlpha(E) << std::endl;

  G4double implant_depth = f1->GetRandom(); //in nm
  delete f1;
  return implant_depth;
}
