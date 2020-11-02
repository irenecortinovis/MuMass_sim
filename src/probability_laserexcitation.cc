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

#include "probability_laserexcitation.hh"
#include "G4PhysicalConstants.hh"
#include "TMath.h"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 The Muonium Yield function as well as the parameters are taken from:
 M. Gonin, R. Kallenbach, P. Bochsler: "Charge exchange of hydrogen atoms
 in carbon foils at 0.4 - 120 keV", Rev.Sci.Instrum. 65 (3), March 1994
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

Prob_LaserExc:: Prob_LaserExc(){;}
Prob_LaserExc::~Prob_LaserExc(){;}


//default values

double Prob_LaserExc::betarabi  = 1.01457928 * 3.68111*pow(10,-5)*CLHEP::hertz*pow(CLHEP::watt/pow(CLHEP::meter,2),-1);
double Prob_LaserExc::betaioni = 1.01457928 * 1.20208*pow(10,-4)*CLHEP::hertz*pow(CLHEP::watt/pow(CLHEP::meter,2),-1);
double Prob_LaserExc::betastark = 1.01457928 * (1.39927*pow(10,-4)-(-2.67827*pow(10,-5)))*CLHEP::hertz*pow(CLHEP::watt/pow(CLHEP::meter,2),-1);

double Prob_LaserExc::gdecay = pow(2.1969811*CLHEP::microsecond,-1);
double Prob_LaserExc::edecay = pow(2.1969811*CLHEP::microsecond,-1);

double Prob_LaserExc::detuning = 0;

//Only interesting for pulsed:
double Prob_LaserExc::LaserTimeOffset = 0;
double Prob_LaserExc::DistanceSourceLaserX = 0; 
double Prob_LaserExc::LaserFWHM = 0;

double Prob_LaserExc::lambda = 243*CLHEP::nanometer;
double Prob_LaserExc::w0 = 0.3*CLHEP::millimeter;
double Prob_LaserExc::Epr = 30*CLHEP::watt;

double Prob_LaserExc::LaserPos_X = 0;
double Prob_LaserExc::LaserPos_Y = 0;
double Prob_LaserExc::LaserPos_Z = 0;
  

G4double Prob_LaserExc::GetTransitionProbability(double TOF)
// time of flight, how long was the particle in the laser volume

{
    
  c[0] = 1.; 
  c[1] = 0.; 
  c[2] = 0.;
  c[3] = 0.;
  c[4] = 0.; 
    
  G4double yscal[5] = {1.,1.,1.,1.,1.};
  G4double t=0;
  dt=1*CLHEP::picosecond;
  dtnext=1*CLHEP::picosecond;
  G4int nnint=int(TOF/dt);
  G4int jk=0;
  nvar = 5;
  dt=dtnext;

  while(jk<nnint && t<TOF){
 
    dt=dtnext;
    cderiv(t,c,dc);
    Rkqs(c, dc, nvar, t, dt, eps, yscal, &dtdid, &dtnext);
    dt=dtdid;t+=dt;
  
    jk++;
  }
  
  ProbMu2S=c[3];
  ProbMu1S=c[0];
  G4double piyld=c[4];
  
  return ProbMu2S;
}


void Prob_LaserExc::cderiv(G4double time, G4double c[], G4double dc[]){
  
{
    //betarabi, betaioni, betastark from macro
    G4double wrabi=4*TMath::Pi()*betarabi*LaserIntensity(Mu_location, Mu_momentum_direction, time, Mu_speed,true);
    G4double wioni=2*TMath::Pi()*betaioni*LaserIntensity(Mu_location, Mu_momentum_direction, time, Mu_speed,false);
    G4double wstark=2*TMath::Pi()*betastark*LaserIntensity(Mu_location, Mu_momentum_direction, time, Mu_speed,false);

    complex<double> i(0.,1.);
    double gg = c[0];
    complex<double> ge(c[1],c[2]);
    double ee = c[3];
    double dgg,dee;
    complex<double> dge;

    //gdecay, edecay, detuning from macro
    dgg = -wrabi*imag(ge)-gdecay*gg;
    dge = -i*(detuning-wstark)*ge+i*(wrabi/2.)*(gg-ee)-ge*(wioni+gdecay+edecay)/2.;
    dee = wrabi*imag(ge)-(wioni+edecay)*ee;

    G4double dnpi=wioni*ee;

    dc[0]=dgg;
    dc[1]=real(dge);
    dc[2]=imag(dge);
    dc[3]=dee;
    dc[4]=dnpi;
}
}

G4double Prob_LaserExc::LaserIntensity(G4ThreeVector MuLaserIntersection, G4ThreeVector MuInLaserPropagation, G4double time, G4double MuVelocity, bool twophoton)
{
    
    // Compute the Intensity of the laser at a given point in W/m^2
    // Check out propagation of ops in z? => change components of oPsPos according to that...
    // MuPos in mm => tranform to m...-> *1.e-3
    // MuPos is a vector relative to in Laser beam axis
    // at (0,0,0) to check with TransitionRate.C

   
    G4double MuPos[3];
    G4ThreeVector LaserPos = G4ThreeVector(LaserPos_X*CLHEP::mm, LaserPos_Y*CLHEP::mm, LaserPos_Z*CLHEP::mm);
    G4double r2=0.;
    
    for(int j=0;j<3;j++)
    {
        MuPos[j]=(MuLaserIntersection[j]-LaserPos[j])+MuInLaserPropagation[j]*time*MuVelocity;
    }
    r2=pow(MuPos[1],2)+pow(MuPos[2],2); // y^2 + z^2, since x is laser propagation
    //Rayleigh range [m], w0 and lambda from macro
    G4double zR= TMath::Pi()*pow(w0,2)/lambda;
    
    //Beam radius [m] calculated for the geometry in which the beam axis is in x (it would more correct to write wx...)
    //w0 from macro
    G4double wz= w0*sqrt((1+pow(MuPos[0],2)/pow(zR,2)));
    //Pulsed Laser
    /* ---------- INCLUDE LASER TIME PROFILE (GAUSSIAN) HERE --------- */
    /*
    //LaserFWHM from macro
    G4double sigma = LaserFWHM / (2*sqrt(2*log(2)));
    // time needed for the reflected intensity, DistanceSourceLaserX from macro
    G4double time_reflected = 2*DistanceSourceLaserX/CLHEP::c_light;
    
    // 
    //Epr & LaserTimeOffset from macro
    G4double Pt1=(Epr/(sqrt(2*pi)*sigma))*exp(-pow(time-LaserTimeOffset,2)/(2*pow(sigma,2)));
    // Reflected Pulse
    G4double Pt2=(Epr/(sqrt(2*pi)*sigma))*exp(-pow(time+time_reflected-LaserTimeOffset,2)/(2*pow(sigma,2)));
    
   // Intensity able to excite
    if (twophoton)
    {
        Pt = 2*sqrt(Pt1*Pt2);
    }
    else
    {
        Pt = Pt1 + Pt2;
    }
    */
    Pt = Epr;
    G4double ILaser=2*Pt/(TMath::Pi()*pow(wz,2))*exp(-2*r2/pow(wz,2));
    //cout << Pt << " and " << wz << " and " << r2 << " and " << MuPos[1] << " and " << MuPos[2] << endl;
    //cout << "Laser Power in Watt/cm^2: " << ILaser/(CLHEP::watt/CLHEP::cm2) << endl;
    return ILaser;     
     
}     


//////// Don't touch below here, it's just Runge Kutta for numerical analysis!

void Prob_LaserExc::Rkqs(G4double *y, G4double *dydx, G4int n, G4double x, G4double htry, G4double eps, G4double *yscal, G4double *hdid, G4double *hnext)
{
 G4int i;
 G4double errmax,h,xnew,yerr[5],ytemp[5];
 const G4double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89E-4 ;

 //Fifth-order Runge-Kutta step with monitoring of local truncation error 
 //to ensure accuracy and adjust stepsize. Input are the dependent variable vector y[1..n] 
 //and its derivative dydx[1..n] at the starting value of the independent variable x. 
 //Also input are the stepsize to be attempted htry, 
 //the required accuracy eps, and the 
 //vector yscal[1..n] against which the error is scaled. 
 
 //On output, y and x are replaced by their new values, 
 //hdid is the stepsize that was actually accomplished, 
 //and hnext is the estimated next stepsize. 
 //derivs is the user-supplied routine that computes the right-hand side derivatives.


 h=htry; errmax = 5.0 ;
 while(errmax > 1.0)
 {
//    cout<<"n "<<n<<" x "<<x<<" h "<<h<<endl;
//    for(int gg=0;gg<5;gg++){
//        cout<<"y     "<<gg<<" "<<y[gg]<<endl;
//        cout<<"dxdy  "<<gg<<" "<<dydx[gg]<<endl;
//        cout<<"ytemp "<<gg<<" "<<ytemp[gg]<<endl;
//        cout<<"yerr  "<<gg<<" "<<yerr[gg]<<endl;
//      }
 Rkck(y,dydx,n,x,h,ytemp,yerr);
 errmax=0.0;
  for (i=0;i<n;i++)
  {
  if(errmax<fabs(yerr[i])/ *(yscal+i))errmax=fabs(yerr[i])/ *(yscal+i);
  }
 errmax /= eps;
  if (errmax > 1.0)
  {
   h *= SAFETY*pow(errmax,PSHRNK);
   xnew=x+h*32;
   if (xnew == x)
   {
   printf("stepsize underflow in rkqs\n");
   *hnext = -1. ;
   break ;
   }
  }
  else
  {
   if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
   else *hnext=5.0*h;
   *hdid = h ;
   x += h ;
   for (i=0;i<n;i++) *(y+i) = ytemp[i];
   break;
  }
 }
}


 void Prob_LaserExc::Rkck(G4double *y, G4double *dydx, G4int n, G4double x, G4double h,G4double *yout, G4double *yerr)
{
 G4int i;
 //Cash-Karp parameters for embedded Runge-Kutta 
 static G4double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
  b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
  b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
  c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
  dc5 = -277.0/14336.0;
 G4double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
  dc4=c4-13525.0/55296.0,dc6=c6-0.25;
 G4double ak2[40],ak3[40],ak4[40],ak5[40],ak6[40],ytemp[40];

 for (i=0;i<n;i++)
  ytemp[i]=*(y+i)+b21*h* *(dydx+i);
 cderiv(x+a2*h,ytemp,ak2);
 for (i=0;i<n;i++)
  ytemp[i]=*(y+i)+h*(b31* *(dydx+i)+b32*ak2[i]);
 cderiv(x+a3*h,ytemp,ak3);
 for (i=0;i<n;i++)
  ytemp[i]=*(y+i)+h*(b41* *(dydx+i)+b42*ak2[i]+b43*ak3[i]);
 cderiv(x+a4*h,ytemp,ak4);
 for (i=0;i<n;i++)
  ytemp[i]=*(y+i)+h*(b51* *(dydx+i)+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
 cderiv(x+a5*h,ytemp,ak5);
 for (i=0;i<n;i++)
  ytemp[i]=*(y+i)+h*(b61* *(dydx+i)+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
 cderiv(x+a6*h,ytemp,ak6);
 for (i=0;i<n;i++)
  *(yout+i)=*(y+i)+h*(c1* *(dydx+i)+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
 for (i=0;i<n;i++)
  *(yerr+i)=h*(dc1* *(dydx+i)+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}  
