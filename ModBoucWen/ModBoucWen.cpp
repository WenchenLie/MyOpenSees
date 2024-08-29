/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: Wenchen Lie (666@e.gzhu.edu.cn)
// Created: July 26, 2024
// Revision: A
//
// Description: This file contains the implementation of the
// ModBoucWen class.


#include <ModBoucWen.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>

#include <math.h>
#include <float.h>


#include <elementAPI.h>
#include <OPS_Globals.h>

static int numModBoucWenMaterials = 0;

void *
OPS_ModBoucWen()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[9];
  int    iter;
  int numData = 1;

  // Print information
  if (numModBoucWenMaterials == 0) {
      opserr << "Modified BoucWen unaxial material - Written by Wenchen Lie (July 27, 2024)\n";
      numModBoucWenMaterials++;
  }

  if (OPS_GetIntInput(&numData, iData) != 0) {
      opserr << "WARNING invalid uniaxialMaterial ModBoucWen tag" << endln;
      return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  int num1 = 9;
  int num2 = 1;

  if (numData != 9 && numData != 10) {
      opserr << "Invalid args for uniaxialMaterial ModBoucWen " << iData[0] << " Fy? uy? alpha? n? Q? b? A? beta? gamma? <iter?>" << endln;
      return 0;
  }

  if (OPS_GetDoubleInput(&num1, dData) != 0) {
      opserr << "Invalid args for uniaxialMaterial ModBoucWen " << iData[0] << " Fy? uy? alpha? n? Q? b? A? beta? gamma? <iter?>" << endln;
      return 0;
  }

  if (numData == 9) {
      iter = MODBOUCWEN_DEFAULT_ITER;
  }
  else {
      if (OPS_GetIntInput(&num2, &iter) != 0) {
          opserr << "WARNING cannot get $iter, want int but got double" << endln;
          return 0;
      }
  }

  int tag = iData[0];
  double Fy, uy, alpha, n, Q, b, A, beta, gamma;
  Fy = dData[0];
  uy = dData[1];
  alpha = dData[2];
  n = dData[3];
  Q = dData[4];
  b = dData[5];
  A = dData[6];
  beta = dData[7];
  gamma = dData[8];

  if (Fy <= 0) {
      opserr << "WARNING Fy must be positive" << endln;
      return 0;
  }
  if (uy <= 0) {
      opserr << "WARNING uy must be positive" << endln;
      return 0;
  }
  if (alpha < 0) {
      opserr << "WARNING alpha must not less than 0" << endln;
      return 0;
  }
  if (n <= 0) {
      opserr << "WARNING n must be positive" << endln;
      return 0;
  }
  if (iter <= 0) {
      opserr << "WARNING iter must be positive" << endln;
      return 0;
  }
  //opserr << "Fy = " << Fy << "\n" << endln;
  //opserr << "uy = " << uy << "\n" << endln;
  //opserr << "alpha = " << alpha << "\n" << endln;
  //opserr << "n = " << n << "\n" << endln;
  //opserr << "Q = " << Q << "\n" << endln;
  //opserr << "b = " << b << "\n" << endln;
  //opserr << "beta = " << beta << "\n" << endln;
  //opserr << "gamma = " << gamma << "\n" << endln;
  //opserr << "iter = " << iter << "\n" << endln;


  // Parsing was successful, allocate the material
  theMaterial = new ModBoucWen(tag, Fy, uy, alpha, n, Q, b, A, beta, gamma, iter);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ModBoucWen Material" << endln;
    return 0;
  }

  return theMaterial;
}



ModBoucWen::ModBoucWen(int tag_, double Fy_, double uy_, double alpha_, double n_,double Q_, double b_, double A_, double beta_, double gamma_,int iter_)
    :UniaxialMaterial(tag_, MAT_TAG_ModBoucWen),
    Fy(Fy_), uy(uy_), alpha(alpha_), n(n_), Q(Q_), b(b_), A(A_), beta(beta_), gamma(gamma_), iter(iter_)
{
    Cstrain = 0.0;
    Cstress = 0.0;
    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = Fy_ / uy_;
    Ctangent = Fy_ / uy_;
    Tz = 0.0;
    Cz = 0.0;
    Twp = 0.0;
    Cwp = 0.0;
    Tface = uy;
    Cface = uy;
}

ModBoucWen::ModBoucWen():UniaxialMaterial(0, MAT_TAG_ModBoucWen),
    Fy(0.0), uy(0.0), alpha(0.0), n(0.0), Q(0.0), b(0.0), A(0.0), beta(0.0), gamma(0.0), iter(0)
{

}

ModBoucWen::~ModBoucWen ()
{
    
}

int ModBoucWen::setTrialStrain (double strain, double strainRate)
{
    Tstrain = Cstrain;
    Tstress = Cstress;
    Ttangent = Ctangent;
    Tz = Cz;
    Twp = Cwp;
    Tface = Cface;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

    if (dStrain == 0.0) {
        return 0;
    }
     // Set trial strain
    Tstrain = strain;
    double dStrain_ = dStrain / iter;
    double z_ = Cz;
    double strain_;
    double sgn;
    double m;
    double S1, S2, S3, S4;
    for (int i = 0; i < iter; i++) {
        // substep
        strain_ = Cstrain + dStrain_ * i;
        if (strain_ > Tface) {
            // positive yielding
            Twp = Twp + strain_ - Tface;
            Tface = strain_;
        }
        else if (strain_ < Tface - 2 * uy) {
            // negative yielding
            Twp = Twp + Tface - 2 * uy - strain_;
            Tface = strain_ + 2 * uy;
        }
        if (dStrain_ * z_ < 0.0) {
            sgn = -1.0;
        }
        else if (dStrain_ * z_ == 0.0) {
            sgn = 0.0;
        }
        else {
            sgn = 1.0;
        }
        m = 1 + Q * (1 - pow(b, -Twp / uy));
        S1 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(fabs(z_ / m), n));
        S2 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(fabs(z_ / m + 0.5 * dStrain_ * S1), n));
        S3 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(fabs(z_ / m + 0.5 * dStrain_ * S2), n));
        S4 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(fabs(z_ / m + dStrain_ * S3), n));
        z_ = z_ + 1.0 / 6.0 * dStrain_ * (S1 + S2 + S3 + S4);
        Tstress = alpha * Fy / uy * strain_ + (1 - alpha) * Fy * z_;
    };
    Tz = z_;
    Ttangent = (Tstress - Cstress) / dStrain;
}

double ModBoucWen::getStrain ()
{
   return Tstrain;
}

double ModBoucWen::getStress ()
{
   return Tstress;
}

double ModBoucWen::getTangent ()
{
    return Ttangent;
}

int ModBoucWen::commitState ()
{
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;
   Cz = Tz;
   Cwp = Twp;
   Cface = Tface;
   return 0;
}

int ModBoucWen::revertToLastCommit ()
{
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ttangent;
   Tz = Cz;
   Twp = Cwp;
   Tface = Cface;
   return 0;
}

int ModBoucWen::revertToStart ()
{
   Cstrain = 0.0;
   Cstress = 0.0;
   Ttangent = Fy / uy;
   Ctangent = Fy / uy;
   Tstrain = 0.0;
   Tstress = 0.0;
   Tz = 0.0;
   Twp = 0.0;
   Tface = uy;
   Cz = 0.0;
   Cwp = 0.0;
   Cface = uy;
   return 0;
}

UniaxialMaterial* ModBoucWen::getCopy ()
{
   ModBoucWen* theCopy = new ModBoucWen(this->getTag(), Fy, uy, alpha, n, Q, b, A, beta, gamma, iter);
   return theCopy;
}

int ModBoucWen::sendSelf (int commitTag, Channel& theChannel)
{
    opserr << "Currently ModBoucWen::sendSelf() is not avaiable for ModBoucWen material" << endln;
    return -1;
}

int ModBoucWen::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
    opserr << "Currently ModBoucWen::recvSelf() is not avaiable for ModBoucWen material" << endln;
    return -1;

}

void ModBoucWen::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "ModBoucWen tag:   " << this->getTag() << endln;
    s << "  Fy:         " << Fy << " ";
    s << "  uy:         " << uy << " ";
    s << "  alpha:      " << alpha << " ";
    s << "  n:          " << n << " ";
    s << "  Q:          " << Q << " ";
    s << "  b:          " << b << " ";
    s << "  A:          " << A << " ";
    s << "  beta:       " << beta << " ";
    s << "  gamma:      " << gamma << " ";
    s << "  iter:       " << iter << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"tag\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"ModBoucWen\", ";
    s << "\"Fy\": " << Fy << ", ";
    s << "\"uy\": " << uy << ", ";
    s << "\"alpha\": " << alpha << ", ";
	s << "\"n " << n << ", ";
    s << "\"Q\": " << Q << ", ";
    s << "\"b\": " << b << ", ";
    s << "\"A\": " << A << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"gamma\": " << gamma << ", ";
    s << "\"iter\": " << iter << ", ";
  }
  
}
