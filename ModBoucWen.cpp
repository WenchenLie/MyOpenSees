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

  // Print information
  if (numModBoucWenMaterials == 0) {
      opserr << "Modified BoucWen unaxial material - Written by Wenchen Lie (July 27, 2024)\n";
      numModBoucWenMaterials++;
  }

  int tag = 1;
  int numData = 1;
  int gotFailureCPD = 0;
  int gotMinMax = 0;
  double failureCPD = 0.0;
  double MinMax = 0.0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 10 && numArgs != 12 && numArgs != 14) {
      opserr << "Invalid #args, want: uniaxialMaterial ModBoucWen tag Fy uy alpha n Q b A beta gamma <-failureCPD $failureCPD>" << endln;
      return 0;
  }

  if (OPS_GetIntInput(&numData, &tag) != 0) {
      opserr << "cannot get material tag of ModBoucWen material with tag " << tag << endln;
      return 0;
  }

  numData = 9;
  double data[9];
  if (OPS_GetDoubleInput(&numData, data) != 0) {
      opserr << "cannot get material arguments of ModBoucWen material with tag " << tag << endln;
      return 0;
  }

  int argc = OPS_GetNumRemainingInputArgs();
  while (argc > 1) {
      const char* str = OPS_GetString();
      numData = 1;
      if (strcmp(str, "-failureCPD") == 0) {
          if (OPS_GetDoubleInput(&numData, &failureCPD) != 0) {
              opserr << "WARNING invalid failureCPD value" << endln;
              return 0;
          }
          gotFailureCPD = 1;
      }
      else if (strcmp(str, "-MinMax") == 0) {
          if (OPS_GetDoubleInput(&numData, &MinMax) != 0) {
              opserr << "WARNING invalid MinMax value" << endln;
              return 0;
          }
          if (MinMax <= 0) {
              opserr << "WARNING invalid MinMax value (MinMax should be larger than 0)" << endln;
              return 0;
          }
          gotMinMax = 1;
      }
      else {
          opserr << "unknow command \"" << str << "\" in ModBoucWen material with tag " << tag << endln;
          opserr << "want \"-failureCPD\" or \"-MinMax\"\n" << tag << endln;
          return 0;
      }
      argc = OPS_GetNumRemainingInputArgs();
  }


  double Fy, uy, alpha, n, Q, b, A, beta, gamma;
  Fy = data[0];
  uy = data[1];
  alpha = data[2];
  n = data[3];
  Q = data[4];
  b = data[5];
  A = data[6];
  beta = data[7];
  gamma = data[8];

  //opserr << "Fy = " << Fy << "\n" << endln;
  //opserr << "uy = " << uy << "\n" << endln;
  //opserr << "alpha = " << alpha << "\n" << endln;
  //opserr << "n = " << n << "\n" << endln;
  //opserr << "Q = " << Q << "\n" << endln;
  //opserr << "b = " << b << "\n" << endln;
  //opserr << "A = " << A << "\n" << endln;
  //opserr << "beta = " << beta << "\n" << endln;
  //opserr << "gamma = " << gamma << "\n" << endln;
  //opserr << "gotFailureCPD = " << gotFailureCPD << "\n" << endln;
  //opserr << "failureCPD = " << failureCPD << "\n" << endln;
  //opserr << "gotMinMax = " << gotMinMax << "\n" << endln;
  //opserr << "MinMax = " << MinMax << "\n" << endln;


  // Parsing was successful, allocate the material
  theMaterial = new ModBoucWen(tag, Fy, uy, alpha, n, Q, b, A, beta, gamma, gotFailureCPD, failureCPD, gotMinMax, MinMax);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ModBoucWen Material" << endln;
    return 0;
  }

  return theMaterial;
}



ModBoucWen::ModBoucWen(int tag_, double Fy_, double uy_, double alpha_, double n_,
    double Q_, double b_, double A_, double beta_, double gamma_,
    int gotFailureCPD_, double failureCPD_, int gotMinMax_, double MinMax_)
    :UniaxialMaterial(tag_, MAT_TAG_ModBoucWen),
    Fy(Fy_), uy(uy_), alpha(alpha_), n(n_),
    Q(Q_), b(b_), A(A_), beta(beta_), gamma(gamma_),
    failure(false),
    gotFailureCPD(gotFailureCPD_), failureCPD(failureCPD_),
    gotMinMax(gotMinMax_), MinMax(MinMax_)
{
    Cstrain = 0.0;
    Cstress = 0.0;
    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = Fy_ / uy_;
    Tz = 0.0;
    Cz = 0.0;
    Twp = 0.0;
    Cwp = 0.0;
    Tface = uy;
    Cface = uy;
}

ModBoucWen::ModBoucWen():UniaxialMaterial(0, MAT_TAG_ModBoucWen),
    Fy(0.0), uy(0.0), alpha(0.0), n(0.0),
    Q(0.0), b(0.0), A(0.0), beta(0.0), gamma(0.0),
    failure(false),
    gotFailureCPD(0), failureCPD(0.0),
    gotMinMax(0), MinMax(0.0)
{

}

ModBoucWen::~ModBoucWen ()
{
    
}

int ModBoucWen::setTrialStrain (double strain, double strainRate)
{
   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   return 0;
}

int ModBoucWen::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void ModBoucWen::determineTrialState (double dStrain)
{
    // get trail stress and update wp and z

    double sgn;
    double m;
    double S1, S2, S3, S4;

    if (Tstrain > Cface) {
        Twp = Cwp + Tstrain - Cface;
        Tface = Tstrain;
    }
    else if (Tstrain < Cface - 2 * uy) {
        Twp = Cwp + Cface - 2 * uy - Tstrain;
        Tface = Tstrain + 2 * uy;
    }

    if (gotFailureCPD && Twp / uy >= failureCPD) {
        failure = true;
    }
    if (gotMinMax && abs(Tstrain) >= MinMax) {
        failure = true;
    }
    if (failure) {
        Tstress = 0.0;
    }
    else {
        if (dStrain * Cz < 0.0) {
            sgn = -1.0;
        }
        else if (dStrain * Cz == 0.0) {
            sgn = 0.0;
        }
        else {
            sgn = 1.0;
        }
        m = 1 + Q * (1 - pow(b, -Twp / uy));
        S1 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(abs(Cz / m), n));
        S2 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(abs(Cz / m + 0.5 * dStrain * S1), n));
        S3 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(abs(Cz / m + 0.5 * dStrain * S2), n));
        S4 = 1.0 / uy * (A - (beta * sgn + gamma) * pow(abs(Cz / m + dStrain * S3), n));
        Tz = Cz + 1.0 / 6.0 * dStrain * (S1 + S2 + S3 + S4);  // update z
        Tstress = alpha * Fy / uy * Tstrain + (1 - alpha) * Fy * Tz;
    };
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
    if (fabs(Tstrain - Cstrain) > DBL_EPSILON) {
        Ttangent = (Tstress - Cstress) / (Tstrain - Cstrain);
    }
    return Ttangent;
}

int ModBoucWen::commitState ()
{
   Cstrain = Tstrain;
   Cstress = Tstress;
   Cz = Tz;
   Cwp = Twp;
   Cface = Tface;
   return 0;
}

int ModBoucWen::revertToLastCommit ()
{
   Tstrain = Cstrain;
   Tstress = Cstress;
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
   ModBoucWen* theCopy = new ModBoucWen(this->getTag(), Fy, uy, alpha, n, Q, b, A, beta, gamma, gotFailureCPD, failureCPD, gotMinMax, MinMax);
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
    s << "  failureCPD: " << failureCPD << " ";
    s << "  MinMax:     " << MinMax << " ";
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
    s << "\"failureCPD\": " << failureCPD << ", ";
    s << "\"MinMax\": " << MinMax << ", ";
  }
  
}
