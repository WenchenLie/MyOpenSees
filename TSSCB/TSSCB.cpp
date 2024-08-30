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
// Last update: July 28, 2024
// Revision: A
//
// Description: This file contains the implementation of the
// TSSCB class.


#include <TSSCB.h>
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

static int numTSSCBMaterials = 0;

void *
OPS_TSSCB()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // Print information
  if (numTSSCBMaterials == 0) {
      opserr << "TSSCB unaxial material - Written by Wenchen Lie (July 26, 2024)\n";
      numTSSCBMaterials++;
  }

  int iData[1];
  double dData[7];
  int numData = 1;


  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial TSSCB tag " << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 7) {
      opserr << "Invalid #args, want: uniaxialMaterial TSSCB " << iData[0] << " F1? k0? ugap? F2? k1? k2? beta?" << endln;
      return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial TSSCB " << iData[0] << " F1? k0? ugap? F2? k1? k2? beta?" << endln;
      return 0;
  }

  // Parsing was successful, allocate the material
  double F1;  double k0;  double ugap;  double F2;  double k1;  double k2;  double beta;
  F1 = dData[0];  k0 = dData[1];  ugap = dData[2];  F2 = dData[3];  k1 = dData[4];  k2 = dData[5];  beta = dData[6];
  if (F1 <= 0) {
      opserr << "WARNING Fy should lager than 0" << endln;
      return 0;
  }
  if (k0 <= 0) {
      opserr << "WARNING k0 should lager than 0" << endln;
      return 0;
  }
  if (ugap < 0) {
      opserr << "WARNING ugap should not less than 0" << endln;
      return 0;
  }
  if (F2 <= 0) {
      opserr << "WARNING F2 should larger than 0" << endln;
      return 0;
  }
  if (k1 <= 0) {
      opserr << "WARNING k1 should larger than 0" << endln;
      return 0;
  }
  if (k2 <= 0) {
      opserr << "WARNING k2 should larger than 0" << endln;
      return 0;
  }
  if (beta < 0 || beta > 1) {
      opserr << "WARNING beta should be within [0, 1]" << endln;
      return 0;
  }
  theMaterial = new TSSCB(iData[0], F1, k0, ugap, F2, k1, k2, beta);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type TSSCB Material" << endln;
    return 0;
  }

  return theMaterial;
}



TSSCB::TSSCB(int tag_, double F1_, double k0_, double ugap_, double F2_, double k1_, double k2_, double beta_):
    UniaxialMaterial(tag_, MAT_TAG_TSSCB),
    F1(F1_), k0(k0_), ugap(ugap_), F2(F2_), k1(k1_), k2(k2_), beta(beta_)
{
    Cstrain = 0.0;
    Cstress = 0.0;
    
    Cstage = 1;
    Tstrain = 0.0;
    Tstress = 0.0;
    if (ugap != 0) {
        Ttangent = k0;
        Ctangent = k0;
    }
    else {
        Ttangent = k1;
        Ctangent = k1;
        F1 = F2;
    }
    Tstage = 1;
    ua = ugap - F1 / k0;
    if (ua < 0) {
        ua = 0;
    }
}

TSSCB::TSSCB():UniaxialMaterial(0, MAT_TAG_TSSCB),
    F1(0.0), k0(0.0), ugap(0.0), F2(0.0), k1(0.0), k2(0.0), beta(0.0)
{

}

TSSCB::~TSSCB ()
{
    
}

int TSSCB::setTrialStrain (double strain, double strainRate)
{
    // Reset history variables to last converged state
    Tstrain = Cstrain;
    Tstress = Cstress;
    Ttangent = Ctangent;
    Tstage = Cstage;

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


void TSSCB::determineTrialState (double dStrain)
{
    double F1_;
    double F2_;
    double du1;
    double du2;
    double usc0;

    // determine working stage
    if (-ugap <= Tstrain && Tstrain <= ugap) {
        Tstage = 1;
    }
    else {
        Tstage = 2;
    }

    // updata Tstress
    if (Cstage == 1 && Tstage == 1) {
        // stage-1 -> stage-1
        Tstress = frictionModel(Cstress, dStrain);
    }
    else if (Cstage == 1 && Tstage == 2) {
        // stage-1 -> stage-2
        if (dStrain > 0) {
            du1 = ugap - Cstrain;
            du2 = dStrain - du1;
            usc0 = ugap - ua;
        }
        else {
            du1 = -ugap - Cstrain;
            du2 = dStrain - du1;
            usc0 = ua - ugap;
        }
        F1_ = frictionModel(Cstress, du1);
        F2_ = SCModel(usc0, F1_, du2);
        Tstress = F2_;
    }
    else if (Cstage == 2 && Tstage == 2) {
        // stage-2 -> stage-2
        if (Cstrain >= 0) {
            usc0 = Cstrain - ua;
        }
        else {
            usc0 = Cstrain + ua;
        }
        Tstress = SCModel(usc0, Cstress, dStrain);
    }
    else if (Cstage == 2 && Tstage == 1) {
        // stage-2 -> stage-1
        if (dStrain < 0) {
            du1 = -(Cstrain - ugap);
            du2 = -(ugap - Tstrain);
            usc0 = Cstrain - ua;
        }
        else {
            du1 = -ugap - Cstrain;
            du2 = Tstrain + ugap;
            usc0 = Cstrain + ua;
        }
        F1_ = SCModel(usc0, Cstress, du1);
        F2_ = frictionModel(F1_, du2);
        Tstress = F2_;
    }
    if (dStrain != 0) {
        Ttangent = (Tstress - Cstress) / dStrain;
    }
}

double TSSCB::frictionModel(double F0, double du)
{
    if (du == 0.0) {
        return F0;
    }
    double F_;
    double F;
    F_ = F0 + du * k0;
    if (F_ > F1) {
        F = F1;
    } else if (F_ < -F1) {
        F = -F1;
    } else {
        F = F_;
    };
    return F;
}

double TSSCB::SCModel(double u0, double F0, double du)
{
    double F;
    double u;
    double uy;
    double F_;

    if (du == 0) {
        F = F0;
        return F;
    };
    u = u0 + du;
    uy = F2 / k1;
    F_ = F0 + du * k1;
    if (du > 0) {
        if (u < -F2 * (1 - 2 * beta) / k1 && F_ > k2 * u - F2 * (1 - 2 * beta) * (1 - k2 / k1)) {
            F = k2 * u - F2 * (1 - 2 * beta) * (1 - k2 / k1);
        }
        else if (-F2 * (1 - 2 * beta) / k1 <= u && u <= uy && F_ > k1 * u) {
            F = k1 * u;
        }
        else if (u > F2 * (1 - 2 * beta) / k1 && F_ > k2 * u + F2 - k2 * uy) {
            F = k2 * u + F2 - k2 * uy;
        }
        else {
            F = F_;
        }
    }
    else {
        if (u > F2 * (1 - 2 * beta) / k1 && F_ < k2 * u + F2 * (1 - 2 * beta) * (1 - k2 / k1)) {
            F = k2 * u + F2 * (1 - 2 * beta) * (1 - k2 / k1);
        }
        else if (-uy <= u && u <= F2 * (1 - 2 * beta) / k1 && F_ < k1 * u) {
            F = k1 * u;
        }
        else if (u < -F2 * (1 - 2 * beta) / k1 && F_ < k2 * u - (F2 - k2 * uy)) {
            F = k2 * u - (F2 - k2 * uy);
        }
        else {
            F = F_;
        }
    }
    return F;
}


double TSSCB::getStrain ()
{
   return Tstrain;
}

double TSSCB::getStress ()
{
   return Tstress;
}

double TSSCB::getTangent ()
{
    return Ttangent;
}

int TSSCB::commitState ()
{
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;
   Cstage = Tstage;
   return 0;
}

double TSSCB::getInitialTangent()
{
    if (ugap == 0) {
        return k1;
    }
    else {
        return k0;
    }
}



int TSSCB::revertToLastCommit ()
{
   Tstrain = Cstrain;
   Tstress = Cstress;
   return 0;
}

int TSSCB::revertToStart ()
{
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = 0.0;
   Cstage = 1;
   Tstrain = 0.0;
   Tstress = 0.0;
   if (ugap == 0) {
       Ttangent = k1;
   }
   else {
       Ttangent = k0;
   }
   Tstage = 1;
   return 0;
}

UniaxialMaterial* TSSCB::getCopy ()
{
   TSSCB* theCopy = new TSSCB(this->getTag(), F1, k0, ugap, F2, k1, k2, beta);

   return theCopy;
}

int TSSCB::sendSelf (int commitTag, Channel& theChannel)
{
    opserr << "Currently TSSCB::sendSelf() is not avaiable for TSSCB material" << endln;
    return -1;
}

int TSSCB::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
    opserr << "Currently TSSCB::recvSelf() is not avaiable for TSSCB material" << endln;
    return -1;

}

void TSSCB::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "TSSCB tag:   " << this->getTag() << endln;
    s << "  F1:      " << F1 << " ";
    s << "  k0:      " << k0 << " ";
    s << "  ugap:    " << ugap << " ";
    s << "  F2:      " << F2 << " ";
    s << "  k1:      " << k1 << " ";
    s << "  k2:      " << k2 << " ";
    s << "  beta:    " << beta << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"TSSCB\", ";
	s << "\"F1\": " << F1 << ", ";
	s << "\"k0\": " << k0 << ", ";
    s << "\"ugap\": " << ugap << ", ";
    s << "\"F2\": " << F2 << ", ";
    s << "\"k1\": " << k1 << ", ";
    s << "\"k2\": " << k2 << ", ";
    s << "\"beta\": " << beta << ", ";
  }
  
}
