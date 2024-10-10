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
** file 'COPYRIGHT'  in main directory for information on usage &&   **
** redistribution,  && for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: Wenchen Lie (666@e.gzhu.edu.cn)
// Created: July 26, 2024
// Last update: Oct 10, 2024
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
      opserr << "TSSCB unaxial material - Written by Wenchen Lie (Oct 10, 2024)\n";
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
  if (numData != 7 && numData != 9 && numData != 12 && numData != 14) {
      opserr << "Invalid #args, want 7, 9, 12, or 14 arguments but got " << numData << " (F1 k0 ugap F2 k1 k2 beta <-hardening uh r1 r2 r3> <-minmax uf>)" << endln;
      return 0;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial TSSCB " << iData[0] << " F1 k0 ugap F2 k1 k2 beta <-hardening uh r1 r2 r3> <-minmax uf>" << endln;
      return 0;
  }

  double hardening_paras[4] = { 1.0e16, 1.0, 1.0, 0.0 };
  double uf = 1.0e16;
  int argc = OPS_GetNumRemainingInputArgs();
  while (argc > 1) {
      const char* argvLoc = OPS_GetString();
      if (strcmp(argvLoc, "-hardening") == 0) {
          numData = 4;
          if (OPS_GetDouble(&numData, hardening_paras) != 0) {
              opserr << "Invalid #args, want: uniaxialMaterial TSSCB " << iData[0] << " F1 k0 ugap F2 k1 k2 beta <-hardening uh r1 r2 r3> <-minmax uf>" << endln;
              return 0;
          }
      }
      else if (strcmp(argvLoc, "-minmax") == 0) {
          numData = 1;
          if (OPS_GetDouble(&numData, &uf) != 0) {
              opserr << "Invalid #args, want: uniaxialMaterial TSSCB " << iData[0] << " F1 k0 ugap F2 k1 k2 beta <-hardening uh r1 r2 r3> <-minmax uf>" << endln;
              return 0;
          }
      }
      else {
          opserr << "Invalid #args, want: uniaxialMaterial TSSCB " << iData[0] << " F1 k0 ugap F2 k1 k2 beta <-hardening uh r1 r2 r3> <-minmax uf>" << endln;
          return 0;
      }
      argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  double F1;  double k0;  double ugap;  double F2;  double k1;  double k2;  double beta; double uh; double r1; double r2; double r3;
  F1 = dData[0];  k0 = dData[1];  ugap = dData[2];  F2 = dData[3];  k1 = dData[4];  k2 = dData[5];  beta = dData[6];
  uh = hardening_paras[0]; r1 = hardening_paras[1]; r2 = hardening_paras[2]; r3 = hardening_paras[3];
  if (F1 < 0) {
      opserr << "WARNING Fy should not less than 0" << endln;
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
  if (uh <= 0) {
      opserr << "WARNING uh should be larger than 0" << endln;
      return 0;
  }
  if (r1 < 0) {
      opserr << "WARNING r1 should not less than 0" << endln;
      return 0;
  }
  if (r2 < 0) {
      opserr << "WARNING r2 should not less than 0" << endln;
      return 0;
  }
  if (r3 < 0) {
      opserr << "WARNING r3 should not less than 0" << endln;
      return 0;
  }
  if (uf <= 0) {
      opserr << "WARNING uf should be larger than 0" << endln;
      return 0;
  }
  //opserr << "F1 = " << F1 << endln;
  //opserr << "k0 = " << k0 << endln;
  //opserr << "ugap = " << ugap << endln;
  //opserr << "F2 = " << F2 << endln;
  //opserr << "k1 = " << k1 << endln;
  //opserr << "k2 = " << k2 << endln;
  //opserr << "beta = " << beta << endln;
  //opserr << "uh = " << uh << endln;
  //opserr << "r1 = " << r1 << endln;
  //opserr << "r2 = " << r2 << endln;
  //opserr << "r3 = " << r3 << endln;
  //opserr << "uf = " << uf << endln;
  theMaterial = new TSSCB(iData[0], F1, k0, ugap, F2, k1, k2, beta, uh, r1, r2, r3, uf);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type TSSCB Material" << endln;
    return 0;
  }
  return theMaterial;
}



TSSCB::TSSCB(int tag_, double F1_, double k0_, double ugap_, double F2_, double k1_, double k2_, double beta_,
    double uh_, double r1_, double r2_, double r3_, double uf_) :
    UniaxialMaterial(tag_, MAT_TAG_TSSCB),
    F1(F1_), k0(k0_), ugap(ugap_), F2(F2_), k1(k1_), k2(k2_), beta(beta_), uh(uh_), r1(r1_), r2(r2_), r3(r3_), uf(uf_)
{
    Cstrain = 0.0;
    Tstrain = 0.0;
    Cstress3 = 0.0;
    Tstress3 = 0.0;
    Cstage = 1;
    Tstage = 1;
    if (ugap == 0) {
        Ctangent = k1;
        Ttangent = k1;
        Cstage = 2;
        Tstage = 2;
        // F1 = F2;
    }
    else {
        Ctangent = k0;
        Ttangent = k0;
    }
    Chardening = false;
    Thardening = false;
    Cstress1 = 0.0;
    Tstress1 = 0.0;
    Cstress2 = 0.0;
    Tstress2 = 0.0;
    Cstress4 = 0.0;
    Tstress4 = 0.0;
    CCDD = 0.0;
    TCDD = 0.0;
    Cfracture = false;
    Tfracture = false;
    Cplate1 = ugap;
    Tplate1 = ugap;
    Cplate2 = -ugap;
    Tplate2 = -ugap;
    ua = ugap - F1 / k1;
    if (ua < 0) {
        ua = 0.0;
    }
}

TSSCB::TSSCB() :UniaxialMaterial(0, MAT_TAG_TSSCB),
    F1(0.0), k0(0.0), ugap(0.0), F2(0.0), k1(0.0), k2(0.0), beta(0.0), ua(0.0),
    uh(1.0e16), r1(0.0), r2(0.0), r3(0.0), uf(1.0e16),
    Tstage(1), Cstage(1), Tstrain(0.0), Cstrain(0.0), Tstress3(0.0), Cstress3(0.0), Ttangent(0.0), Ctangent(0.0),
    Chardening(false), Thardening(false), Cstress1(0.0), Tstress1(0.0), Cstress2(0.0), Tstress2(0.0),
    Cstress4(0.0), Tstress4(0.0),
    CCDD(0.0), TCDD(0.0), Cfracture(false), Tfracture(false), Cplate1(ugap), Tplate1(ugap), Cplate2(-ugap), Tplate2(-ugap)
{
    
}

TSSCB::~TSSCB ()
{
    
}

int TSSCB::setTrialStrain (double strain, double strainRate)
{
    // Reset history variables to last converged state
    Tstrain = Cstrain;
    Tstress3 = Cstress3;
    Ttangent = Ctangent;
    Tstage = Cstage;
    Tstage = Cstage;
    Thardening = Chardening;
    Tstress1 = Cstress1;
    Tstress2 = Cstress2;
    Tstress4 = Cstress4;
    TCDD = CCDD;
    Tfracture = Cfracture;
    Tplate1 = Cplate1;
    Tplate2 = Cplate2;
    double dStrain = strain - Cstrain;
    if (fabs(dStrain) > DBL_EPSILON) {
        Tstrain = strain;
        // Determine whether to start hardening
        if (fabs(Tstrain) > uh || Chardening) {
            Thardening = true;
        }
        if (fabs(Tstrain) > uf) {
            Tfracture = true;
        }
        // Update Tstress
        determineTrialState(dStrain);
        // Update endplate position
        if (dStrain > 0) {
            Tplate1 = fmax(Tplate1, Tstrain);
            if (!Tfracture) {
                Tplate2 += dStrain;
            }
        }
        else {
            Tplate2 = fmin(Tplate2, Tstrain);
            if (!Tfracture) {
                Tplate1 += dStrain;
            }
        }
        if (Tplate1 < ugap) {
            Tplate1 = ugap;
        }
        if (Tplate2 > -ugap) {
            Tplate2 = -ugap;
        }
        // Calculate tangent stiffness
        Ttangent = (Tstress4 - Cstress4) / dStrain;
   }
   return 0;
}


void TSSCB::determineTrialState (double dStrain)
{
    double F1_;
    double F2_;
    double F1_ideal1;
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
    if (ugap == 0) {
        Tstage = 2;  // If ugap is zero, always in stage-2
    }
    if (Tfracture) {
        if (Tplate2 <= Tstrain && Tstrain <= Tplate1) {
            Tstress4 = 0.0;
        }
        else {
            if (dStrain > 0 && Tstrain > 0) {
                Tstress4 = F1;
            }
            else if (dStrain < 0 && Tstrain > 0) {
                Tstress4 = 0.0;
            }
            else if (dStrain > 0 && Tstrain < 0) {
                Tstress4 = 0.0;
            }
            else {
                Tstress4 = -F1;
            }
        }
        return;
    }
    // updata Tstress
    if (Cstage == 1 && Tstage == 1) {
        // stage-1 -> stage-1
        Tstress1 = frictionModel(Cstress3, dStrain);
        Tstress2 = Tstress1;
        Tstress3 = Tstress1;
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
        if (Thardening) {
            TCDD = CCDD + fabs(du2) / (uh - ugap);
        }
        F1_ = frictionModel(Cstress1, du1);
        F2_ = SCModel(usc0, F1_, du2);
        Tstress3 = F2_;
        // Apply degradation
        Tstress2 = Tstress1;
        if (Thardening && Tstrain > 0) {
            Tstress2 = Tstress1 - (F2 - F1 / 2) * TCDD * (r1 - r2 * (fabs(Tstrain) - ugap) / (uh - ugap));
        }
        else if (Thardening && Tstrain < 0) {
            Tstress2 = Tstress1 + (F2 - F1 / 2) * TCDD * (r1 - r2 * (fabs(Tstrain) - ugap) / (uh - ugap));
        }
        // Apply modifiction
        Tstress3 = Tstress2;
        if (dStrain > 0 && Tstress2 < F1) {
            Tstress3 = F1;
        }
        else if (dStrain < 0 && Tstress2 > -F1) {
            Tstress3 = -F1;
        }
    }
    else if (Cstage == 2 && Tstage == 2) {
        // stage-2 -> stage-2
        if (Thardening) {
            TCDD = CCDD + fabs(dStrain) / (uh - ugap);
        }
        if (Tstrain >= 0) {
            usc0 = Cstrain - ua;
        }
        else {
            usc0 = Cstrain + ua;
        }
        Tstress1 = SCModel(usc0, Cstress1, dStrain);
        // Apply degradation
        Tstress2 = Tstress1;
        if (Thardening && Tstrain > 0) {
            Tstress2 = Tstress1 - (F2 - F1 / 2) * TCDD * (r1 - r2 * (fabs(Tstrain) - ugap) / (uh - ugap));
        }
        else if (Thardening && Tstrain < 0) {
            Tstress2 = Tstress1 + (F2 - F1 / 2) * TCDD * (r1 - r2 * (fabs(Tstrain) - ugap) / (uh - ugap));
        }
        // Apply modifiction
        Tstress3 = Tstress2;
        if (dStrain > 0 && Tstrain > 0 && Tstress2 < F1 && ugap > 0 && Cstress3 == F1) {
            Tstress3 = F1;
        }
        else if (dStrain < 0 && Tstrain < 0 && Tstress2 > -F1 && ugap > 0 && Cstress3 == -F1) {
            Tstress3 = -F1;
        }
        else if (Tstrain > 0 && Tstress2 < 0) {
            Tstress3 = 0;  // Prevent positive compressive stress in SMA cables
        }
        else if (Tstrain < 0 && Tstress2 > 0) {
            Tstress3 = 0;  // Prevent negative compressive stress in SMA cables;
        }
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
        if (Thardening) {
            TCDD = CCDD + fabs(du1) / (uh - ugap);
        }
        F1_ = SCModel(usc0, Cstress1, du1);
        // Apply degradation
        F1_ideal1 = F1_;
        if (Thardening && Tstrain > 0) {
            F1_ideal1 = F1_ - (F2 - F1 / 2) * TCDD * (r1 - r2 * (fabs(Tstrain) - ugap) / (uh - ugap));
        }
        else if (Thardening && Tstrain < 0) {
            F1_ideal1 = F1_ + (F2 - F1 / 2) * TCDD * (r1 - r2 * (fabs(Tstrain) - ugap) / (uh - ugap));
        }
        // Apply modifiction
        F1_ = F1_ideal1;
        if (dStrain > 0 && Tstrain > 0 && F1_ideal1 < F1 && ugap > 0 && Cstress3 == F1) {
            F1_ = F1;
        }
        else if (dStrain < 0 && Tstrain < 0 && F1_ideal1 > -F1 && ugap > 0 && Cstress3 == -F1) {
            F1_ = -F1;
        }
        else if (Tstrain > 0 && F1_ideal1 < 0) {
            F1_ = 0;  // Prevent positive compressive stress in SMA cables
        }
        else if (Tstrain < 0 && F1_ideal1 > 0) {
            F1_ = 0;  // Prevent negative compressive stress in SMA cables
        }
        F2_ = frictionModel(F1_, du2);
        Tstress1 = F2_;
        Tstress2 = Tstress1;
        Tstress3 = Tstress1;
    }
    double F_hardening = fmax(fabs(Tstrain) - uh, 0) * k2 * r3;  // Strength enhancement due to hardening
    Tstress4 = Tstress3;
    if (Tstrain > 0) {
        Tstress4 += F_hardening;
    }
    else {
        Tstress4 -= F_hardening;
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
    if (Tfracture) {
        return 0.0;
    }
    u = u0 + du;
    uy = F2 / k1;
    F_ = F0 + du * k1;
    if (du > 0) {
        if (u < -F2 * (1 - beta) / k1 && F_ > k2 * u - F2 * (1 - beta) * (1 - k2 / k1)) {
            F = k2 * u - F2 * (1 - beta) * (1 - k2 / k1);
        }
        else if (-F2 * (1 - beta) / k1 <= u && u <= uy && F_ > k1 * u) {
            F = k1 * u;
        }
        else if (u > F2 * (1 - beta) / k1 && F_ > k2 * u + F2 - k2 * uy) {
            F = k2 * u + F2 - k2 * uy;
        }
        else {
            F = F_;
        }
    }
    else {
        if (u > F2 * (1 - beta) / k1 && F_ < k2 * u + F2 * (1 - beta) * (1 - k2 / k1)) {
            F = k2 * u + F2 * (1 - beta) * (1 - k2 / k1);
        }
        else if (-uy <= u && u <= F2 * (1 - beta) / k1 && F_ < k1 * u) {
            F = k1 * u;
        }
        else if (u < -F2 * (1 - beta) / k1 && F_ < k2 * u - (F2 - k2 * uy)) {
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
   return Tstress4;
}

double TSSCB::getTangent ()
{
    return Ttangent;
}

int TSSCB::commitState ()
{
   Cstrain = Tstrain;
   Ctangent = Ttangent;
   Cstage = Tstage;
   Chardening = Thardening;
   Cstress1 = Tstress1;
   Cstress2 = Tstress2;
   Cstress3 = Tstress3;
   Cstress4 = Tstress4;
   CCDD = TCDD;
   Cfracture = Tfracture;
   Cplate1 = Tplate1;
   Cplate2 = Tplate2;
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
    Tstress3 = Cstress3;
    Ttangent = Ctangent;
    Tstage = Cstage;
    Tstage = Cstage;
    Thardening = Chardening;
    Tstress1 = Cstress1;
    Tstress2 = Cstress2;
    Tstress4 = Cstress4;
    TCDD = CCDD;
    Tfracture = Cfracture;
    Tplate1 = Cplate1;
    Tplate2 = Cplate2;
    return 0;
}

int TSSCB::revertToStart ()
{
   Cstrain = 0.0;
   Cstress3 = 0.0;
   Cstage = 1;
   Tstrain = 0.0;
   Tstress3 = 0.0;
   Tstage = 1;
   if (ugap == 0) {
       Ttangent = k1;
       Ctangent = k1;
       Cstage = 2;
       Tstage = 2;
   }
   else {
       Ttangent = k0;
       Ctangent = k0;
   }
   Chardening = false;
   Thardening = false;
   Cstress1 = 0.0;
   Tstress1 = 0.0;
   Cstress2 = 0.0;
   Tstress2 = 0.0;
   Cstress4 = 0.0;
   Tstress4 = 0.0;
   CCDD = 0;  // Dimensionless cumulative damage deformation
   TCDD = 0;
   Cfracture = false;  // Whether the cable is fractured
   Tfracture = false;
   Cplate1 = ugap;  // Position of left end plate
   Tplate1 = ugap;
   Cplate2 = -ugap;  // Position of right end plate
   Tplate2 = -ugap;
   return 0;
}

UniaxialMaterial* TSSCB::getCopy ()
{
   TSSCB* theCopy = new TSSCB(this->getTag(), F1, k0, ugap, F2, k1, k2, beta, uh, r1, r2, r3, uf);

   return theCopy;
}

int TSSCB::sendSelf (int commitTag, Channel& theChannel)
{
    opserr << "TSSCB::sendSelf() is not avaiable for TSSCB material" << endln;
    return -1;
}

int TSSCB::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
    opserr << "TSSCB::recvSelf() is not avaiable for TSSCB material" << endln;
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
    s << "  uh:      " << uh << " ";
    s << "  r1:      " << r1 << " ";
    s << "  r2:      " << r2 << " ";
    s << "  r3:      " << r3 << " ";
    s << "  uf:      " << uf << " ";

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
    s << "\"uh\": " << uh << ", ";
    s << "\"r1\": " << r1 << ", ";
    s << "\"r2\": " << r2 << ", ";
    s << "\"r3\": " << r3 << ", ";
    s << "\"uf\": " << uf << ", ";
  }
  
}
