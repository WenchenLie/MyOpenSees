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
// TSB class.


#include <TSB.h>
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

static int numTSBMaterials = 0;

void *
OPS_TSB()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  // Print information
  if (numTSBMaterials == 0) {
      opserr << "TSB unaxial material - Written by Wenchen Lie (July 26, 2024)\n";
      numTSBMaterials++;
  }

  int tag;
  double friction_args[3];
  int N = 0;
  double sc_args[60];
  int numData = 1;


  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "cannot get material tag of TSB material with tag " << tag << endln;
    return 0;
  }

  int numArgs = OPS_GetNumRemainingInputArgs();
  if ((numArgs - 4) % 6 != 0) {
      opserr << "TSB material " << tag << " has a wrong number of args" << endln;
      return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, friction_args) != 0) {
    opserr << "cannot get material arguments of TSB material with tag " << tag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &N) != 0) {
      opserr << "cannot get material arguments of TSB material " << tag << endln;
      return 0;
  }
  if (N == 0 || N > 10) {
      opserr << "N of TSB material with tag " << tag << " should within 1-10 (" << N << ")" << endln;
      return 0;
  }

  numData = N * 6;
  if (OPS_GetDoubleInput(&numData, sc_args) != 0) {
      opserr << "cannot get material arguments of TSB material " << tag << endln;
      return 0;
  }
  Vector sc_args_(sc_args, numData);

  // Parsing was successful, allocate the material
  theMaterial = new TSB(tag, friction_args[0], friction_args[1], friction_args[2], N, sc_args_);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type TSB Material" << endln;
    return 0;
  }

  return theMaterial;
}



TSB::TSB(int tag_, double Fslip_, double k_, double ugap_, int N_, const Vector &sc_args_)
    :UniaxialMaterial(tag_, MAT_TAG_TSB),
    Fslip(Fslip_), k(k_), ugap(ugap_), N(N_), sc_args(sc_args_),
    Tstress(0.0), Ttangent(0.0), fracture(false), F0_sc(Vector (N_)),
    stage(1)
{
    Cstrain = 0.0;
    Cstress = 0.0;
    Tstrain = 0.0;
    Tstress = 0.0;
    Ttangent = k;

}

TSB::TSB():UniaxialMaterial(0, MAT_TAG_TSB),
    Fslip(0.0), k(0.0), ugap(0.0), N(1), sc_args(2),
    Tstress(0.0), Ttangent(0.0), fracture(false), F0_sc(Vector(1)),
    stage(1)
{

}

TSB::~TSB ()
{
    
}

int TSB::setTrialStrain (double strain, double strainRate)
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

int TSB::setTrial (double strain, double &stress, double &tangent, double strainRate)
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

void TSB::determineTrialState (double dStrain)
{
    double k_total = 0.0;  // total stiffness of self-centering components
    double ua;
    double dStrain_f;
    double dStrain_sc;
    double usc0;
    double Fsc0;
    double F_sc;
    double Tstress_sc_i;
    double FTstress;

    std::vector<double> Fy_ls(N);
    for (int i = 0; i < N; ++i) {
        Fy_ls[i] = sc_args[6 * i];
    }
    std::vector<double> k1_ls(N);
    for (int i = 0; i < N; ++i) {
        k1_ls[i] = sc_args[6 * i + 1];
    }
    std::vector<double> k2_ls(N);
    for (int i = 0; i < N; ++i) {
        k2_ls[i] = sc_args[6 * i + 2];
    }
    std::vector<double> beta_ls(N);
    for (int i = 0; i < N; ++i) {
        beta_ls[i] = sc_args[6 * i + 3];
    }
    std::vector<double> ubear_ls(N);
    for (int i = 0; i < N; ++i) {
        ubear_ls[i] = sc_args[6 * i + 4];
    }
    std::vector<double> kbear_ls(N);
    for (int i = 0; i < N; ++i) {
        kbear_ls[i] = sc_args[6 * i + 5];
    }

    for (int i = 0; i < k1_ls.size(); ++i) {
        k_total += k1_ls[i];
    }
    
    ua = ugap - Fslip / k_total;
    if (ua < 0) {
        ua = 0;
    }
    if (fabs(Cstrain) <= ugap && fabs(Tstrain) <= ugap) {
        // stage-1 -> stage-1
        Tstress = frictionModel(Cstrain, Cstress, dStrain);
        stage = 1;
        std::vector<double> F0_sc(N, 0.0);
    }
    else if (fabs(Cstrain) <= ugap && fabs(Tstrain) > ugap) {
        // stage-1 -> stage-2
        if (dStrain > 0) {
            dStrain_f = ugap - Cstrain;
            dStrain_sc = Tstrain - ugap;
            usc0 = ugap - ua;
        }
        else {
            dStrain_f = Cstrain - ugap;
            dStrain_sc = Tstrain + ugap;
            usc0 = ua - ugap;
        }
        Fsc0 = frictionModel(Cstrain, Cstress, dStrain_f);
        Tstress = 0.0;
        for (int i = 0; i < N; ++i) {
            F_sc = SCmodel(usc0, Fsc0, dStrain_sc, Fy_ls[i], k1_ls[i], k2_ls[i], beta_ls[i], ubear_ls[i] - ua, kbear_ls[i]);
            F0_sc[i] = F_sc;
            Tstress += F_sc;
        }
        stage = 2;
    }
    else if (fabs(Cstrain) > ugap && fabs(Tstrain) > ugap) {
        // stage-2 -> stage-2
        if (Cstrain >= 0) {
            usc0 = Cstrain - ua;
        }
        else {
            usc0 = Cstrain + ua;
        }
        Tstress = 0.0;
        for (int i = 0; i < N; ++i) {
            Tstress_sc_i = F0_sc[i];
            F_sc = SCmodel(usc0, Tstress_sc_i, dStrain, Fy_ls[i], k1_ls[i], k2_ls[i], beta_ls[i], ubear_ls[i] - ua, kbear_ls[i]);
            F0_sc[i] = F_sc;
            Tstress += F_sc;
        }
        stage = 2;
    }
    else if (fabs(Cstrain) > ugap && fabs(Tstrain) <= ugap) {
        // stage-2 -> stage-1
        if (dStrain < 0) {
            dStrain_sc = -(Cstrain - ugap);
            dStrain_f = -(ugap - Tstrain);
            usc0 = Cstrain - ua;
        }
        else {
            dStrain_sc = -ugap - Cstrain;
            dStrain_f = Tstrain + ugap;
            usc0 = Cstrain + ua;
        }
        FTstress = 0;
        for (int i = 0; i < N; ++i) {
            Tstress_sc_i = F0_sc[i];
            F_sc = SCmodel(usc0, Tstress_sc_i, dStrain_sc, Fy_ls[i], k1_ls[i], k2_ls[i], beta_ls[i], ubear_ls[i] - ua, kbear_ls[i]);
            F0_sc[i] = F_sc;
            FTstress += F_sc;
        }
        Tstress = frictionModel(ugap, FTstress, dStrain_f);
        stage = 1;
        std::vector<double> F0_sc(N, 0.0);
    }
}

double TSB::frictionModel(double u0, double F0, double du)
{
    if (fabs(du) <= DBL_EPSILON) {
        return F0;
    }
    double F_;
    double F;
    F_ = F0 + du * k;
    if (F_ > Fslip) {
        F = Fslip;
    } else if (F_ < -Fslip) {
        F = -Fslip;
    } else {
        F = F_;
    };
    return F;
}

double TSB::SCmodel(double u0, double F0, double du, double Fy, double k1, double k2, double beta, double ubear, double kbear)
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
    uy = Fy / k1;
    F_ = F0 + du * k1;
    if (u >= ubear) {
        F = kbear * (u - ubear) + Fy + (ubear - uy) * k2;
        return F;
    }
    else if (u <= -ubear) {
        F = kbear * (u + ubear) - Fy - (ubear - uy) * k2;
        return F;
    }
    if (du > 0) {
        if (u < -Fy * (1 - 2 * beta) / k1 && F_ > k2 * u - Fy * (1 - 2 * beta) * (1 - k2 / k1)) {
            F = k2 * u - Fy * (1 - 2 * beta) * (1 - k2 / k1);
        }
        else if (-Fy * (1 - 2 * beta) / k1 <= u && u <= uy && F_ > k1 * u) {
            F = k1 * u;
        }
        else if (u > Fy * (1 - 2 * beta) / k1 && F_ > k2 * u + Fy - k2 * uy) {
            F = k2 * u + Fy - k2 * uy;
        }
        else {
            F = F_;
        }
    }
    else {
        if (u > Fy * (1 - 2 * beta) / k1 && F_ < k2 * u + Fy * (1 - 2 * beta) * (1 - k2 / k1)) {
            F = k2 * u + Fy * (1 - 2 * beta) * (1 - k2 / k1);
        }
        else if (-uy <= u && u <= Fy * (1 - 2 * beta) / k1 && F_ < k1 * u) {
            F = k1 * u;
        }
        else if (u < -Fy * (1 - 2 * beta) / k1 && F_ < k2 * u - (Fy - k2 * uy)) {
            F = k2 * u - (Fy - k2 * uy);
        }
        else {
            F = F_;
        }
    }
    return F;
}



double TSB::getStrain ()
{
   return Tstrain;
}

double TSB::getStress ()
{
   return Tstress;
}

double TSB::getTangent ()
{
    if (fabs(Tstrain - Cstrain) > DBL_EPSILON) {
        Ttangent = (Tstress - Cstress) / (Tstrain - Cstrain);
    }
    return Ttangent;
}

int TSB::commitState ()
{
   Cstrain = Tstrain;
   Cstress = Tstress;
   return 0;
}

int TSB::revertToLastCommit ()
{
   Tstrain = Cstrain;
   Tstress = Cstress;
   return 0;
}

int TSB::revertToStart ()
{
   Cstrain = 0.0;
   Cstress = 0.0;
   Ttangent = k;
   Tstrain = 0.0;
   Tstress = 0.0;
   return 0;
}

UniaxialMaterial* TSB::getCopy ()
{
   TSB* theCopy = new TSB(this->getTag(), Fslip, k, ugap, N, sc_args);

   return theCopy;
}

int TSB::sendSelf (int commitTag, Channel& theChannel)
{
    opserr << "Currently TSB::sendSelf() is not avaiable for TSB material" << endln;
    return -1;
}

int TSB::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
    opserr << "Currently TSB::recvSelf() is not avaiable for TSB material" << endln;
    return -1;

}

void TSB::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "TSB tag:   " << this->getTag() << endln;
    s << "  Fslip:   " << Fslip << " ";
    s << "  k:       " << k << " ";
    s << "  ugap:    " << ugap << " ";
    s << "  N:       " << N << " ";
    s << "  sc_args: " << sc_args << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"TSB\", ";
	s << "\"Fslip\": " << Fslip << ", ";
	s << "\"k\": " << k << ", ";
    s << "\"ugap\": " << ugap << ", ";
    s << "\"sc_args\": " << sc_args << ", ";
  }
  
}
