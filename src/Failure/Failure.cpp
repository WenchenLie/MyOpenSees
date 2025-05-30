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
                                                                        
// $Revision: 1.0 $
// $Date: 2025-01-22 13:50:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Failure.cpp,v $

// Written: Wenchen Lie
// Created: Jan 22, 2025
//
// Description: This file contains the class definition for Failure.h

#include <stdlib.h>
#include <string.h>

#include <Failure.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

static int numFailureMaterials = 0;

OPS_Export void *
OPS_Failure(void)
{
  // Print information
  if (numFailureMaterials == 0) {
      opserr << "Failure unaxial wrapper material - Written by Wenchen Lie (Jan 22, 2025)\n";
      numFailureMaterials++;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  UniaxialMaterial *theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  double minStress = -1.0e16;
  double maxStress = 1.0e16;
  double uy = 1.0e16;
  double maxCPD = 1.0e16;
  int    iData[2];
  double jData[2];

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial Failure $tag $otherTag <-minStrain $minStrain> <-maxStrain $maxStrain> <-minStress $minStress> <-maxStress $maxStress> <-maxCPD $uy $maxCPD>" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial Failure $tag $otherTag" << endln;
    return 0;
  }

  theOtherMaterial = OPS_GetUniaxialMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag uniaxialMaterial Failure tag: " << iData[0] << endln;
    return 0;	
  }

  argc = OPS_GetNumRemainingInputArgs();  
  while (argc > 1) {
    const char *argvLoc = OPS_GetString();
    numData = 1;
    int numData2 = 2;
    if (strcmp(argvLoc, "-minStrain") == 0) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {      
	    opserr << "WARNING invalid minStrain value for uniaxialMaterial Failure tag: " << iData[0] << "\n<-minStrain $minStrain>" << endln;
	    return 0;
      }
    }
    else if (strcmp(argvLoc, "-maxStrain") == 0) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {      
	    opserr << "WARNING invalid maxStrain value for uniaxialMaterial Failure tag: " << iData[0] << "\n<-maxStrain $maxStrain>" << endln;
	    return 0;
      }
    }
    else if (strcmp(argvLoc, "-minStress") == 0) {
      if (OPS_GetDouble(&numData, &minStress) != 0) {
          opserr << "WARNING invalid minStress value for uniaxialMaterial Failure tag: " << iData[0] << "\n<-minStress $minStress>" << endln;
          return 0;
      }
    }
    else if (strcmp(argvLoc, "-maxStress") == 0) {
      if (OPS_GetDouble(&numData, &maxStress) != 0) {
          opserr << "WARNING invalid maxStress value for uniaxialMaterial Failure tag: " << iData[0] << "\n<-maxStress $maxStress>" << endln;
          return 0;
      }
    }
    else if (strcmp(argvLoc, "-maxCPD") == 0) {
        if (OPS_GetDouble(&numData2, jData) != 0) {
            opserr << "WARNING invalid uy or maxCPD value for uniaxialMaterial Failure tag: " << iData[0] << "\n<-maxCPD $uy $maxCPD>" << endln;
            return 0;
        }
        uy = jData[0];
        maxCPD = jData[1];
    }
    else {
      opserr << "WARNING invalid option:" << argvLoc << " uniaxialMaterial Failure tag: " << iData[0] << "\nFailure $tag $otherTag <-minStrain $minStrain> <-maxStrain $maxStrain> <-minStress $minStress> <-maxStress $maxStress> <-maxCPD $uy $maxCPD>" << endln;
      return 0;
    }
    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial = new Failure(iData[0], *theOtherMaterial, minStrain, maxStrain, minStress, maxStress, uy, maxCPD);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Failure\n";
    return 0;
  }

  return theMaterial;
}

Failure::Failure(int tag, UniaxialMaterial &material, double minStrain, double maxStrain,
    double minStress, double maxStress, double uy, double maxCPD)
  :UniaxialMaterial(tag,MAT_TAG_Failure), theMaterial(0),
   minStrain(minStrain), maxStrain(maxStrain), minStress(minStress), maxStress(maxStress), uy(uy), maxCPD(maxCPD),
   Twp(0.0), Cwp(0.0), Tyieldface(uy), Cyieldface(uy), Tfailed(false), Cfailed(false)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "Failure::Failure -- failed to get copy of material\n";
    exit(-1);
  }
}

Failure::Failure()
  :UniaxialMaterial(0,MAT_TAG_Failure), theMaterial(0),
   minStrain(0.0), maxStrain(0.0), minStress(0.0), maxStress(0.0), uy(0.0), maxCPD(0.0),
   Tyieldface(0.0), Cyieldface(0.0), Twp(0.0), Cwp(0.0), Tfailed(false), Cfailed(false)
{

}

Failure::~Failure()
{
  if (theMaterial)
    delete theMaterial;
}

int Failure::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed)
    return 0;
  double stress = getStress();
  if (strain > Cyieldface)
  {
    // Posotive yielding
    Twp = Cwp + strain - Cyieldface;
    Tyieldface = strain;
  }
  else if (strain < Cyieldface - 2 * uy)
  {
    // Negative yielding
    Twp = Cwp + Cyieldface - 2 * uy - strain;
    Tyieldface = strain + 2 * uy;
  }
  if (strain >= maxStrain || strain <= minStrain || stress >= maxStress || stress <= minStress) {
    Tfailed = true;  // Ductile failure
    return 0;
  }
  else if (Twp / uy >= maxCPD) {
    Tfailed = true;  // Failure due to cumulative plastic displacement
    return 0;
  }
  else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}

double Failure::getStress(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getStress();
}

double Failure::getTangent(void)
{
  if (Tfailed)
    //return 0.0;
    return 1.0e-8*theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double Failure::getDampTangent(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}

double Failure::getStrain(void)
{
  return theMaterial->getStrain();
}

double Failure::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int Failure::commitState(void)
{	
  Cfailed = Tfailed;
  Cyieldface = Tyieldface;
  Cwp = Twp;

  // Check if failed at current step
  if (Tfailed)
    return 0;
  else
    return theMaterial->commitState();
}

int Failure::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int Failure::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;
  Cyieldface = 0.0;
  Tyieldface = 0.0;
  Cwp = 0.0;
  Twp = 0.0;
  return theMaterial->revertToStart();
}

UniaxialMaterial *
Failure::getCopy(void)
{
  Failure *theCopy = 
    new Failure(this->getTag(), *theMaterial, minStrain, maxStrain, minStress, maxStress, uy, maxCPD);
        
  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  theCopy->Cyieldface = Cyieldface;
  theCopy->Tyieldface = Tyieldface;
  theCopy->Cwp = Cwp;
  theCopy->Twp = Twp;
  
  return theCopy;
}

int Failure::sendSelf(int cTag, Channel &theChannel)
{
    opserr << "Failure::sendSelf() is not avaiable for Failure material" << endln;
    return -1;
}

int Failure::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
    opserr << "Failure::recvSelf() is not avaiable for Failure material" << endln;
    return -1;
}

void Failure::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Failure, tag: " << this->getTag() << endln;
        s << "  material: " << theMaterial->getTag() << endln;
        s << "  min strain: " << minStrain << endln;
        s << "  max strain: " << maxStrain << endln;
        s << "  min stress: " << minStress << endln;
        s << "  max stress: " << maxStress << endln;
        s << "  uy        : " << uy << endln;
        s << "  maxCPD    : " << maxCPD << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"Failure\", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\", ";
        s << "\"minStrain\": " << minStrain << ", ";
        s << "\"maxStrain\": " << maxStrain << ", ";
        s << "\"minStress\": " << minStress << ", ";
        s << "\"maxStress\": " << maxStress << ", ";
        s << "\"uy\": " << uy << ", ";
        s << "\"maxCPD\": " << maxCPD << "}";
    }
}
