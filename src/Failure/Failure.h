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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Failure.h,v $
                                                      
// Written: Wenchen Lie
// Created: Jan 22, 2025
//
// Description: This file contains the class definition for Failure.h

#ifndef Failure_h
#define Failure_h

#include <UniaxialMaterial.h>

class Failure : public UniaxialMaterial
{
  public:
    Failure(int tag, UniaxialMaterial &material, double minStrain, double maxStrain,
        double minStress, double maxStress, double uy, double maxCPD); 
    Failure();
    ~Failure();
    
    const char *getClassType(void) const {return "Failure";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStrainRate(void);
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void);
    double getInitialTangent(void) {return theMaterial->getInitialTangent();}

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    bool hasFailed(void) {return Cfailed;}

  protected:
    
  private:
	UniaxialMaterial *theMaterial;

	double minStrain;
	double maxStrain;
    double minStress;
    double maxStress;
    double uy;
    double maxCPD;
    double Tyieldface;
    double Cyieldface;
    double Twp;
    double Cwp;
	bool Tfailed;
	bool Cfailed;
};


#endif

