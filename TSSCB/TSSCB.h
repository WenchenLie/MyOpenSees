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
// $Date: 2024-07-26 00:40:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TSSCB.h,v $
                                                                        
                                                                        
#ifndef TSSCB_h
#define TSSCB_h

// Written: Wenchen Lie 
// Created: July 26, 2024
// Last update: Sep 22, 2024
// Revision: A
//
// Description: This file contains the class definition for 
// TSSCB.h
// 
//
//
// What: "@(#) TSSCB.h, revA"


#include <UniaxialMaterial.h>


class TSSCB : public UniaxialMaterial
{
  public:
    TSSCB(int tag, double F1, double k0, double ugap, double F2, double k1, double k2, double beta, double uh, double r1, double r2, double r3, double uf);
    TSSCB();
    ~TSSCB();

    const char *getClassType(void) const {return "TSSCB";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    

 protected:
    
 private:
	/*** Material Properties ***/
    double F1;      // Friction slipping force
    double k0;      // Initial stiffness
    double ugap;    // Gap length
    double F2;      // Self-centering force at stage-II
    double k1;      // First stiffness at stage-II
    double k2;      // Second stiffness at stage-II
    double beta;    // Energy dissipation coeffience
    double uh;      // Displacement where hardening starts
    double r1;      // Coeffience to control strength degradation at the beginning of stage-2
    double r2;      // Coeffience to control strength degradation at the end of stage-2
    double r3;      // Coeffience to control stiffness enhancement due to hardening
    double uf;      // Fracture deformation

    /*** State variables ***/
    double ua;
    double Tstrain;   // strain
    double Ttangent;  // tangent stiffness
    double Tstage;    // working stage
    bool Thardening;  // whether hardening is triggered
    double Tstress1;  // Stress without degradation and modification
    double Tstress2;  // Stress without modification but with degradation
    double Tstress3;  // stress with degradation and modification but without strength enhancement due to hardening
    double Tstress4;  // Stress with strength enhancement due to hardening
    double TCDD;  // Cumulative damage deformation
    bool Tfracture;  // Whether SMA cables fracture
    double Tplate1;  // Position of left end plate
    double Tplate2;  // Position of right end plate

    double Cstrain;
    double Ctangent;
    double Cstage;
    bool Chardening;
    double Cstress1;
    double Cstress2;
    double Cstress3;
    double Cstress4;
    double CCDD;
    bool Cfracture;
    double Cplate1;
    double Cplate2;

    void determineTrialState(double dStrain);
    double frictionModel(double F0, double du);
    double SCModel(double u0, double F0, double du);

};

#endif
