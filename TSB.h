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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TSB.h,v $
                                                                        
                                                                        
#ifndef TSB_h
#define TSB_h

// Written: Wenchen Lie 
// Created: July 26, 2024
// Revision: A
//
// Description: This file contains the class definition for 
// TSB.h
// 
//
//
// What: "@(#) TSB.h, revA"


#include <UniaxialMaterial.h>


class TSB : public UniaxialMaterial
{
  public:
    TSB(int tag, double Fslip, double k, double ugap, int N, const Vector &sc_args);
    TSB();
    ~TSB();

    const char *getClassType(void) const {return "TSB";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return k;};

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
    double Fslip;     // Friction slipping force
    double k;         // Initial stiffness
    double ugap;      // Gap length
    double N;         // Number of groups of parameters modeling self-centering behavior
    Vector sc_args;   // Parameters for self-centering behavior

    /*** State variables ***/
    double Tstrain;   // strain at current step
    double Tstress;   // stress at current step
    double Ttangent;  // tangent stiffness at current step
    double Cstrain;   // strain at previous step
    double Cstress;   // stress at previous step
    bool fracture;    // whether self-centering components are fracture
    int stage;        // stage (1 or 2)
    Vector F0_sc;     // trial stree of each self-centering component

    void determineTrialState(double dStrain);
    double frictionModel(double u0, double F0, double du);
    double SCmodel(double u0, double F0, double du, double Fy, double k1, double k2, double beta, double ubear, double kbear);

};

#endif
