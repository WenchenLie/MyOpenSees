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
// $Date: 2024-07-27 11:58:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ModBoucWen.h,v $
                                                                        
                                                                        
#ifndef ModBoucWen_h
#define ModBoucWen_h

// Written: Wenchen Lie 
// Created: July 26, 2024
// Revision: A
//
// Description: This file contains the class definition for 
// ModBoucWen.h
// 
//
//
// What: "@(#) ModBoucWen.h, revA"


#include <UniaxialMaterial.h>


class ModBoucWen : public UniaxialMaterial
{
  public:
    ModBoucWen(int tag, double Fy, double uy, double alpha, double n,
        double Q, double b, double A, double beta, double gamma,
        int gotFailureCPD, double failureCPD, int gotMinMax, double MinMax);
    ModBoucWen();
    ~ModBoucWen();

    const char *getClassType(void) const {return "ModBoucWen";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return Fy/ uy;};

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
    double Fy;
    double uy;
    double alpha;
    double n;
    double Q;
    double b;
    double A;
    double beta;
    double gamma;
    int gotFailureCPD;
    int gotMinMax;
    double failureCPD;
    double MinMax;

    /*** State variables ***/
    double Tstrain;   // strain at current step
    double Tstress;   // stress at current step
    double Ttangent;  // tangent stiffness at current step
    double Cstrain;   // strain at previous step
    double Cstress;   // stress at previous step
    double Tz;     // varialbe z
    double Twp;    // cumulative plastic displacement
    double Tface;  // yielding face
    double Cz;
    double Cwp;
    double Cface;
    bool failure;  // whether PCD exceed limit

    void determineTrialState(double dStrain);

};

#endif
