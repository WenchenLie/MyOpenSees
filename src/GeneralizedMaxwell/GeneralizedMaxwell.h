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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: Dimitrios G. Lignos, PhD, Assistant Professor, McGill University
// Created: February 2011
// Revision: A
//
// Description: This file contains the class interface for 
// Maxwell Model.  Maxwell defines a force(F)-velocity(v)
// relationship of the form F = K*u + C*pow(V,a)
// Variables:
// K: axial stiffness of a damper
// C: Velocity constant of a damper
// Alpha: Exponent of velocity of a damper
// L: Length of the Damper

// Olsson, A.K., and Austrell, P-E., (2001), "A fitting procedure for viscoelastic-elastoplastic material models," 
// Proceedings of the Second European Conference on Constitutive Models for Rubber, Germany, 2001.
// Ottosen, N.S., and Ristinmaa, M., (1999). "The mechanics of constitutive modelling, (Numerical and thermodynamical topics)," 
// Lund University,Division of Solid Mechanics, Sweden, 1999. 


#ifndef GeneralizedMaxwell_h
#define GeneralizedMaxwell_h

#include <UniaxialMaterial.h>
#include <Vector.h>

class GeneralizedMaxwell : public UniaxialMaterial
{
public:
    // 构造函数
    GeneralizedMaxwell(int tag, double k0, int n_layer,
        const Vector& k, const Vector& c, const Vector& alpha,
        int n_iter = 10);
    GeneralizedMaxwell();
    ~GeneralizedMaxwell();

    const char* getClassType(void) const { return "GeneralizedMaxwell"; };

    // 核心力学方法
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStrainRate(void);
    double getStress(void);

    double getTangent(void);
    double getInitialTangent(void);

    // 状态提交与回滚
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // 并行及数据流
    UniaxialMaterial* getCopy(void);
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);

    void Print(OPS_Stream& s, int flag = 0);

private:
    // 计算内部状态导数的辅助函数
    void compute_dS(const Vector& S, Vector& dS, double d_eps, double dt);

    // 材料参数
    double k0;
    int n_layer;
    int n_iter;

    Vector k_arr;
    Vector c_arr;
    Vector alpha_arr;

    // 历史变量与当前试验变量 (总应力应变)
    double Cstrain;
    double Cstress;
    double Ctangent;

    double Tstrain;
    double Tstress;
    double Ttangent;

    // 记录每个分支应力的向量
    Vector Cstress_i; // Committed stress in each layer
    Vector Tstress_i; // Trial stress in each layer
};

#endif

