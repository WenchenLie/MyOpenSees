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
// Created: February, 2011
// Revision: A
//
// Description: This file contains the class implementation for Maxwell Model

#include <GeneralizedMaxwell.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <OPS_Globals.h>

// 引入全局时间步以供积分使用
extern double ops_Dt;

void*
OPS_GeneralizedMaxwell()
{
    // 命令格式: GeneralizedMaxwell $matTag $k0 $k1 $c1 $alpha1 <$k2 $c2 $alpha2> ... <-iter $n_iter>
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 5) {
        opserr << "WARNING: Invalid #args, want: uniaxialMaterial GeneralizedMaxwell $tag $k0 $k1 $c1 $alpha1 ...\n";
        return 0;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
        opserr << "WARNING: invalid GeneralizedMaxwell tag\n";
        return 0;
    }

    double k0;
    if (OPS_GetDoubleInput(&numData, &k0) != 0) {
        opserr << "WARNING: invalid k0 for GeneralizedMaxwell tag " << tag << "\n";
        return 0;
    }

    std::vector<double> ks, cs, alphas;
    int n_iter = 10; // 默认迭代次数

    // 动态解析剩余的输入参数
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();

        if (strcmp(opt, "-iter") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                OPS_GetIntInput(&numData, &n_iter);
            }
            break; // 读取到 -iter 则认为分支参数读取完毕
        }
        else {
            // 如果不是 -iter 标签，则它必须是一个分层弹簧的 k 值
            double ki = atof(opt);

            // 还需要接着读入对应的 c 和 alpha
            if (OPS_GetNumRemainingInputArgs() < 2) {
                opserr << "WARNING: Missing c or alpha for branch in GeneralizedMaxwell tag " << tag << "\n";
                return 0;
            }

            double ci, alphai;
            if (OPS_GetDoubleInput(&numData, &ci) != 0) return 0;
            if (OPS_GetDoubleInput(&numData, &alphai) != 0) return 0;

            ks.push_back(ki);
            cs.push_back(ci);
            alphas.push_back(alphai);
        }
    }

    int n_layer = ks.size();
    if (n_layer == 0) {
        opserr << "WARNING: At least one Maxwell branch (k1, c1, alpha1) must be provided.\n";
        return 0;
    }

    // 转换为 OpenSees 的 Vector 对象
    Vector k_arr(n_layer), c_arr(n_layer), alpha_arr(n_layer);
    for (int i = 0; i < n_layer; i++) {
        k_arr(i) = ks[i];
        c_arr(i) = cs[i];
        alpha_arr(i) = alphas[i];
    }

    // 实例化材料对象，在外部库中通常使用0作为类标签(MAT_TAG)
    UniaxialMaterial* theMaterial = new GeneralizedMaxwell(tag, k0, n_layer, k_arr, c_arr, alpha_arr, n_iter);

    if (theMaterial == 0) {
        opserr << "WARNING: could not create uniaxialMaterial of type GeneralizedMaxwell\n";
        return 0;
    }

    return theMaterial;
}

GeneralizedMaxwell::GeneralizedMaxwell(int tag, double _k0, int _n_layer,
    const Vector& _k, const Vector& _c, const Vector& _alpha,
    int _n_iter)
    : UniaxialMaterial(tag, 0), k0(_k0), n_layer(_n_layer), n_iter(_n_iter),
    k_arr(_k), c_arr(_c), alpha_arr(_alpha),
    Cstress_i(_n_layer), Tstress_i(_n_layer)
{
    // 初始化状态变量
    Cstrain = 0.0;
    Tstrain = 0.0;
    Cstress = 0.0;
    Tstress = 0.0;

    // 初始化瞬态极限刚度 (玻璃态刚度)
    double k_sum = 0.0;
    for (int i = 0; i < n_layer; i++) {
        k_sum += k_arr(i);
        Cstress_i(i) = 0.0;
        Tstress_i(i) = 0.0;
    }

    Ctangent = k0 + k_sum;
    Ttangent = Ctangent;
}

GeneralizedMaxwell::GeneralizedMaxwell()
    : UniaxialMaterial(0, 0), k0(0.0), n_layer(0), n_iter(10),
    Cstress_i(0), Tstress_i(0)
{
    Cstrain = 0.0; Tstrain = 0.0;
    Cstress = 0.0; Tstress = 0.0;
    Ctangent = 0.0; Ttangent = 0.0;
}

GeneralizedMaxwell::~GeneralizedMaxwell()
{
    // Vector 类在销毁时会自动释放内存
}

// ======================= 核心的控制方程导数计算 =======================
void
GeneralizedMaxwell::compute_dS(const Vector& S, Vector& dS, double d_eps, double dt)
{
    for (int i = 0; i < n_layer; i++) {
        double S_val = S(i);
        double k_val = k_arr(i);
        double c_val = c_arr(i);
        double alpha_val = alpha_arr(i);

        double sign_S = 0.0;
        if (S_val > 0.0) sign_S = 1.0;
        else if (S_val < 0.0) sign_S = -1.0;

        double abs_S = fabs(S_val);
        // 阻尼器项，确保防除零报错 (alpha_val != 0 在构建前已默认由用户保证)
        double dashpot_vel = sign_S * pow(abs_S / c_val, 1.0 / alpha_val);

        dS(i) = k_val * d_eps - dt * k_val * dashpot_vel;
    }
}

int
GeneralizedMaxwell::setTrialStrain(double strain, double strainRate)
{
    Tstrain = strain;
    double d_eps = Tstrain - Cstrain;

    // 使用OpenSees全局时间步(对于动力时程)
    double dt = ops_Dt;

    Vector S = Cstress_i; // 拷贝上一步的收敛状态作为RK4起点
    double h = 1.0 / n_iter;

    Vector k1(n_layer), k2(n_layer), k3(n_layer), k4(n_layer);
    Vector S_tmp(n_layer);

    // 四阶龙格库塔(RK4)子步循环
    for (int iter = 0; iter < n_iter; iter++) {

        compute_dS(S, k1, d_eps, dt);

        for (int i = 0; i < n_layer; i++) S_tmp(i) = S(i) + 0.5 * h * k1(i);
        compute_dS(S_tmp, k2, d_eps, dt);

        for (int i = 0; i < n_layer; i++) S_tmp(i) = S(i) + 0.5 * h * k2(i);
        compute_dS(S_tmp, k3, d_eps, dt);

        for (int i = 0; i < n_layer; i++) S_tmp(i) = S(i) + h * k3(i);
        compute_dS(S_tmp, k4, d_eps, dt);

        for (int i = 0; i < n_layer; i++) {
            S(i) += (h / 6.0) * (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i));
        }
    }

    // 更新各分支与整体试验应力
    Tstress_i = S;
    Tstress = k0 * Tstrain;
    for (int i = 0; i < n_layer; i++) {
        Tstress += Tstress_i(i);
    }

    // 更新伪隐式算法切线刚度(Algorithmic Tangent)
    Ttangent = k0;
    for (int i = 0; i < n_layer; i++) {
        double Si = Tstress_i(i);
        double k_val = k_arr(i);
        double c_val = c_arr(i);
        double alpha_val = alpha_arr(i);

        double g = 0.0;
        if (fabs(Si) < 1.0e-14) { // 避免严格的 0.0
            if (alpha_val < 1.0) g = 0.0;
            else if (alpha_val == 1.0) g = k_val * dt / c_val;
            else g = 1.0e12; // 相当于无穷大
        }
        else {
            double power = 1.0 / alpha_val - 1.0;
            g = (k_val * dt / (alpha_val * c_val)) * pow(fabs(Si) / c_val, power);
        }

        if (g < 1.0e10) {
            Ttangent += k_val / (1.0 + g);
        }
    }

    return 0;
}

double GeneralizedMaxwell::getStress(void) { return Tstress; }
double GeneralizedMaxwell::getStrain(void) { return Tstrain; }
double GeneralizedMaxwell::getStrainRate(void) { return 0.0; }
double GeneralizedMaxwell::getTangent(void) { return Ttangent; }

double GeneralizedMaxwell::getInitialTangent(void)
{
    double K_init = k0;
    for (int i = 0; i < n_layer; i++) K_init += k_arr(i);
    return K_init;
}

int GeneralizedMaxwell::commitState(void)
{
    Cstrain = Tstrain;
    Cstress = Tstress;
    Ctangent = Ttangent;
    Cstress_i = Tstress_i; // 提交各层内力历史
    return 0;
}

int GeneralizedMaxwell::revertToLastCommit(void)
{
    Tstrain = Cstrain;
    Tstress = Cstress;
    Ttangent = Ctangent;
    Tstress_i = Cstress_i;
    return 0;
}

int GeneralizedMaxwell::revertToStart(void)
{
    Cstrain = 0.0; Tstrain = 0.0;
    Cstress = 0.0; Tstress = 0.0;

    double K_init = getInitialTangent();
    Ctangent = K_init;
    Ttangent = K_init;

    for (int i = 0; i < n_layer; i++) {
        Cstress_i(i) = 0.0;
        Tstress_i(i) = 0.0;
    }
    return 0;
}

UniaxialMaterial* GeneralizedMaxwell::getCopy(void)
{
    GeneralizedMaxwell* theCopy = new GeneralizedMaxwell(this->getTag(), k0, n_layer, k_arr, c_arr, alpha_arr, n_iter);

    theCopy->Cstrain = Cstrain;
    theCopy->Cstress = Cstress;
    theCopy->Ctangent = Ctangent;
    theCopy->Tstrain = Tstrain;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;
    theCopy->Cstress_i = Cstress_i;
    theCopy->Tstress_i = Tstress_i;

    return theCopy;
}

int GeneralizedMaxwell::sendSelf(int cTag, Channel& theChannel)
{
    int dataSize = 7 + 4 * n_layer;
    Vector data(dataSize);

    data(0) = this->getTag();
    data(1) = k0;
    data(2) = n_layer;
    data(3) = n_iter;
    data(4) = Cstrain;
    data(5) = Cstress;
    data(6) = Ctangent;

    int ptr = 7;
    for (int i = 0; i < n_layer; i++) data(ptr++) = k_arr(i);
    for (int i = 0; i < n_layer; i++) data(ptr++) = c_arr(i);
    for (int i = 0; i < n_layer; i++) data(ptr++) = alpha_arr(i);
    for (int i = 0; i < n_layer; i++) data(ptr++) = Cstress_i(i);

    int res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) opserr << "GeneralizedMaxwell::sendSelf() - failed to send data\n";
    return res;
}

int GeneralizedMaxwell::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    // 首先我们需要接收前几个参数，用于确定数组的大小
    Vector initData(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, initData) < 0) {
        opserr << "GeneralizedMaxwell::recvSelf() - failed to receive initData\n";
        return -1;
    }

    this->setTag((int)initData(0));
    k0 = initData(1);
    n_layer = (int)initData(2);
    n_iter = (int)initData(3);

    int dataSize = 7 + 4 * n_layer;
    Vector data(dataSize);

    // 重新获取包含全数据的Vector
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "GeneralizedMaxwell::recvSelf() - failed to receive data\n";
        return -1;
    }

    Cstrain = data(4);
    Cstress = data(5);
    Ctangent = data(6);

    k_arr.resize(n_layer);
    c_arr.resize(n_layer);
    alpha_arr.resize(n_layer);
    Cstress_i.resize(n_layer);
    Tstress_i.resize(n_layer);

    int ptr = 7;
    for (int i = 0; i < n_layer; i++) k_arr(i) = data(ptr++);
    for (int i = 0; i < n_layer; i++) c_arr(i) = data(ptr++);
    for (int i = 0; i < n_layer; i++) alpha_arr(i) = data(ptr++);
    for (int i = 0; i < n_layer; i++) Cstress_i(i) = data(ptr++);

    Tstrain = Cstrain;
    Tstress = Cstress;
    Ttangent = Ctangent;
    Tstress_i = Cstress_i;

    return 0;
}

void GeneralizedMaxwell::Print(OPS_Stream& s, int flag)
{
    s << "Generalized Maxwell Material, tag: " << this->getTag() << endln;
    s << "  k0: " << k0 << endln;
    s << "  Number of branches: " << n_layer << endln;
    for (int i = 0; i < n_layer; i++) {
        s << "  Branch " << i + 1 << " -> k: " << k_arr(i)
            << ", c: " << c_arr(i) << ", alpha: " << alpha_arr(i) << endln;
    }
}