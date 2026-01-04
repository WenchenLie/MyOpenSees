import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class ModBoucWen(UniaxialMaterial):
    
    def __init__(self,
            tag: int,
            Fy: float,
            uy: float,
            alpha: float,
            n: float,
            Q: float,
            b: float,
            A: float,
            beta: float,
            gamma: float,
            iter: int=10
        ):
        self.tag = tag
        # Materail parameters
        self.Fy = Fy
        self.uy = uy
        self.alpha = alpha
        self.n = n
        self.Q = Q
        self.b = b
        self.A = A
        self.beta = beta
        self.gamma = gamma
        self.iter = iter
        # Response history
        self.Cstrain: float  # 上一步的应变
        self.Tstrain: float  # 当前步的应变
        self.Cstress: float  # 上一步的应力
        self.Tstress: float # 当前步的应力
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self.Cz: float  # 滞回变量z
        self.Tz: float
        self.Cwp: float  # 累积塑性应变
        self.Twp: float
        self.Cface: float  # 正向屈服面
        self.Tface: float
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        assert self.Fy > 0, "Fy must be positive"
        assert self.uy > 0, "uy must be positive"
        assert self.alpha >= 0, "alpha must not less than 0"
        assert self.n > 0, "n must be positive"

    def _init_paras(self):
        self.k0 = self.Fy / self.uy
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Ctangent = self.k0
        self.Ttangent = self.k0
        self.Cz = 0
        self.Tz = 0
        self.Cwp = 0
        self.Twp = 0
        self.Cface = self.uy
        self.Tface = self.uy

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = self.Cstrain
        self.Tstress = self.Cstress
        self.Ttangent = self.Ctangent
        self.Tz = self.Cz
        self.Twp = self.Cwp
        self.Tface = self.Cface
        dStrain = strain - self.Cstrain
        if dStrain == 0:
            return
        else:
            self.Tstrain = strain
        # Calculate stress and tangent
        dStrain_ = dStrain / self.iter
        z_ = self.Cz  # 临时的z值
        for i in range(self.iter):
            # 迭代iter个子步
            strain_ = self.Cstrain + dStrain_ * i  # 每一步应变
            if strain_ > self.Tface:
                # 正向屈服
                self.Twp += strain_ - self.Tface
                self.Tface = strain_
            elif strain_ < self.Tface - 2 * self.uy:
                # 负向屈服
                self.Twp += self.Tface - 2 * self.uy - strain_
                self.Tface = strain_ + 2 * self.uy
            if dStrain_ * z_ < 0.0:
                sgn = -1.0
            elif dStrain_ * z_ == 0.0:
                sgn = 0.0
            else:
                sgn = 1.0
            m = 1 + self.Q * (1 - pow(self.b, -self.Twp / self.uy))
            S1 = 1.0 / self.uy * (self.A - (self.beta * sgn + self.gamma) * pow(abs(z_ / m), self.n))
            S2 = 1.0 / self.uy * (self.A - (self.beta * sgn + self.gamma) * pow(abs(z_ / m + 0.5 * dStrain_ * S1), self.n))
            S3 = 1.0 / self.uy * (self.A - (self.beta * sgn + self.gamma) * pow(abs(z_ / m + 0.5 * dStrain_ * S2), self.n))
            S4 = 1.0 / self.uy * (self.A - (self.beta * sgn + self.gamma) * pow(abs(z_ / m + dStrain_ * S3), self.n))
            z_ = z_ + 1.0 / 6.0 * dStrain_ * (S1 + S2 + S3 + S4)
            self.Tstress = self.alpha * self.Fy / self.uy * strain_ + (1 - self.alpha) * self.Fy * z_
        else:
            # 将临时变量赋给实例变量
            self.Tz = z_
        self.Ttangent = (self.Tstress - self.Cstress) / dStrain

    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cz = self.Tz
        self.Cwp = self.Twp
        self.Cface = self.Tface

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent
    

