import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class Steel01Material(UniaxialMaterial):
    
    def __init__(self, tag: int,
                 Fy: float,
                 k: float,
                 b: float
    ):
        self.tag = tag
        # Materail parameters
        self.Fy = Fy
        self.k = k
        self.b = b
        # Response history
        self.Cstrain: float  # 上一步的应变
        self.Tstrain: float  # 当前步的应变
        self.Cstress: float  # 上一步的应力
        self.Tstress: float # 当前步的应力
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        assert self.Fy > 0
        assert self.k > 0

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Ctangent = 0
        self.Ttangent = 0
        self.uy = self.Fy / self.k

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = strain
        dStrain = strain - self.Cstrain
        if abs(dStrain) > sys.float_info.epsilon:
            f = self.Tstress + dStrain * self.k
            if f > self.b * self.k * (self.Tstrain - self.uy) + self.Fy:
                f = self.b * self.k * (self.Tstrain - self.uy) + self.Fy
            elif f < self.b * self.k * (self.Tstrain + self.uy) - self.Fy:
                f = self.b * self.k * (self.Tstrain + self.uy) - self.Fy
            self.Tstress = f
            self.Ttangent = (self.Tstress - self.Cstress) / dStrain
        else:
            self.Tstress = self.Cstress
            self.Ttangent = self.Ctangent

    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent
    

