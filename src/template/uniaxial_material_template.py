import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class MyMaterial(UniaxialMaterial):
    
    def __init__(self, tag: int, *args, **kwargs):
        self.tag = tag
        # Materail parameters
        ...
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
        ...

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Ctangent = 0
        self.Ttangent = 0
        ...

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = strain
        self.Tstress = ...
        self.Ttangent = ...

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
    

