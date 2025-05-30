import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class Failure(UniaxialMaterial):
    
    def __init__(self,
            tag: int,
            other_tag: int,
            minStrain: float=None,
            maxStrain: float=None,
            minForce: float=None,
            maxForce: float=None,
            uy: float=None,
            maxCPD: float=None,
        ):
        self.tag = tag
        self.other_tag = other_tag
        # Materail parameters
        self.minStrain = minStrain
        self.maxStrain = maxStrain
        self.minForce = minForce
        self.maxForce = maxForce
        self.maxCPD = maxCPD
        self.uy = uy
        # Response history
        self.Cfailure = False  # 是否破坏
        self.Tfailure = False
        self.Cstrain: float  # 上一步的应变
        self.Tstrain: float  # 当前步的应变
        self._check_paras()
        self._init_paras()
        # other material

    def _check_paras(self):
        if self.minStrain is None:
            self.minStrain = -sys.float_info.max
        if self.maxStrain is None:
            self.maxStrain = sys.float_info.max
        if self.minForce is None:
            self.minForce = -sys.float_info.max
        if self.maxForce is None:
            self.maxForce = sys.float_info.max
        if self.maxCPD is not None and self.uy is None:
            raise ValueError("uy must be provided when maxCPD is provided")
        if self.maxCPD is None:
            self.maxCPD = sys.float_info.max
        if self.uy is None:
            self.uy = sys.float_info.max
        assert self.minStrain <= 0, "minStrain must be negative"
        assert self.maxStrain >= 0, "maxStrain must be positive"
        assert self.minForce <= 0, "minForce must be negative"
        assert self.maxForce >= 0, "maxForce must be positive"
        assert self.maxCPD > 0, "maxCPD must be positive"
        assert self.uy > 0, "uy must be positive"

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cfailure = False
        self.Tfailure = False
        self.Cyieldface = self.uy
        self.Tyieldface = self.uy
        self.Cwp = 0
        self.Twp = 0
        self.material: UniaxialMaterial = UniaxialMaterial.getUniaxialMaterial(self.other_tag)

    def setTrainStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = self.Cstrain
        self.Tfailure = self.Cfailure
        self.Tyieldface = self.Cyieldface
        self.Twp = self.Cwp
        dStrain = strain - self.Cstrain
        if dStrain == 0:
            return
        else:
            self.Tstrain = strain
        self.material._setTrainStrain(strain, strainRate)
        if self.Tstrain > self.Cyieldface:
            # 正向屈服
            self.Twp = self.Cwp + self.Tstrain - self.Cyieldface
            self.Tyieldface = self.Tstrain
        elif self.Tstrain < self.Cyieldface - 2 * self.uy:
            # 负向屈服
            self.Twp = self.Cwp + self.Cyieldface - 2 * self.uy - self.Tstrain
            self.Tyieldface = self.Tstrain + 2 * self.uy
        # 1 判断是否延性破坏
        if self.Tstrain < self.minStrain:
            # 负向破坏
            self.Tfailure = True
        elif self.Tstrain > self.maxStrain:
            # 正向破坏
            self.Tfailure = True
        # 2 判断是否累积塑性应变破坏
        if self.Twp > self.maxCPD * self.uy:
            self.Tfailure = True
        # 3 判断是否承载力破坏
        if self.material.getStress() < self.minForce:
            self.Tfailure = True
        elif self.material.getStress() > self.maxForce:
            self.Tfailure = True

    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cfailure = self.Tfailure
        self.Cyieldface = self.Tyieldface
        self.Cwp = self.Twp
        self.material.commitState()

    def getStrain(self):
        return self.material.getStrain()

    def getStress(self):
        if self.Tfailure:
            return 0
        else:
            return self.material.getStress()

    def getTangent(self):
        if self.Tfailure:
            return 0
        else:
            return self.material.getTangent()

