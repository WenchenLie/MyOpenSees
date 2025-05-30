import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class TwoStage(UniaxialMaterial):
    
    def __init__(self,
            tag: int,
            F1: float,
            k1: float,
            kp1: float,
            F2: float,
            k2: float,
            kp2: float,
            ua: float
        ):
        self.tag = tag
        # Materail parameters
        self.F1 = F1
        self.k1 = k1
        self.kp1 = kp1
        self.F2 = F2
        self.k2 = k2
        self.kp2 = kp2
        self.ua = ua
        # Response history
        self.Cstrain: float  # 应变
        self.Tstrain: float
        self.Cstress: float  # 应力
        self.Tstress: float
        self.Ctangent: float  # 切线刚度
        self.Ttangent: float
        self.Cstage: int  # 阶段(1 or 2)
        self.Tstage: int
        self.Cu_pos: float  # 正向转换面位移
        self.Tu_pos: float
        self.Cu_neg: float  # 负向转换面位移
        self.Tu_neg: float
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        assert self.F1 > 0, "Fy must be positive"
        assert self.k1 > 0, "k1 must be positive"
        assert self.kp1 >= 0, "kp1 must be non-negative"
        assert self.F2 > 0, "F2 must be positive"
        assert self.k2 > 0, "k2 must be positive"
        assert self.kp2 >= 0, "kp2 must be non-negative"
        assert self.ua >= 0, "ua must be non-negative"

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Ctangent = self.k1
        self.Ttangent = self.k1
        self.Cstage = 1
        self.Tstage = 1
        self.Cu_pos = self.ua
        self.Tu_pos = self.ua
        self.Cu_neg = -self.ua
        self.Tu_neg = -self.ua

    def setTrainStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        dStrain = strain - self.Cstrain
        self.Tstrain = strain
        if abs(dStrain) <= sys.float_info.epsilon:
            return
        # Determine stage
        if self.Tstrain > self.Cu_pos:
            self.Tstage = 2
            self.Tu_pos = self.Tstrain
            self.Tu_neg = self.Tu_pos - 2 * self.ua
        elif self.Tstrain < self.Cu_neg:
            self.Tstage = 2
            self.Tu_neg = self.Tstrain
            self.Tu_pos = self.Tu_neg + 2 * self.ua
        else:
            self.Tstage = 1
        # Calculate stress and tangent
        if self.Cstage == 1 and self.Tstage == 1:
            # Remain stage-1
            self.Tstress = self.bilinear(self.Cstress, self.Cstrain, dStrain, self.F1, self.k1, self.kp1)
        elif self.Cstage == 2 and self.Tstage == 2:
            # Remain stage-2
            self.Tstress = self.bilinear(self.Cstress, self.Cstrain, dStrain, self.F2, self.k2, self.kp2)
        elif self.Cstage == 1 and self.Tstage == 2:
            # Transition stage-1 to stage-2
            du1: float  # Strain component in stage-1
            du2: float  # Strain component in stage-2
            F_trans: float
            if dStrain > 0:
                du2 = self.Tstrain - self.Cu_pos
            else:
                du2 = self.Tstrain - self.Cu_neg
            du1 = dStrain - du2
            F_trans = self.bilinear(self.Cstress, self.Cstrain, du1, self.F1, self.k1, self.kp1)
            self.Tstress = self.bilinear(F_trans, self.Cstrain + du1, du2, self.F2, self.k2, self.kp2)
        elif self.Cstage == 2 and self.Tstage == 1:
            # Transition stage-2 to stage-1
            du2: float  # Strain component in stage-2
            du1: float  # Strain component in stage-1
            F_trans: float
            if dStrain < 0:
                du1 = self.Tstrain - self.Cu_pos
            else:
                du1 = self.Tstrain - self.Cu_neg
            du2 = dStrain - du1
            F_trans = self.bilinear(self.Cstress, self.Cstrain, du2, self.F2, self.k2, self.kp2)
            self.Tstress = self.bilinear(F_trans, self.Cstrain + du2, du1, self.F1, self.k1, self.kp1)
        else:
            assert False, "Invalid stage"
        self.Ttangent = (self.Tstress - self.Cstress) / dStrain

    @staticmethod
    def bilinear(
            F_prev: float,
            u_prev: float,
            du: float,
            Fy: float,
            k: float,
            kp: float,
        ) -> float:
        """双线性模型

        Args:
            F_prev (float): 上一步应力
            u_prev (float): 上一步应变
            du (float): 位移增量
            Fy (float): 屈服力
            k (float): 弹性刚度
            kp (float): 塑性模量

        Returns:
            float: 当前步应力
        """
        F_next = F_prev + du * k
        if F_next > kp * (u_prev + du) + (1 - kp / k) * Fy and du > 0:
            F_next = kp * (u_prev + du) + (1 - kp / k) * Fy
        elif F_next < kp * (u_prev + du) - (1 - kp / k) * Fy and du < 0:
            F_next = kp * (u_prev + du) - (1 - kp / k) * Fy
        return F_next
    
    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cu_pos = self.Tu_pos
        self.Cu_neg = self.Tu_neg
        self.Cstage = self.Tstage

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent
    

