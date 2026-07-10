import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class Degradation(UniaxialMaterial):
    def __init__(self,
        tag: int,
        other_tag: int,
        uy: float,
        u_th: float,
        a1: float=None,
        a2: float=None,
    ):
        self.tag = tag
        self.other_tag = other_tag
        # Materail parameters
        self.uy = uy
        self.u_th = u_th
        self.a1 = a1
        self.a2 = a2
        # Response history
        self.Cstate: bool  # 是否进入退化阶段
        self.Tstate: bool
        self.Cyieldface: float  # 屈服面
        self.Tyieldface: float
        self.Cuc: float  # 累积塑性位移
        self.Tuc: float
        self.Cuc2: float  # 进入退化状态后的累积塑形位移
        self.Tuc2: float
        self.Cumax: float  # 历史最大位移
        self.Tumax: float
        self.Cstrain: float  # 应变
        self.Tstrain: float
        self.Cstress: float
        self.Tstress: float
        self.Ctangent: float
        self.Ttangent: float
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        assert self.uy >= 0, "uy must be >= 0"
        assert self.u_th >= 0, "u_th must be >= 0"
        assert self.a1 >= 0, "a1 must be >= 0"
        assert self.a2 >= 0, "a2 must be >= 0"

    def _init_paras(self):
        self.Cstate = False
        self.Tstate = False
        self.Cyieldface = self.uy
        self.Tyieldface = self.uy
        self.Cuc = 0
        self.Tuc = 0
        self.Cuc2 = 0
        self.Tuc2 = 0
        self.Cumax = 0
        self.Tumax = 0
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Ctangent = 0
        self.Ttangent = 0
        self.material: UniaxialMaterial = UniaxialMaterial.getUniaxialMaterial(self.other_tag)

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        dStrain = strain - self.Cstrain
        if dStrain == 0:
            return
        else:
            self.Tstrain = strain
        self.material.setTrialStrain(strain, strainRate)
        self.Tumax = max(self.Tumax, abs(self.Tstrain))

        # 判断退化状态
        if self.Cstate:
            self.Tstate = True
        elif abs(self.Tstrain) >= self.u_th:
            self.Tstate = True  # 首次进入损伤状态
        else:
            self.Tstate = False
        
        # 更新屈服面，计算累积塑形位移
        if self.Tstrain > self.Cyieldface:
            # 正向屈服
            self.Tuc = self.Cuc + self.Tstrain - self.Cyieldface
            self.Tyieldface = self.Tstrain
            if self.Tstate:
                self.Tuc2 = self.Cuc2 + self.Tstrain - self.Cyieldface
        elif self.Tstrain < self.Cyieldface - 2 * self.uy:
            # 负向屈服
            self.Tuc = self.Cuc + self.Cyieldface - 2 * self.uy - self.Tstrain
            self.Tyieldface = self.Tstrain + 2 * self.uy
            if self.Tstate:
                self.Tuc2 = self.Cuc2 + self.Cyieldface - 2 * self.uy - self.Tstrain

        # 计算循环退化
        if self.Tstate:
            lamda1: float = 1 - self.a1 * self.Tuc2 / self.u_th
            lamda1 = max(lamda1, 0)
        else:
            lamda1: float = 1
        
        # 计算最大位移退化
        if self.Tstate:
            lamda2: float = 1 - self.a2 * (self.Tuc - self.Tuc2) / self.Tuc2
            lamda2 = max(lamda2, 0)
        else:
            lamda2: float = 1
        
        # 修正力响应
        lamda: float = min(lamda1, lamda2)
        self.material.setTrialStrain(self.Tstrain, strainRate)
        self.Tstress = self.material.getStress() * lamda
        self.Ttangent = self.material.getTangent() * lamda

    def commitState(self):
        self.Cstate = self.Tstate
        self.Cuc = self.Tuc
        self.Cuc2 = self.Tuc2
        self.Cumax = self.Tumax
        self.Cyieldface = self.Tyieldface
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ttangent = self.Ttangent
        self.material.commitState()

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent

