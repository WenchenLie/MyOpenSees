import sys
from utils.uniaxialmaterial import UniaxialMaterial


class TSSCBMaterial(UniaxialMaterial):
    
    def __init__(self,
        tag: int,
        F1: float,
        k0: float,
        ugap: float,
        F2: float,
        k1: float,
        k2: float,
        beta: float,
    ):
        """双阶自复位本构

        Args:
            tag (int): 独立的材料编号
            F1 (float): 滑动摩擦力
            k0 (float): 摩擦初始刚度
            ugap (float): 间隙
            F2 (float): 第二阶段启动力
            k1 (float): 第二阶段的第一刚度
            k2 (float): 第二阶段的第二刚度
            beta (float): 耗能系数
        """
        self.tag = tag
        # Materail parameters
        self.F1 = F1
        self.k0 = k0
        self.ugap = ugap
        self.F2 = F2
        self.k1 = k1
        self.k2 = k2
        self.beta = beta
        # Calculated parameters
        self.ua = self.ugap - self.F1 / self.k1
        self.ua = max(0, self.ua)  # 第一阶段进入第二阶段时自复位分量的初始应变
        # Response history
        self.Cstage: int  # 上一步所处阶段
        self.Tstage: int  # 当前步所处阶段
        self.Cstrain: float  # 上一步的应变
        self.Tstrain: float  # 当前步的应变
        self.Cstress: float  # 上一步的应力
        self.Tstress: float # 当前步的应力
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self._check_paras()
        self._init_paras()

    
    def _init_paras(self):
        self.Cstage = 1
        self.Tstage = 1
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        if self.ugap == 0:
            self.Ctangent = self.k0
            self.Ttangent = self.k0
        else:
            self.Ctangent = self.k1
            self.Ttangent = self.k1


    def _check_paras(self):
        assert self.F1 > 0
        assert self.k0 > 0
        assert self.ugap >= 0
        assert self.F2 > 0
        assert self.k1 > 0
        assert self.k2 > 0
        assert 0 <= self.beta <= 1


    def _setTrainStrain(self, strain, strainRate: float=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = self.Cstrain
        self.Tstress = self.Cstress
        self.Ttangent = self.Ctangent
        self.Tstage = self.Cstage
        dStrain = strain - self.Cstrain
        if abs(dStrain) > sys.float_info.epsilon:
            self.Tstrain = strain
            self._determineTrialState(dStrain)
        else:
            pass


    def _determineTrialState(self, dStrain: float):
        """Determine the trial state based on the current strain increment"""
        # Determine stage
        if -self.ugap <= self.Tstrain <= self.ugap:
            self.Tstage = 1  # stage-1
        else:
            self.Tstage = 2  # stage-2
        # Determine trial stress
        if self.Cstage == 1 and self.Tstage == 1:
            # stage-1 -> stage-1
            self.Tstress = self._frictionModel(self.Cstress, dStrain)
        elif self.Cstage == 1 and self.Tstage == 2:
            # stage-1 -> stage-2
            if dStrain > 0:
                du1 = self.ugap - self.Cstrain  # Strain increment in stage-1
                du2 = dStrain - du1  # Strain increment in stage-2
                usc0 = self.ugap - self.ua
            else:
                du1 = -self.ugap - self.Cstrain  # Strain increment in stage-1
                du2 = dStrain - du1  # Strain increment in stage-2
                usc0 = self.ua - self.ugap
            F1_ = self._frictionModel(self.Cstress, du1)
            F2_ = self._SCModel(usc0, F1_, du2)
            self.Tstress = F2_
        elif self.Cstage == 2 and self.Tstage == 2:
            # stage-2 -> stage-2
            if self.Cstrain >= 0:
                usc0 = self.Cstrain - self.ua
            else:
                usc0 = self.Cstrain + self.ua
            self.Tstress = self._SCModel(usc0, self.Cstress, dStrain)
        elif self.Cstage == 2 and self.Tstage == 1:
            # stage-2 -> stage-1
            if dStrain < 0:
                du1 = -(self.Cstrain - self.ugap)  # Strain increment in stage-2
                du2 = -(self.ugap - self.Tstrain)  # Strain increment in stage-1
                usc0 = self.Cstrain - self.ua
            else:
                du1 = -self.ugap - self.Cstrain  # Strain increment in stage-2
                du2 = self.Tstrain + self.ugap  # Strain increment in stage-1
                usc0 = self.Cstrain + self.ua
            F1_ = self._SCModel(usc0, self.Cstress, du1)
            F2_ = self._frictionModel(F1_, du2)
            self.Tstress = F2_
        else:
            assert False, f"Invalid state transition, Cstage={self.Cstage}, Tstage={self.Tstage}"
        self.Ttangent = (self.Tstress - self.Cstress) / dStrain
    

    def _frictionModel(self,
            F0: float,
            du: float,
        ) -> float:
        """滑动摩擦力模型

        Args:
            F0 (float): 上一步摩擦力
            du (float): 位移增量

        Returns:
            float: 当前步摩擦力
        """
        if du == 0:
            return F0
        F_ = F0 + du * self.k0
        if F_ > self.F1:
            F = self.F1
        elif F_ < -self.F1:
            F = -self.F1
        else:
            F = F_
        return F


    def _SCModel(self,
            u0: float,
            F0: float,
            du: float,
        ) -> float:
        """旗帜型自复位模型

        Args:
            u0 (float): 上一步位移(无滑移位移)
            F0 (float): 上一步力
            du (float): 位移增量

        Returns:
            float: 当前步力
        """
        if du == 0:
            F = F0
            return F
        u = u0 + du
        uy = self.F2 / self.k1
        F_ = F0 + du * self.k1
        if du > 0:
            if u < -self.F2 * (1 - 2 * self.beta) / self.k1 and F_ > self.k2 * u - self.F2 * (1 - 2 * self.beta) * (1 - self.k2 / self.k1):
                F = self.k2 * u - self.F2 * (1 - 2 * self.beta) * (1 - self.k2 / self.k1)
            elif -self.F2 * (1 - 2 * self.beta) / self.k1 <= u <= uy and F_ > self.k1 * u:
                F = self.k1 * u
            elif u > self.F2 * (1 - 2 * self.beta) / self.k1 and F_ > self.k2 * u + self.F2 - self.k2 * uy:
                F = self.k2 * u + self.F2 - self.k2 * uy
            else:
                F = F_
        else:
            if u > self.F2 * (1 - 2 * self.beta) / self.k1 and F_ < self.k2 * u + self.F2 * (1 - 2 * self.beta) * (1 - self.k2 / self.k1):
                F = self.k2 * u + self.F2 * (1 - 2 * self.beta) * (1 - self.k2 / self.k1)
            elif -uy <= u <= self.F2 * (1 - 2 * self.beta) / self.k1 and F_ < self.k1 * u:
                F = self.k1 * u
            elif u < -self.F2 * (1 - 2 * self.beta) / self.k1 and F_ < self.k2 * u - (self.F2 - self.k2 * uy):
                F = self.k2 * u - (self.F2 - self.k2 * uy)
            else:
                F = F_
        return F

    def _commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cstage = self.Tstage

    def getStrain(self):
        return self.Tstrain


    def getStress(self):
        return self.Tstress


    def getTangent(self):
        return self.Ttangent
    

