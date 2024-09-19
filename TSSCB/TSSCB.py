import sys
from typing import Literal
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
        hardening: Literal['-hardening', None]=None,
        uh: float=1e16,
        r1: float=1,
        r2: float=1,
        r3: float=0,
        minmax: Literal['-minmax', None]=None,
        uf: float=1e16
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
            hardening (Literal['-hardening', None]): 是否考虑硬化
            uh (float): 硬化开始位移
            r1 (float): 退化系数，控制第二阶段开始时的承载力退化程度
            r2 (float): 退化系数，控制第二阶段结束时的承载力退化程度
            r3 (float): 硬化系数，控制硬化后的硬化刚度
            minmax (Literal['-minmax', None]): 是否SMA断裂
            uf (float): SMA线缆断裂时的位移
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
        if hardening is not None:
            if hardening != '-hardening':
                raise ValueError('`hardening` should be `-hardening` if it is given')
            self.hasHardening = True  # Whether to consider hardening
        else:
            self.hasHardening = False
        if minmax is not None:
            if minmax != '-minmax':
                raise ValueError('`minmax` should be `-minmax` if it is given')
        self.uh = uh
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.uf = uf
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
        self.Chardening: bool  # 是否开始硬化
        self.Thardening: bool
        self.Cstress_ideal: float  # 无退化，无修正应力
        self.Tstress_ideal: float
        self.Cstress_ideal1: float  # 有退化，无修正应力
        self.Tstress_ideal1: float
        self.CCDD: float
        self.TCDD: float
        self.Cfracture: bool  # 是否断裂
        self.Tfracture: bool
        self.Cplate1: float  # 左端板位置
        self.Tplate1: float
        self.Cplate2: float  # 右端板位置
        self.Tplate2: float
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
            self.Ctangent = self.k1
            self.Ttangent = self.k1
            self.Cstage = 2
            self.Tstage = 2
        else:
            self.Ctangent = self.k0
            self.Ttangent = self.k0
        self.Chardening = False
        self.Thardening = False
        self.Cstress_ideal = 0
        self.Tstress_ideal = 0
        self.Cstress_ideal1 = 0
        self.Tstress_ideal1 = 0
        self.CCDD = 0  # Dimensionless cumulative damage deformation
        self.TCDD = 0
        self.Cfracture = False  # Whether the cable is fractured
        self.Tfracture = False
        self.Cplate1 = self.ugap  # Position of left end plate
        self.Tplate1 = self.ugap
        self.Cplate2 = -self.ugap  # Position of right end plate
        self.Tplate2 = -self.ugap
        self.i = 0  # For test


    def _check_paras(self):
        assert self.F1 >= 0
        assert self.k0 > 0
        assert self.ugap >= 0
        assert self.F2 > 0
        assert self.k1 > 0
        assert self.k2 > 0
        assert 0 <= self.beta <= 1
        assert self.uh > 0
        assert self.r1 >= 0
        assert self.r2 >= 0
        assert self.r3 >= 0
        assert self.uf > 0


    def _setTrainStrain(self, strain, strainRate: float=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = self.Cstrain
        self.Tstress = self.Cstress
        self.Ttangent = self.Ctangent
        self.Tstage = self.Cstage
        self.Thardening = self.Chardening
        self.Tstress_ideal = self.Cstress_ideal
        self.Tstress_ideal1 = self.Cstress_ideal1
        self.TCDD = self.CCDD
        self.Tfracture = self.Cfracture
        self.Tplate1 = self.Cplate1
        self.Tplate2 = self.Cplate2
        dStrain = strain - self.Cstrain
        if abs(dStrain) > sys.float_info.epsilon:
            self.Tstrain = strain
            # Determine whether to start hardening
            if abs(self.Tstrain) > self.uh:
                self.Thardening = True
            if self.Thardening and abs(self.Tstrain) > self.ugap:
                self.TCDD = self.CCDD + abs(dStrain) / self.uh
            if abs(self.Tstrain) > self.uf:
                self.Tfracture = True
            # Update Tstress
            self._determineTrialState(dStrain)
            # Update endplate position
            if dStrain > 0:
                self.Tplate1 = max(self.Tplate1, self.Tstrain)
                if not self.Tfracture:
                    self.Tplate2 += dStrain
            else:
                self.Tplate2 = min(self.Tplate2, self.Tstrain)
                if not self.Tfracture:
                    self.Tplate1 += dStrain
            if self.Tplate1 < self.ugap:
                self.Tplate1 = self.ugap
            if self.Tplate2 > -self.ugap:
                self.Tplate2 = -self.ugap
            # Calculate tangent stiffness
            self.Ttangent = (self.Tstress - self.Cstress) / dStrain
        else:
            pass


    def _determineTrialState(self, dStrain: float):
        """Determine the trial state based on the current strain increment"""
        # Determine stage
        if -self.ugap <= self.Tstrain <= self.ugap:
            self.Tstage = 1  # stage-1
        else:
            self.Tstage = 2  # stage-2
        if self.ugap == 0:
            self.Tstage = 2  # If ugap is zero, always in stage-2
        if self.Tfracture:
            if self.Tplate2 <= self.Tstrain <= self.Tplate1:
                self.Tstress = 0
            else:
                if dStrain > 0 and self.Tstrain > 0:
                    self.Tstress = self.F1
                elif dStrain < 0 and self.Tstrain > 0:
                    self.Tstress = 0
                elif dStrain > 0 and self.Tstrain < 0:
                    self.Tstress = 0
                else:
                    self.Tstress = -self.F1
            return
        if self.Cstage == 1 and self.Tstage == 1:
            # NOTE: stage-1 -> stage-1
            self.Tstress_ideal = self._frictionModel(self.Cstress, dStrain)
            self.Tstress_ideal1 = self.Tstress_ideal
            self.Tstress = self.Tstress_ideal
        elif self.Cstage == 1 and self.Tstage == 2:
            # NOTE: stage-1 -> stage-2
            if dStrain > 0:
                du1 = self.ugap - self.Cstrain  # Strain increment in stage-1
                du2 = dStrain - du1  # Strain increment in stage-2
                usc0 = self.ugap - self.ua
            else:
                du1 = -self.ugap - self.Cstrain  # Strain increment in stage-1
                du2 = dStrain - du1  # Strain increment in stage-2
                usc0 = self.ua - self.ugap
            F1_ = self._frictionModel(self.Cstress_ideal, du1)
            F2_ = self._SCModel(usc0, F1_, du2)
            self.Tstress_ideal = F2_
            # Apply degradation
            self.Tstress_ideal1 = self.Tstress_ideal
            if self.Thardening and self.Tstrain > 0:
                self.Tstress_ideal1 = self.Tstress_ideal - self.F2 * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / self.uh)
            elif self.Thardening and self.Tstrain < 0:
                self.Tstress_ideal1 = self.Tstress_ideal + self.F2 * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / self.uh)
            # Apply modifiction
            self.Tstress = self.Tstress_ideal1
            if dStrain > 0 and self.Tstress_ideal1 < self.F1:
                self.Tstress = self.F1
            elif dStrain < 0 and self.Tstress_ideal1 > -self.F1:
                self.Tstress = -self.F1
        elif self.Cstage == 2 and self.Tstage == 2:
            # NOTE: stage-2 -> stage-2
            if self.Tstrain >= 0:
                usc0 = self.Cstrain - self.ua
            else:
                usc0 = self.Cstrain + self.ua
            self.Tstress_ideal = self._SCModel(usc0, self.Cstress_ideal, dStrain)
            # Apply degradation
            self.Tstress_ideal1 = self.Tstress_ideal
            if self.Thardening and self.Tstrain > 0:
                self.Tstress_ideal1 = self.Tstress_ideal - self.F2 * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / self.uh)
            elif self.Thardening and self.Tstrain < 0:
                self.Tstress_ideal1 = self.Tstress_ideal + self.F2 * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / self.uh)
            # Apply modifiction
            self.Tstress = self.Tstress_ideal1
            if dStrain > 0 and self.Tstrain > 0 and self.Tstress_ideal1 < self.F1 and self.ugap > 0 and self.Cstress == self.F1:
                self.Tstress = self.F1
            elif dStrain < 0 and self.Tstrain < 0 and self.Tstress_ideal1 > -self.F1 and self.ugap > 0 and self.Cstress == -self.F1:
                self.Tstress = -self.F1
            elif self.Tstrain > 0 and self.Tstress_ideal1 < 0:
                self.Tstress = 0  # Prevent positive compressive stress in SMA cables
            elif self.Tstrain < 0 and self.Tstress_ideal1 > 0:
                self.Tstress = 0  # Prevent negative compressive stress in SMA cables
        elif self.Cstage == 2 and self.Tstage == 1:
            # NOTE: stage-2 -> stage-1
            if dStrain < 0:
                du1 = -(self.Cstrain - self.ugap)
                du2 = -(self.ugap - self.Tstrain)
                usc0 = self.Cstrain - self.ua
            else:
                du1 = -self.ugap - self.Cstrain
                du2 = self.Tstrain + self.ugap
                usc0 = self.Cstrain + self.ua
            F1_ = self._SCModel(usc0, self.Cstress_ideal, du1)
            # Apply degradation
            F1_ideal1 = F1_
            if self.Thardening and self.Tstrain > 0:
                F1_ideal1 = F1_ - self.F2 * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / self.uh)
            elif self.Thardening and self.Tstrain < 0:
                F1_ideal1 = F1_ + self.F2 * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / self.uh)
            # Apply modifiction
            F1_ = F1_ideal1
            if dStrain > 0 and self.Tstrain > 0 and F1_ideal1 < self.F1 and self.ugap > 0 and self.Cstress == self.F1:
                F1_ = self.F1
            elif dStrain < 0 and self.Tstrain < 0 and F1_ideal1 > -self.F1 and self.ugap > 0 and self.Cstress == -self.F1:
                F1_ = -self.F1
            elif self.Tstrain > 0 and F1_ideal1 < 0:
                F1_ = 0  # Prevent positive compressive stress in SMA cables
            elif self.Tstrain < 0 and F1_ideal1 > 0:
                F1_ = 0  # Prevent negative compressive stress in SMA cables
            F2_ = self._frictionModel(F1_, du2)
            self.Tstress_ideal = F2_
            self.Tstress_ideal1 = self.Tstress_ideal
            self.Tstress = self.Tstress_ideal
            # if dStrain < 0 and self.Tstress < 0:
            #     self.Tstress = 0
            # elif dStrain > 0 and self.Tstress > 0:
            #     self.Tstress = 0
        else:
            assert False, f"Invalid state transition, Cstage={self.Cstage}, Tstage={self.Tstage}"
        F_hardening = max(abs(self.Tstrain) - self.uh, 0) * self.k2 * self.r3  # Strength enhancement due to hardening
        if self.Tstrain > 0:
            self.Tstress += F_hardening
        else:
            self.Tstress -= F_hardening
        # self.Tstress = self.Tstress_ideal  # HACK: Test for ideal stress
    

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
        if self.Tfracture:
            return 0
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
        self.Chardening = self.Thardening
        self.Cstress_ideal = self.Tstress_ideal
        self.Cstress_ideal1 = self.Tstress_ideal1
        self.CCDD = self.TCDD
        self.Cfracture = self.Tfracture
        self.Cplate1 = self.Tplate1
        self.Cplate2 = self.Tplate2
        self.i += 1


    def getStrain(self):
        return self.Tstrain


    def getStress(self):
        return self.Tstress


    def getTangent(self):
        return self.Ttangent
    

