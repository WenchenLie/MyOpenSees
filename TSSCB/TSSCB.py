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
        uf: float=1e16,
        configType: Literal['-configType', None]=None,
        configTypeVal: Literal[1, 2]=1
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
            configType (Literal['-configType', None]): 退化类型，默认1
            configTypeVal (Literal[1, 2]): 退化类型, 1: 第二阶段只有一半摩擦片滑动, 2: 第二阶段所有摩擦片均滑动
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
        if configType is not None:
            if configType != '-configType':
                raise ValueError('`configType` should be `-configType` if it is given')
        self.uh = uh
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.uf = uf
        self.configTypeVal = configTypeVal
        # Calculated parameters
        self.ua = self.ugap - self.F1 / self.k1
        self.ua = max(0, self.ua)  # 第一阶段进入第二阶段时自复位分量的初始应变
        # Response history
        self.Cstage: int  # 上一步所处阶段
        self.Tstage: int  # 当前步所处阶段
        self.Cstress: float  # 上一步的应力
        self.Tstress: float # 当前步的应力
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self.Chardening: bool  # 是否开始硬化
        self.Thardening: bool
        self.Cstress1: float  # 无退化，无修正应力
        self.Tstress1: float
        self.Cstress2: float  # 有退化，无修正应力
        self.Tstress2: float
        self.Cstrain3: float  # 考虑退化及修正后应力，无硬化强度提升
        self.Tstrain3: float
        self.Tstress4: float  # 考虑硬化提升的强度(最终返回的应力)
        self.Cstress4: float
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
        self.Cstrain3 = 0
        self.Tstrain3 = 0
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
        self.Cstress1 = 0
        self.Tstress1 = 0
        self.Cstress2 = 0
        self.Tstress2 = 0
        self.Cstress4 = 0
        self.Tstress4 = 0
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
        self.Tstrain3 = self.Cstrain3
        self.Tstress = self.Cstress
        self.Ttangent = self.Ctangent
        self.Tstage = self.Cstage
        self.Thardening = self.Chardening
        self.Tstress1 = self.Cstress1
        self.Tstress2 = self.Cstress2
        self.Tstress4 = self.Cstress4
        self.TCDD = self.CCDD
        self.Tfracture = self.Cfracture
        self.Tplate1 = self.Cplate1
        self.Tplate2 = self.Cplate2
        dStrain = strain - self.Cstrain3
        if abs(dStrain) > sys.float_info.epsilon:
            self.Tstrain3 = strain
            # Determine whether to start hardening
            if abs(self.Tstrain3) > self.uh or self.Chardening:
                self.Thardening = True
            if abs(self.Tstrain3) > self.uf:
                self.Tfracture = True
            # Update Tstress
            self._determineTrialState(dStrain)
            # Update endplate position
            if dStrain > 0:
                self.Tplate1 = max(self.Tplate1, self.Tstrain3)
                if not self.Tfracture:
                    self.Tplate2 += dStrain
            else:
                self.Tplate2 = min(self.Tplate2, self.Tstrain3)
                if not self.Tfracture:
                    self.Tplate1 += dStrain
            if self.Tplate1 < self.ugap:
                self.Tplate1 = self.ugap
            if self.Tplate2 > -self.ugap:
                self.Tplate2 = -self.ugap
            # Calculate tangent stiffness
            self.Ttangent = (self.Tstress4 - self.Cstress4) / dStrain
        else:
            pass


    def _determineTrialState(self, dStrain: float):
        """Determine the trial state based on the current strain increment"""
        # Determine stage
        if -self.ugap <= self.Tstrain3 <= self.ugap:
            self.Tstage = 1  # stage-1
        else:
            self.Tstage = 2  # stage-2
        if self.ugap == 0:
            self.Tstage = 2  # If ugap is zero, always in stage-2
        if self.Tfracture:
            # SMA cable fracture
            if self.configTypeVal == 1:
                uy = self.F1 / self.k0
                if self.Tplate2 + uy <= self.Tstrain3 <= self.Tplate1 - uy:
                    self.Tstress4 = 0
                else:
                    self.Tstress4 = self._frictionModel(self.Cstress4, dStrain)
                    if dStrain < 0 and self.Tstrain3 > 0 and self.Tstress4 <= 0:
                        self.Tstress4 = 0
                    elif dStrain > 0 and self.Tstrain3 < 0 and self.Tstress4 >= 0:
                        self.Tstress4 = 0
            elif self.configTypeVal == 2:
                self.Tstress4 = self._frictionModel(self.Cstress4, dStrain)
            return
        if self.Cstage == 1 and self.Tstage == 1:
            # NOTE: stage-1 -> stage-1
            self.Tstress1 = self._frictionModel(self.Cstress, dStrain)
            self.Tstress2 = self.Tstress1
            self.Tstress = self.Tstress1
        elif self.Cstage == 1 and self.Tstage == 2:
            # NOTE: stage-1 -> stage-2
            if dStrain > 0:
                du1 = self.ugap - self.Cstrain3  # Strain increment in stage-1
                du2 = dStrain - du1  # Strain increment in stage-2
                usc0 = self.ugap - self.ua
            else:
                du1 = -self.ugap - self.Cstrain3  # Strain increment in stage-1
                du2 = dStrain - du1  # Strain increment in stage-2
                usc0 = self.ua - self.ugap
            if self.Thardening:
                self.TCDD = self.CCDD + abs(du2) / (self.uh - self.ugap)
            F1_ = self._frictionModel(self.Cstress1, du1)
            F2_ = self._SCModel(usc0, F1_, du2)
            self.Tstress1 = F2_
            # Apply degradation
            self.Tstress2 = self.Tstress1
            if self.Thardening and self.Tstrain3 > 0:
                self.Tstress2 = self.Tstress1 - (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain3) - self.ugap) / (self.uh - self.ugap))
            elif self.Thardening and self.Tstrain3 < 0:
                self.Tstress2 = self.Tstress1 + (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain3) - self.ugap) / (self.uh - self.ugap))
            # Apply modifiction
            self.Tstress = self.Tstress2
            if dStrain > 0 and self.Tstress2 < self.F1:
                self.Tstress = self.F1
            elif dStrain < 0 and self.Tstress2 > -self.F1:
                self.Tstress = -self.F1
        elif self.Cstage == 2 and self.Tstage == 2:
            # NOTE: stage-2 -> stage-2
            if self.Thardening:
                self.TCDD = self.CCDD + abs(dStrain) / (self.uh - self.ugap)
            if self.Tstrain3 >= 0:
                usc0 = self.Cstrain3 - self.ua
            else:
                usc0 = self.Cstrain3 + self.ua
            self.Tstress1 = self._SCModel(usc0, self.Cstress1, dStrain)
            # Apply degradation
            self.Tstress2 = self.Tstress1
            if self.Thardening and self.Tstrain3 > 0:
                self.Tstress2 = self.Tstress1 - (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain3) - self.ugap) / (self.uh - self.ugap))
            elif self.Thardening and self.Tstrain3 < 0:
                self.Tstress2 = self.Tstress1 + (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain3) - self.ugap) / (self.uh - self.ugap))
            # Apply modifiction
            if self.configTypeVal == 1:
                F_bound = 0  # Only half of the friction pads are sliding at stage-2
            else:
                F_bound = self.F1  # All friction pads are sliding at stage-2
            self.Tstress = self.Tstress2
            if dStrain > 0 and self.Tstrain3 > 0 and self.Tstress2 < self.F1 and self.ugap > 0 and self.Cstress == self.F1:
                self.Tstress = self.F1
            elif dStrain < 0 and self.Tstrain3 < 0 and self.Tstress2 > -self.F1 and self.ugap > 0 and self.Cstress == -self.F1:
                self.Tstress = -self.F1
            elif self.Tstrain3 > 0 and self.Tstress2 < -F_bound:
                self.Tstress = -F_bound  # Prevent positive compressive stress in SMA cables
            elif self.Tstrain3 < 0 and self.Tstress2 > F_bound:
                self.Tstress = F_bound  # Prevent negative compressive stress in SMA cables
        elif self.Cstage == 2 and self.Tstage == 1:
            # NOTE: stage-2 -> stage-1
            if dStrain < 0:
                du1 = -(self.Cstrain3 - self.ugap)
                du2 = -(self.ugap - self.Tstrain3)
                usc0 = self.Cstrain3 - self.ua
            else:
                du1 = -self.ugap - self.Cstrain3
                du2 = self.Tstrain3 + self.ugap
                usc0 = self.Cstrain3 + self.ua
            if self.Thardening:
                self.TCDD = self.CCDD + abs(du1) / (self.uh - self.ugap)
            F1_ = self._SCModel(usc0, self.Cstress1, du1)
            # Apply degradation
            F1_ideal1 = F1_
            if self.Thardening and self.Tstrain3 > 0:
                F1_ideal1 = F1_ - (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain3) - self.ugap) / (self.uh - self.ugap))
            elif self.Thardening and self.Tstrain3 < 0:
                F1_ideal1 = F1_ + (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain3) - self.ugap) / (self.uh - self.ugap))
            # Apply modifiction
            if self.configTypeVal == 1:
                F_bound = 0  # Only half of the friction pads are sliding at stage-2
            else:
                F_bound = self.F1  # All friction pads are sliding at stage-2
            F1_ = F1_ideal1
            if dStrain > 0 and self.Tstrain3 > 0 and F1_ideal1 < self.F1 and self.ugap > 0 and self.Cstress == self.F1:
                F1_ = self.F1
            elif dStrain < 0 and self.Tstrain3 < 0 and F1_ideal1 > -self.F1 and self.ugap > 0 and self.Cstress == -self.F1:
                F1_ = -self.F1
            elif self.Tstrain3 > 0 and F1_ideal1 < -F_bound:
                F1_ = -F_bound  # Prevent positive compressive stress in SMA cables
            elif self.Tstrain3 < 0 and F1_ideal1 > F_bound:
                F1_ = F_bound  # Prevent negative compressive stress in SMA cables
            F2_ = self._frictionModel(F1_, du2)
            self.Tstress1 = F2_
            self.Tstress2 = self.Tstress1
            self.Tstress = self.Tstress1
        else:
            assert False, f"Invalid state transition, Cstage={self.Cstage}, Tstage={self.Tstage}"
        F_hardening = max(abs(self.Tstrain3) - self.uh, 0) * self.k2 * self.r3  # Strength enhancement due to hardening
        self.Tstress4 = self.Tstress
        if self.Tstrain3 > 0:
            self.Tstress4 += F_hardening
        else:
            self.Tstress4 -= F_hardening
    

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
            if u < -self.F2 * (1 - self.beta) / self.k1 and F_ > self.k2 * u - self.F2 * (1 - self.beta) * (1 - self.k2 / self.k1):
                F = self.k2 * u - self.F2 * (1 - self.beta) * (1 - self.k2 / self.k1)
            elif -self.F2 * (1 - self.beta) / self.k1 <= u <= uy and F_ > self.k1 * u:
                F = self.k1 * u
            elif u > self.F2 * (1 - self.beta) / self.k1 and F_ > self.k2 * u + self.F2 - self.k2 * uy:
                F = self.k2 * u + self.F2 - self.k2 * uy
            else:
                F = F_
        else:
            if u > self.F2 * (1 - self.beta) / self.k1 and F_ < self.k2 * u + self.F2 * (1 - self.beta) * (1 - self.k2 / self.k1):
                F = self.k2 * u + self.F2 * (1 - self.beta) * (1 - self.k2 / self.k1)
            elif -uy <= u <= self.F2 * (1 - self.beta) / self.k1 and F_ < self.k1 * u:
                F = self.k1 * u
            elif u < -self.F2 * (1 - self.beta) / self.k1 and F_ < self.k2 * u - (self.F2 - self.k2 * uy):
                F = self.k2 * u - (self.F2 - self.k2 * uy)
            else:
                F = F_
        return F

    def _commitState(self):
        self.Cstrain3 = self.Tstrain3
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cstage = self.Tstage
        self.Chardening = self.Thardening
        self.Cstress1 = self.Tstress1
        self.Cstress2 = self.Tstress2
        self.Cstress4 = self.Tstress4
        self.CCDD = self.TCDD
        self.Cfracture = self.Tfracture
        self.Cplate1 = self.Tplate1
        self.Cplate2 = self.Tplate2
        self.i += 1


    def getStrain(self):
        return self.Tstrain3


    def getStress(self):
        return self.Tstress4


    def getTangent(self):
        return self.Ttangent
    

