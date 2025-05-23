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
        _configType: Literal['-configType', None]=None,
        configType: Literal[1, 2]=1,
        _up: Literal['-up', None]=None,
        up: float=0
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
            _configType (Literal['-configType', None]): 退化类型，默认1
            configType (Literal[1, 2]): 退化类型, 1: 第二阶段只有一半摩擦片滑动, 
                2: 第二阶段所有摩擦片均滑动
            _up (Literal['-up', None]): 是否允许断裂起始到完全断裂间平滑过渡(可改善收敛性)
            up (float): 断裂起始到完全断裂间允许的位移量，默认0，即断裂瞬间承载力会突降，
                当`up`不为0时，断裂开始后，在累积变形达到`up`时承载力才会完全下降
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
        if _configType is not None:
            if _configType != '-configType':
                raise ValueError('`_configType` should be `-configType` if it is given')
        if _up is not None:
            if _up != '-up':
                raise ValueError('`_up` should be `-up` if it is given')
        self.uh = uh
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.uf = uf
        self.configType = configType
        self.up = up
        # Calculated parameters
        self.ua = self.ugap - self.F1 / self.k1
        self.ua = max(0, self.ua)  # 第一阶段进入第二阶段时自复位分量的初始应变
        # Response history
        self.Cstage: int  # 上一步所处阶段
        self.Tstage: int  # 当前步所处阶段
        self.Cstrain: float  # 应变
        self.Tstrain: float
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self.Chardening: bool  # 是否开始硬化
        self.Thardening: bool
        self.Cstress1: float  # 无退化，无修正应力
        self.Tstress1: float
        self.Cstress2: float  # 有退化，无修正应力
        self.Tstress2: float
        self.Cstress3: float  # 有退化，有修正后应力
        self.Tstress3: float
        self.Tstress4: float  # 考虑硬化提升的强度，最终返回的应力
        self.Cstress4: float
        self.CCDD: float
        self.TCDD: float
        self.Cfracture: bool  # 是否已经断裂
        self.Tfracture: bool
        self.Cplate1: float  # 左端板位置
        self.Tplate1: float
        self.Cplate2: float  # 右端板位置
        self.Tplate2: float
        self.Cfracturing: bool  # 是否正处于断裂过程中
        self.Tfracturing: bool
        self.CfractureFore: float  # 断裂起始时的强度
        self.TfractureFore: float
        self.Crp: float  # 断裂进度(0-1，0代表断裂起始，1代表完全断裂)
        self.Trp: float
        self._check_paras()
        self._init_paras()
    
    def _init_paras(self):
        self.Cstage = 1
        self.Tstage = 1
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress3 = 0
        self.Tstress3 = 0
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
        self.Cfracturing = False
        self.Tfracturing = False
        self.CfractureFore = 0
        self.TfractureFore = 0
        self.Crp = 0
        self.Trp = 0
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
        assert self.up >= 0

    def _setTrainStrain(self, strain, strainRate: float=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = self.Cstrain
        self.Tstress3 = self.Cstress3
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
        dStrain = strain - self.Cstrain
        if abs(dStrain) > sys.float_info.epsilon:
            self.Tstrain = strain
            # Determine whether to start hardening
            if abs(self.Tstrain) > self.uh or self.Chardening:
                self.Thardening = True
            if abs(self.Tstrain) > self.uf:
                if self.up == 0:
                    self.Tfracture = True
                else:
                    self.Tfracturing = True
            if self.Tfracturing:
                self.Trp = self.Crp + abs(dStrain) / self.up
                if self.Trp >= 1:
                    self.Tfracture = True
                    self.Trp = 1
                if self.Crp == 0:
                    self.TfractureFore = self.Cstress4  # 记录首次进入断裂状态时的力
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
            self.Ttangent = (self.Tstress4 - self.Cstress4) / dStrain
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
            # SMA cable fracture completed
            if self.configType == 1:
                uy = self.F1 / self.k0
                if self.Tplate2 + uy <= self.Tstrain <= self.Tplate1 - uy:
                    self.Tstress4 = 0
                else:
                    self.Tstress4 = self._frictionModel(self.Cstress4, dStrain)
                    if dStrain < 0 and self.Tstrain > 0 and self.Tstress4 <= 0:
                        self.Tstress4 = 0
                    elif dStrain > 0 and self.Tstrain < 0 and self.Tstress4 >= 0:
                        self.Tstress4 = 0
            elif self.configType == 2:
                self.Tstress4 = self._frictionModel(self.Cstress4, dStrain)
            return
        if self.Tfracturing:
            # SMA cable fracture starts
            if self.Cstress4 >= 0:
                self.Tstress4 = self.TfractureFore - (self.TfractureFore - self.F1) * self.Trp
            else:
                self.Tstress4 = self.TfractureFore + (-self.F1 - self.TfractureFore) * self.Trp
            return
        if self.Cstage == 1 and self.Tstage == 1:
            # NOTE: stage-1 -> stage-1
            self.Tstress1 = self._frictionModel(self.Cstress3, dStrain)
            self.Tstress2 = self.Tstress1
            self.Tstress3 = self.Tstress1
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
            if self.Thardening:
                self.TCDD = self.CCDD + abs(du2) / (self.uh - self.ugap)
            F1_ = self._frictionModel(self.Cstress1, du1)
            F2_ = self._SCModel(usc0, F1_, du2)
            self.Tstress1 = F2_
            # Apply degradation
            self.Tstress2 = self.Tstress1
            Fd = (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            if self.Thardening and self.Tstrain > 0:
                self.Tstress2 = self.Tstress1 - Fd
            elif self.Thardening and self.Tstrain < 0:
                self.Tstress2 = self.Tstress1 + Fd
            # Apply modifiction
            self.Tstress3 = self.Tstress2
            if dStrain > 0 and self.Tstress2 < self.F1:
                self.Tstress3 = self.F1
            elif dStrain < 0 and self.Tstress2 > -self.F1:
                self.Tstress3 = -self.F1
        elif self.Cstage == 2 and self.Tstage == 2:
            # NOTE: stage-2 -> stage-2
            if self.Thardening:
                self.TCDD = self.CCDD + abs(dStrain) / (self.uh - self.ugap)
            if self.Tstrain >= 0:
                usc0 = self.Cstrain - self.ua
            else:
                usc0 = self.Cstrain + self.ua
            self.Tstress1 = self._SCModel(usc0, self.Cstress1, dStrain)
            # Apply degradation
            self.Tstress2 = self.Tstress1
            Fd = (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            if self.Thardening and self.Tstrain > 0:
                self.Tstress2 = self.Tstress1 - Fd
            elif self.Thardening and self.Tstrain < 0:
                self.Tstress2 = self.Tstress1 + Fd
            # Apply modifiction
            if self.configType == 1:
                F_bound = 0  # Only half of the friction pads are sliding at stage-2
            else:
                F_bound = self.F1  # All friction pads are sliding at stage-2
            self.Tstress3 = self.Tstress2
            if dStrain > 0 and self.Tstrain > 0 and self.Tstress2 < self.F1 and self.ugap > 0 and self.Cstress3 == self.F1:
                self.Tstress3 = self.F1
            elif dStrain < 0 and self.Tstrain < 0 and self.Tstress2 > -self.F1 and self.ugap > 0 and self.Cstress3 == -self.F1:
                self.Tstress3 = -self.F1
            elif self.Tstrain > 0 and self.Tstress2 < -F_bound:
                self.Tstress3 = -F_bound  # Prevent positive compressive stress in SMA cables
            elif self.Tstrain < 0 and self.Tstress2 > F_bound:
                self.Tstress3 = F_bound  # Prevent negative compressive stress in SMA cables
            if dStrain > 0 and self.Tstress3 <= self.Cstress3:
                self.Tstress3 = self.Cstress4
            elif dStrain < 0 and self.Tstress3 >= self.Cstress3:
                self.Tstress3 = self.Cstress4
            if self.configType == 1 and self.Tstrain >= 0 and dStrain > 0 and self.Cstress3 == 0 and self.Thardening and self.Tstress2 < self.F1:
                self.Tstress3 = self._frictionModel(self.Cstress3, dStrain)
            elif self.configType == 1 and self.Tstrain <= 0 and dStrain < 0 and self.Cstress3 == 0 and self.Thardening and self.Tstress2 > self.F1:
                self.Tstress3 = self._frictionModel(self.Cstress3, dStrain)
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
            if self.Thardening:
                self.TCDD = self.CCDD + abs(du1) / (self.uh - self.ugap)
            F1_ = self._SCModel(usc0, self.Cstress1, du1)
            # Apply degradation
            F1_ideal1 = F1_
            if self.Thardening and self.Tstrain > 0:
                F1_ideal1 = F1_ - (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            elif self.Thardening and self.Tstrain < 0:
                F1_ideal1 = F1_ + (self.F2 - self.F1 / 2) * self.TCDD * (self.r1 - self.r2 * (abs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            # Apply modifiction
            if self.configType == 1:
                F_bound = 0  # Only half of the friction pads are sliding at stage-2
            else:
                F_bound = self.F1  # All friction pads are sliding at stage-2
            F1_ = F1_ideal1
            if dStrain > 0 and self.Tstrain > 0 and F1_ideal1 < self.F1 and self.ugap > 0 and self.Cstress3 == self.F1:
                F1_ = self.F1
            elif dStrain < 0 and self.Tstrain < 0 and F1_ideal1 > -self.F1 and self.ugap > 0 and self.Cstress3 == -self.F1:
                F1_ = -self.F1
            elif self.Tstrain > 0 and F1_ideal1 < -F_bound:
                F1_ = -F_bound  # Prevent positive compressive stress in SMA cables
            elif self.Tstrain < 0 and F1_ideal1 > F_bound:
                F1_ = F_bound  # Prevent negative compressive stress in SMA cables
            F2_ = self._frictionModel(F1_, du2)
            self.Tstress1 = F2_
            self.Tstress2 = self.Tstress1
            self.Tstress3 = self.Tstress1
        else:
            assert False, f"Invalid state transition, Cstage={self.Cstage}, Tstage={self.Tstage}"
        F_hardening = max(abs(self.Tstrain) - self.uh, 0) * self.k2 * self.r3  # Strength enhancement due to hardening
        self.Tstress4 = self.Tstress3
        if self.Tstrain > 0:
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
        if self.Tfracturing:
            return F * self.Trp
        else:
            return F

    def _commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress3 = self.Tstress3
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
        self.Cfracturing = self.Tfracturing
        self.CfractureFore = self.TfractureFore
        self.Crp = self.Trp
        self.i += 1

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress4

    def getTangent(self):
        return self.Ttangent
    

