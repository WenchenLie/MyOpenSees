# cython: language_level=3, boundscheck=False, wraparound=False
from libc.math cimport fabs
from libc.float cimport DBL_EPSILON


cdef class TSSCBMaterial:
    cdef public int tag
    # 材料参数（均用 double）
    cdef double F1, k0, ugap, F2, k1, k2, beta
    cdef bint hasHardening
    cdef double uh, r1, r2, r3, uf
    cdef int configType
    cdef double up
    # 计算变量
    cdef double ua
    # 历史状态变量
    cdef int Cstage, Tstage
    cdef double Cstrain, Tstrain, Ctangent, Ttangent
    cdef bint Chardening, Thardening
    cdef double Cstress1, Tstress1, Cstress2, Tstress2, Cstress3, Tstress3, Cstress4, Tstress4
    cdef double CCDD, TCDD
    cdef bint Cfracture, Tfracture, Cfracturing, Tfracturing
    cdef double Cplate1, Tplate1, Cplate2, Tplate2, CfractureFore, TfractureFore
    cdef double Crp, Trp
    cdef int i

    def __cinit__(self,
                  int tag,
                  double F1, double k0, double ugap, double F2, double k1, double k2, double beta,
                  object hardening=None, double uh=1e16, double r1=1, double r2=1, double r3=0,
                  object minmax=None, double uf=1e16,
                  object _configType=None, int configType=1,
                  object _up=None, double up=0):
        self.tag = tag
        # 初始化材料参数
        self.F1 = F1;   self.k0 = k0;   self.ugap = ugap
        self.F2 = F2;   self.k1 = k1;   self.k2 = k2;   self.beta = beta
        if hardening is not None:
            if hardening != b'-hardening' and hardening != "-hardening":
                raise ValueError("`hardening` should be `-hardening` if given")
            self.hasHardening = True
        else:
            self.hasHardening = False
        if minmax is not None:
            if minmax != b'-minmax' and minmax != "-minmax":
                raise ValueError("`minmax` should be `-minmax` if given")
        if _configType is not None:
            if _configType != b'-configType' and _configType != "-configType":
                raise ValueError("`_configType` should be `-configType` if given")
        if _up is not None:
            if _up != b'-up' and _up != "-up":
                raise ValueError("`_up` should be `-up` if given")
        self.uh = uh; self.r1 = r1; self.r2 = r2; self.r3 = r3; self.uf = uf
        self.configType = configType
        self.up = up
        # 计算参数 ua
        self.ua = self.ugap - self.F1 / self.k1
        if self.ua < 0:
            self.ua = 0.0
        # 验证参数有效性并初始化状态
        self._check_paras()
        self._init_paras()

    cdef void _check_paras(self):
        if not (self.F1 >= 0 and self.k0 > 0 and self.ugap >= 0 and
                self.F2 > 0 and self.k1 > 0 and self.k2 > 0 and
                0 <= self.beta <= 2 and self.uh > 0 and
                self.r1 >= 0 and self.r2 >= 0 and self.r3 >= 0 and
                self.uf > 0 and self.up >= 0):
            raise ValueError("Invalid parameters for TSSCBMaterial")

    cdef void _init_paras(self):
        # 初始化历史变量和当前变量
        self.Cstage = 1; self.Tstage = 1
        self.Cstrain = 0.0; self.Tstrain = 0.0
        self.Cstress3 = 0.0; self.Tstress3 = 0.0
        if self.ugap == 0.0:
            self.Ctangent = self.k1; self.Ttangent = self.k1
            self.Cstage = 2; self.Tstage = 2
        else:
            self.Ctangent = self.k0; self.Ttangent = self.k0
        self.Chardening = False; self.Thardening = False
        self.Cstress1 = 0.0; self.Tstress1 = 0.0
        self.Cstress2 = 0.0; self.Tstress2 = 0.0
        self.Cstress4 = 0.0; self.Tstress4 = 0.0
        self.CCDD = 0.0; self.TCDD = 0.0
        self.Cfracture = False; self.Tfracture = False
        self.Cplate1 = self.ugap; self.Tplate1 = self.ugap
        self.Cplate2 = -self.ugap; self.Tplate2 = -self.ugap
        self.Cfracturing = False; self.Tfracturing = False
        self.CfractureFore = 0.0; self.TfractureFore = 0.0
        self.Crp = 0.0; self.Trp = 0.0
        self.i = 0

    cpdef setTrialStrain(self, double strain, double strainRate=0):
        cdef double dStrain
        # 将试验态变量设为上一步的收敛态
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
        # 计算应变增量
        dStrain = strain - self.Cstrain
        if fabs(dStrain) > DBL_EPSILON:
            self.Tstrain = strain
            # 判断是否开始硬化或断裂
            if fabs(self.Tstrain) > self.uh or self.Chardening:
                self.Thardening = True
            if fabs(self.Tstrain) > self.uf:
                if self.up == 0.0:
                    self.Tfracture = True
                else:
                    self.Tfracturing = True
            if self.Tfracturing:
                self.Trp = self.Crp + fabs(dStrain) / self.up
                if self.Trp >= 1.0:
                    self.Tfracture = True
                    self.Trp = 1.0
                if self.Crp == 0.0:
                    self.TfractureFore = self.Cstress4
            # 确定当前试验状态（阶段、应力等）
            self._determineTrialState(dStrain)
            # 更新端板位置
            if dStrain > 0:
                if self.Tplate1 < self.Tstrain:
                    self.Tplate1 = self.Tstrain
                if not self.Tfracture:
                    self.Tplate2 += dStrain
            else:
                if self.Tplate2 > self.Tstrain:
                    self.Tplate2 = self.Tstrain
                if not self.Tfracture:
                    self.Tplate1 += dStrain
            if self.Tplate1 < self.ugap:
                self.Tplate1 = self.ugap
            if self.Tplate2 > -self.ugap:
                self.Tplate2 = -self.ugap
            # 更新切线刚度
            self.Ttangent = (self.Tstress4 - self.Cstress4) / dStrain
        else:
            pass  # 无变化则保持状态

    cdef void _determineTrialState(self, double dStrain):
        # 根据当前试验应变确定阶段和应力
        if -self.ugap <= self.Tstrain <= self.ugap:
            self.Tstage = 1
        else:
            self.Tstage = 2
        if self.ugap == 0.0:
            self.Tstage = 2
        # 如果已经完全断裂，直接使用摩擦模型
        if self.Tfracture:
            if self.configType == 1:
                uy = self.F1 / self.k0
                if (self.Tplate2 + uy <= self.Tstrain <= self.Tplate1 - uy):
                    self.Tstress4 = 0.0
                else:
                    self.Tstress4 = self._frictionModel(self.Cstress4, dStrain, 0.5)
                    if dStrain < 0 and self.Tstrain > 0 and self.Tstress4 <= 0:
                        self.Tstress4 = 0.0
                    elif dStrain > 0 and self.Tstrain < 0 and self.Tstress4 >= 0:
                        self.Tstress4 = 0.0
            else:
                self.Tstress4 = self._frictionModel(self.Cstress4, dStrain)
            return
        # 如果正在断裂过程中
        if self.Tfracturing:
            if self.Cstress4 >= 0.0:
                self.Tstress4 = self.TfractureFore - (self.TfractureFore - self.F1) * self.Trp
            else:
                self.Tstress4 = self.TfractureFore + (-self.F1 - self.TfractureFore) * self.Trp
            return
        # 各种状态转移逻辑
        if self.Cstage == 1 and self.Tstage == 1:
            # 1->1 级：摩擦模型
            self.Tstress1 = self._frictionModel(self.Cstress3, dStrain)
            self.Tstress2 = self.Tstress1
            self.Tstress3 = self.Tstress1
        elif self.Cstage == 1 and self.Tstage == 2:
            # 1->2 级：进入第二阶段
            if dStrain > 0:
                du1 = self.ugap - self.Cstrain
                du2 = dStrain - du1
                usc0 = self.ugap - self.ua
            else:
                du1 = -self.ugap - self.Cstrain
                du2 = dStrain - du1
                usc0 = self.ua - self.ugap
            if self.Thardening:
                self.TCDD = self.CCDD + fabs(du2) / (self.uh - self.ugap)
            F1_ = self._frictionModel(self.Cstress1, du1)
            F2_ = self._SCModel(usc0, F1_, du2)
            self.Tstress1 = F2_
            self.Tstress2 = self.Tstress1
            Fd = (self.F2 - self.F1 / 2.0) * self.TCDD * (self.r1 - self.r2 * (fabs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            if self.Thardening and self.Tstrain > 0:
                self.Tstress2 = self.Tstress1 - Fd
            elif self.Thardening and self.Tstrain < 0:
                self.Tstress2 = self.Tstress1 + Fd
            self.Tstress3 = self.Tstress2
            if dStrain > 0 and self.Tstress2 < self.F1:
                self.Tstress3 = self.F1
            elif dStrain < 0 and self.Tstress2 > -self.F1:
                self.Tstress3 = -self.F1
        elif self.Cstage == 2 and self.Tstage == 2:
            # 2->2 级：双阶自复位
            if self.Thardening:
                self.TCDD = self.CCDD + fabs(dStrain) / (self.uh - self.ugap)
            if self.Tstrain >= 0:
                usc0 = self.Cstrain - self.ua
            else:
                usc0 = self.Cstrain + self.ua
            self.Tstress1 = self._SCModel(usc0, self.Cstress1, dStrain)
            self.Tstress2 = self.Tstress1
            Fd = (self.F2 - self.F1 / 2.0) * self.TCDD * (self.r1 - self.r2 * (fabs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            if self.Thardening and self.Tstrain > 0:
                self.Tstress2 = self.Tstress1 - Fd
            elif self.Thardening and self.Tstrain < 0:
                self.Tstress2 = self.Tstress1 + Fd
            if self.configType == 1:
                F_bound = 0.0
            else:
                F_bound = self.F1
            self.Tstress3 = self.Tstress2
            if dStrain > 0 and self.Tstrain > 0 and self.Tstress2 < self.F1 and self.ugap > 0 and self.Cstress3 == self.F1:
                self.Tstress3 = self.F1
            elif dStrain < 0 and self.Tstrain < 0 and self.Tstress2 > -self.F1 and self.ugap > 0 and self.Cstress3 == -self.F1:
                self.Tstress3 = -self.F1
            elif self.Tstrain > 0 and self.Tstress2 < -F_bound:
                self.Tstress3 = -F_bound
            elif self.Tstrain < 0 and self.Tstress2 > F_bound:
                self.Tstress3 = F_bound
            if dStrain > 0 and self.Tstress3 <= self.Cstress3:
                self.Tstress3 = self.Cstress4
            elif dStrain < 0 and self.Tstress3 >= self.Cstress3:
                self.Tstress3 = self.Cstress4
            if self.configType == 1 and self.Tstrain >= 0 and dStrain > 0 and self.Thardening and self.Tstress2 < self.F1:
                self.Tstress3 = self._frictionModel(self.Cstress3, dStrain)
            elif self.configType == 1 and self.Tstrain <= 0 and dStrain < 0 and self.Thardening and self.Tstress2 > self.F1:
                self.Tstress3 = self._frictionModel(self.Cstress3, dStrain)
        elif self.Cstage == 2 and self.Tstage == 1:
            # 2->1 级：从第二阶段返回第一阶段
            if dStrain < 0:
                du1 = -(self.Cstrain - self.ugap)
                du2 = -(self.ugap - self.Tstrain)
                usc0 = self.Cstrain - self.ua
            else:
                du1 = -self.ugap - self.Cstrain
                du2 = self.Tstrain + self.ugap
                usc0 = self.Cstrain + self.ua
            if self.Thardening:
                self.TCDD = self.CCDD + fabs(du1) / (self.uh - self.ugap)
            F1_ = self._SCModel(usc0, self.Cstress1, du1)
            F1_ideal1 = F1_
            if self.Thardening and self.Tstrain > 0:
                F1_ideal1 = F1_ - (self.F2 - self.F1 / 2.0) * self.TCDD * (self.r1 - self.r2 * (fabs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            elif self.Thardening and self.Tstrain < 0:
                F1_ideal1 = F1_ + (self.F2 - self.F1 / 2.0) * self.TCDD * (self.r1 - self.r2 * (fabs(self.Tstrain) - self.ugap) / (self.uh - self.ugap))
            if self.configType == 1:
                F_bound = 0.0
            else:
                F_bound = self.F1
            F1_ = F1_ideal1
            if dStrain > 0 and self.Tstrain > 0 and F1_ideal1 < self.F1 and self.ugap > 0 and self.Cstress3 == self.F1:
                F1_ = self.F1
            elif dStrain < 0 and self.Tstrain < 0 and F1_ideal1 > -self.F1 and self.ugap > 0 and self.Cstress3 == -self.F1:
                F1_ = -self.F1
            elif self.Tstrain > 0 and F1_ideal1 < -F_bound:
                F1_ = -F_bound
            elif self.Tstrain < 0 and F1_ideal1 > F_bound:
                F1_ = F_bound
            F2_ = self._frictionModel(F1_, du2)
            self.Tstress1 = F2_
            self.Tstress2 = self.Tstress1
            self.Tstress3 = self.Tstress1
        else:
            raise AssertionError(f"Invalid state transition: Cstage={self.Cstage}, Tstage={self.Tstage}")
        # 硬化部分力学修正
        cdef double F_hardening = 0.0
        if self.Tstrain > self.uh:
            F_hardening = (self.Tstrain - self.uh) * self.k2 * self.r3
        elif -self.Tstrain > self.uh:
            F_hardening = ( -self.Tstrain - self.uh) * self.k2 * self.r3
        if self.Tstrain >= 0:
            self.Tstress4 = self.Tstress3 + F_hardening
        else:
            self.Tstress4 = self.Tstress3 - F_hardening

    cdef double _frictionModel(self, double F0, double du, double half=1.0):
        # 滑动摩擦力模型
        cdef double F_, F
        if du == 0.0:
            return F0
        F_ = F0 + du * self.k0
        if F_ > self.F1 * half:
            F = self.F1 * half
        elif F_ < -self.F1 * half:
            F = -self.F1 * half
        else:
            F = F_
        return F

    cdef double _SCModel(self, double u0, double F0, double du):
        # 双阶自复位模型
        cdef double u, uy, F_, F
        if du == 0.0:
            return F0
        if self.Tfracture:
            return 0.0
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

    cpdef commitState(self):
        # 提交当前状态为收敛状态
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

    cpdef double getStress(self):
        return self.Tstress4

    cpdef double getTangent(self):
        return self.Ttangent

    cpdef setStrain(self, double strain):
        self.setTrialStrain(strain)
        self.commitState()