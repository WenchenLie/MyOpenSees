# cython: language_level=3
# cython: boundscheck=False, wraparound=False, cdivision=True, nonecheck=False

from libc.math cimport fabs
from libc.float cimport DBL_EPSILON

cdef inline double dmax(double a, double b) noexcept:
    return a if a > b else b

cdef inline double dmin(double a, double b) noexcept:
    return a if a < b else b

cdef class TwoStage:
    cdef public int tag

    # Material parameters
    cdef double F1, k1, kp1
    cdef double F2, k2, kp2
    cdef double ua

    # Response history
    cdef double Cstrain, Tstrain
    cdef double Cstrain2, Tstrain2
    cdef double Cstress, Tstress
    cdef double Cstress1, Tstress1
    cdef double Cstress2, Tstress2
    cdef double Ctangent, Ttangent
    cdef double Chookgap, Thookgap

    def __init__(self,
                 int tag,
                 double F1, double k1, double kp1,
                 double F2, double k2, double kp2,
                 double ua):
        self.tag = tag
        self.F1 = F1
        self.k1 = k1
        self.kp1 = kp1
        self.F2 = F2
        self.k2 = k2
        self.kp2 = kp2
        self.ua = ua
        self._check_paras()
        self._init_paras()

    cdef void _check_paras(self) except *:
        if self.ua < 0.0: raise ValueError("ua must be non-negative")

    cdef void _init_paras(self):
        self.Cstrain = 0.0
        self.Tstrain = 0.0
        self.Cstrain2 = 0.0
        self.Tstrain2 = 0.0
        self.Cstress = 0.0
        self.Tstress = 0.0
        self.Cstress1 = 0.0
        self.Tstress1 = 0.0
        self.Cstress2 = 0.0
        self.Tstress2 = 0.0
        self.Ctangent = self.k1
        self.Ttangent = self.k1
        self.Chookgap = 0.0
        self.Thookgap = 0.0

    @staticmethod
    cdef inline double bilinear(double F_prev, double u_prev, double du,
                                double Fy, double k, double kp) noexcept:
        cdef double F_next
        if Fy == 0.0:
            return 0.0
        F_next = F_prev + du * k
        if du > 0.0:
            if F_next > kp * (u_prev + du) + (1.0 - kp / k) * Fy:
                F_next = kp * (u_prev + du) + (1.0 - kp / k) * Fy
        elif du < 0.0:
            if F_next < kp * (u_prev + du) - (1.0 - kp / k) * Fy:
                F_next = kp * (u_prev + du) - (1.0 - kp / k) * Fy
        return F_next

    cpdef void setTrialStrain(self, double strain, double strainRate=0.0):
        """传入当前步的应变值 strain"""
        cdef double dStrain, dstrain2
        cdef double tmp, denom

        dStrain = strain - self.Cstrain
        self.Tstrain = strain
        if fabs(dStrain) <= DBL_EPSILON:
            return

        # 一阶单元（总是变形）
        self.Tstress1 = TwoStage.bilinear(self.Cstress1, self.Cstrain, dStrain,
                                          self.F1, self.k1, self.kp1)
        # 二阶单元（钩距触发）
        if (-self.ua < self.Thookgap) and (self.Thookgap < self.ua):
            if dStrain > 0.0:
                tmp = self.Thookgap + dStrain - self.ua
                dstrain2 = tmp if tmp > 0.0 else 0.0
                self.Thookgap = dmin(self.Thookgap + dStrain, self.ua)
            else:
                tmp = self.Thookgap + dStrain + self.ua
                dstrain2 = tmp if tmp < 0.0 else 0.0
                self.Thookgap = dmax(self.Thookgap + dStrain, -self.ua)
            self.Tstress2 = TwoStage.bilinear(self.Cstress2, self.Cstrain2, dstrain2,
                                              self.F2, self.k2, self.kp2)
            self.Tstrain2 = self.Cstrain2 + dstrain2
        elif self.Thookgap == self.ua:
            self.Tstress2 = TwoStage.bilinear(self.Cstress2, self.Cstrain2, dStrain,
                                              self.F2, self.k2, self.kp2)
            self.Tstrain2 = self.Cstrain2 + dStrain
            if dStrain <= 0.0 and self.Tstress2 < 0.0:
                denom = fabs(self.Tstress2) + fabs(self.Cstress2)
                if denom > DBL_EPSILON:
                    dstrain2 = dStrain * fabs(self.Cstress2) / denom
                    self.Thookgap = self.ua + (dStrain - dstrain2)
                    if self.Thookgap < -self.ua:
                        self.Thookgap = -self.ua
                self.Tstress2 = 0.0
        elif self.Thookgap == -self.ua:
            self.Tstress2 = TwoStage.bilinear(self.Cstress2, self.Cstrain2, dStrain,
                                              self.F2, self.k2, self.kp2)
            self.Tstrain2 = self.Cstrain2 + dStrain
            if dStrain >= 0.0 and self.Tstress2 > 0.0:
                denom = fabs(self.Tstress2) + fabs(self.Cstress2)
                if denom > DBL_EPSILON:
                    dstrain2 = dStrain * fabs(self.Cstress2) / denom
                    self.Thookgap = -self.ua + (dStrain - dstrain2)
                    if self.Thookgap > self.ua:
                        self.Thookgap = self.ua
                self.Tstress2 = 0.0
        else:
            raise AssertionError(f"Should not reach here (Thookgap = {self.Thookgap}, ua = {self.ua})")
        # 总应力
        self.Tstress = self.Tstress1 + self.Tstress2
        self.Ttangent = (self.Tstress - self.Cstress) / dStrain

    cpdef void commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstrain2 = self.Tstrain2
        self.Cstress = self.Tstress
        self.Cstress1 = self.Tstress1
        self.Cstress2 = self.Tstress2
        self.Ctangent = self.Ttangent
        self.Chookgap = self.Thookgap

    cpdef double getStrain(self):
        return self.Tstrain

    cpdef double getStress(self):
        return self.Tstress

    cpdef double getTangent(self):
        return self.Ttangent

    cpdef double setStrain(self, double strain, double strainRate=0.0):
        self.setTrialStrain(strain, strainRate)
        self.commitState()
