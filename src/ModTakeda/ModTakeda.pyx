# cython: language_level=3
import sys

cdef class ModTakeda:
    """
    Modified Takeda model (Cython version)
    """
    cdef public int tag
    cdef public double Fy, k0, r, alpha, beta
    cdef public double uy
    cdef public double Cstrain, Tstrain, Cstress, Tstress, Ctangent, Ttangent
    cdef public double Cdm_pos, Tdm_pos, Cdm_neg, Tdm_neg
    cdef public double CFm_pos, TFm_pos, CFm_neg, TFm_neg

    def __init__(self, int tag,
                 double Fy,
                 double k0,
                 double r,
                 double alpha,
                 double beta):
        self.tag = tag
        self.Fy = Fy
        self.k0 = k0
        self.r = r
        self.alpha = alpha
        self.beta = beta
        self._check_paras()
        self._init_paras()

    cdef void _check_paras(self):
        assert self.Fy >= 0, "Fy must be non-negative"
        assert self.k0 >= 0, "k0 must be non-negative"
        assert self.r >= 0, "r must be non-negative"
        assert self.alpha >= 0, "alpha must be non-negative"
        assert self.beta >= 0, "beta must be non-negative"

    cdef void _init_paras(self):
        self.Cstrain = 0.0
        self.Tstrain = 0.0
        self.Cstress = 0.0
        self.Tstress = 0.0
        self.Ctangent = self.k0
        self.Ttangent = self.k0
        self.uy = self.Fy / self.k0
        self.Cdm_pos = self.Fy / self.k0
        self.Tdm_pos = self.Fy / self.k0
        self.Cdm_neg = -self.Fy / self.k0
        self.Tdm_neg = -self.Fy / self.k0
        self.CFm_pos = self.Fy
        self.TFm_pos = self.Fy
        self.CFm_neg = -self.Fy
        self.TFm_neg = -self.Fy

    cpdef void setTrialStrain(self, double strain, double strainRate=0):
        """
        Update the state with the given trial strain.
        """
        # Reset history variables to last converged state
        self.Tstrain = strain
        cdef double dStrain = strain - self.Cstrain
        cdef double dStrain1, dStrain2, u0, ku, kr, u_flag, F_flag

        if abs(dStrain) > sys.float_info.epsilon:
            if dStrain > 0:
                u_flag = max(self.uy, self.Cdm_pos - self.beta * (self.Cdm_pos - self.uy))
                F_flag = max(self.Fy, self.CFm_pos - self.beta * (self.Cdm_pos - self.uy) * self.r * self.k0)
                if self.Cstress < 0:
                    ku = self.k0 * abs(self.uy / self.Cdm_pos) ** self.alpha  # 卸载刚度
                    if self.Cstress + ku * dStrain > 0:
                        dStrain1 = -self.Cstress / ku
                        dStrain2 = dStrain - dStrain1
                        u0 = self.Cstrain + dStrain1  # 滞回曲线与x轴交点横坐标
                        if strain < u_flag:
                            kr = F_flag / (u_flag - u0)  # 再加载刚度
                        else:
                            kr = self.k0
                        self.Tstress = kr * dStrain2
                    else:
                        self.Tstress = self.Cstress + ku * dStrain
                else:
                    if u_flag > strain:
                        kr = (F_flag - self.Cstress) / (u_flag - self.Cstrain)  # 再加载刚度
                        self.Tstress = self.Cstress + kr * dStrain
                    else:
                        self.Tstress = self.Cstress + dStrain * self.k0
                if self.Tstress > self.r * self.k0 * (strain - self.uy) + self.Fy:
                    self.Tstress = self.r * self.k0 * (strain - self.uy) + self.Fy
            else:  # dStrain < 0
                u_flag = min(-self.uy, self.Cdm_neg - self.beta * (self.Cdm_neg + self.uy))
                F_flag = min(-self.Fy, self.CFm_neg - self.beta * (self.Cdm_neg + self.uy) * self.r * self.k0)
                if self.Cstress > 0:
                    ku = self.k0 * abs(self.uy / self.Cdm_neg) ** self.alpha  # 卸载刚度
                    if self.Cstress + ku * dStrain < 0:
                        dStrain1 = -self.Cstress / ku
                        dStrain2 = dStrain - dStrain1
                        u0 = self.Cstrain + dStrain1  # 滞回曲线与x轴交点横坐标
                        if strain > u_flag:
                            kr = F_flag / (u_flag - u0)  # 再加载刚度
                        else:
                            kr = self.k0
                        self.Tstress = kr * dStrain2
                    else:
                        self.Tstress = self.Cstress + ku * dStrain
                else:
                    if u_flag < strain:
                        kr = (F_flag - self.Cstress) / (u_flag - self.Cstrain)  # 再加载刚度
                        self.Tstress = self.Cstress + kr * dStrain
                    else:
                        self.Tstress = self.Cstress + dStrain * self.k0
                if self.Tstress < self.r * self.k0 * (strain + self.uy) - self.Fy:
                    self.Tstress = self.r * self.k0 * (strain + self.uy) - self.Fy

            # 更新Flag点
            if dStrain > 0:
                if self.Tstrain > self.Tdm_pos:
                    self.Tdm_pos = self.Tstrain
                if self.Tstress > self.TFm_pos:
                    self.TFm_pos = self.Tstress
            else:
                if self.Tstrain < self.Tdm_neg:
                    self.Tdm_neg = self.Tstrain
                if self.Tstress < self.TFm_neg:
                    self.TFm_neg = self.Tstress

            # 更新切线刚度
            self.Ttangent = (self.Tstress - self.Cstress) / dStrain
        else:
            self.Tstress = self.Cstress
            self.Ttangent = self.Ctangent

    cpdef void commitState(self):
        """Commit the current state as converged."""
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cdm_pos = self.Tdm_pos
        self.Cdm_neg = self.Tdm_neg
        self.CFm_pos = self.TFm_pos
        self.CFm_neg = self.TFm_neg
    
    cpdef void setStrain(self, double strain, double strainRate=0):
        self.setTrialStrain(strain, strainRate)
        self.commitState()

    cpdef double getStrain(self):
        return self.Tstrain

    cpdef double getStress(self):
        return self.Tstress

    cpdef double getTangent(self):
        return self.Ttangent