# cython: language_level=3, boundscheck=False, wraparound=False
from libc.math cimport fabs
from libc.float cimport DBL_EPSILON


cdef class Steel01Material:
    cdef public int tag
    # 材料参数（均用 double）
    cdef double Fy, k, b
    # 计算变量
    cdef double uy
    # 历史状态变量
    cdef double Cstrain, Tstrain, Ctangent, Ttangent, Cstress, Tstress

    def __cinit__(self,
                 int tag,
                 double Fy,
                 double k,
                 double b):
        self.tag = tag
        # 初始化材料参数
        self.Fy = Fy
        self.k = k
        self.b = b
        self._check_paras()
        self._init_paras()

    cdef void _check_paras(self):
        if not self.Fy > 0:
            raise ValueError("Fy must be positive")
        if not self.k > 0:
            raise ValueError("k must be positive")

    cdef void _init_paras(self):
        self.Cstrain = 0.0
        self.Tstrain = 0.0
        self.Cstress = 0.0
        self.Tstress = 0.0
        self.Ctangent = 0.0
        self.Ttangent = 0.0
        self.uy = self.Fy / self.k

    cpdef setTrialStrain(self, double strain, double strainRate=0):
        cdef double dStrain, f
        self.Tstrain = strain
        dStrain = strain - self.Cstrain
        if fabs(dStrain) > DBL_EPSILON:
            f = self.Cstress + dStrain * self.k
            if f > self.b * self.k * (self.Tstrain - self.uy) + self.Fy:
                f = self.b * self.k * (self.Tstrain - self.uy) + self.Fy
            elif f < self.b * self.k * (self.Tstrain + self.uy) - self.Fy:
                f = self.b * self.k * (self.Tstrain + self.uy) - self.Fy
            self.Tstress = f
            self.Ttangent = (self.Tstress - self.Cstress) / dStrain
        else:
            self.Tstress = self.Cstress
            self.Ttangent = self.Ctangent

    cpdef commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent

    cpdef double getStrain(self):
        return self.Tstrain

    cpdef double getStress(self):
        return self.Tstress

    cpdef double getTangent(self):
        return self.Ttangent
