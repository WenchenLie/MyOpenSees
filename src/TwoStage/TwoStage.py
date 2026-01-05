import sys
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
        self.Cstrain2: float  # 二阶单元应变
        self.Tstrain2: float
        self.Cstress: float  # 总应力
        self.Tstress: float
        self.Cstress1: float  # 一阶单元应力(总是变形)
        self.Tstress1: float
        self.Cstress2: float  # 二阶单元应力(触发变形)
        self.Tstress2: float
        self.Ctangent: float  # 切线刚度
        self.Ttangent: float
        self.Chookgap: float  # 钩距(始终处于[-ua, ua]范围内)
        self.Thookgap: float
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        assert self.ua >= 0, "ua must be non-negative"

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstrain2 = 0
        self.Tstrain2 = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Cstress1 = 0
        self.Tstress1 = 0
        self.Cstress2 = 0
        self.Tstress2 = 0
        self.Ctangent = self.k1
        self.Ttangent = self.k1
        self.Chookgap = 0
        self.Thookgap = 0

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        dstrain2: float
        dStrain = strain - self.Cstrain
        self.Tstrain = strain
        if abs(dStrain) <= sys.float_info.epsilon:
            return
        # 计算一阶单元应力
        self.Tstress1 = self.bilinear(self.Cstress1, self.Cstrain, dStrain, self.F1, self.k1, self.kp1)
        # 计算二阶单元应力
        if -self.ua < self.Thookgap < self.ua:
            # 原本在钩距内
            if dStrain > 0:
                dstrain2 = max(self.Thookgap + dStrain - self.ua, 0)
                self.Thookgap = min(self.Thookgap + dStrain, self.ua)
            else:
                dstrain2 = min(self.Thookgap + dStrain + self.ua, 0)
                self.Thookgap = max(self.Thookgap + dStrain, -self.ua)
            self.Tstress2 = self.bilinear(self.Cstress2, self.Cstrain2, dstrain2, self.F2, self.k2, self.kp2)
            self.Tstrain2 = self.Cstrain2 + dstrain2
        elif self.Thookgap == self.ua:
            # 钩距已经到达最大值
            if dStrain > 0:
                self.Tstress2 = self.bilinear(self.Cstress2, self.Cstrain2, dStrain, self.F2, self.k2, self.kp2)
                self.Tstrain2 = self.Cstrain2 + dStrain
            else:
                self.Tstress2 = self.bilinear(self.Cstress2, self.Cstrain2, dStrain, self.F2, self.k2, self.kp2)
                self.Tstrain2 = self.Cstrain2 + dStrain
                if self.Tstress2 < 0:
                    dstrain2 = dStrain * abs(self.Cstress2) / (abs(self.Tstress2) + abs(self.Cstress2))
                    self.Thookgap = self.ua + (dStrain - dstrain2)
                    if self.Thookgap < -self.ua:
                        self.Thookgap = -self.ua
                    self.Tstress2 = 0
        elif self.Thookgap == -self.ua:
            # 钩距已经到达最小值
            if dStrain < 0:
                self.Tstress2 = self.bilinear(self.Cstress2, self.Cstrain2, dStrain, self.F2, self.k2, self.kp2)
                self.Tstrain2 = self.Cstrain2 + dStrain
            else:
                self.Tstress2 = self.bilinear(self.Cstress2, self.Cstrain2, dStrain, self.F2, self.k2, self.kp2)
                self.Tstrain2 = self.Cstrain2 + dStrain
                if self.Tstress2 > 0:
                    dstrain2 = dStrain * abs(self.Cstress2) / (abs(self.Tstress2) + abs(self.Cstress2))
                    self.Thookgap = -self.ua + (dStrain - dstrain2)
                    if self.Thookgap > self.ua:
                        self.Thookgap = self.ua
                    self.Tstress2 = 0
        else:
            assert False, f"Should not reach here (Thookgap = {self.Thookgap}, ua = {self.ua})"
        # 总应力
        self.Tstress = self.Tstress1 + self.Tstress2
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
            u (float): 下一步位移位移
            du (float): 位移增量
            Fy (float): 屈服力
            k (float): 弹性刚度
            kp (float): 塑性模量

        Returns:
            float: 当前步应力
        """
        if Fy == 0:
            return 0
        F_next = F_prev + du * k
        if F_next > kp * (u_prev + du) + (1 - kp / k) * Fy and du > 0:
            F_next = kp * (u_prev + du) + (1 - kp / k) * Fy
        elif F_next < kp * (u_prev + du) - (1 - kp / k) * Fy and du < 0:
            F_next = kp * (u_prev + du) - (1 - kp / k) * Fy
        return F_next
    
    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstrain2 = self.Tstrain2
        self.Cstress = self.Tstress
        self.Cstress1 = self.Tstress1
        self.Cstress2 = self.Tstress2
        self.Ctangent = self.Ttangent
        self.Chookgap = self.Thookgap

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent
    

