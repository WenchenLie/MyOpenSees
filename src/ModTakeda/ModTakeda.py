import sys
from ..UniaxialMaterial import UniaxialMaterial


class ModTakeda(UniaxialMaterial):
    """Modified Tekada model"""
    
    def __init__(self, tag: int,
                 Fy: float,
                 k0: float,
                 r: float,
                 alpha: float,
                 beta: float
    ):
        self.tag = tag
        # Materail parameters
        self.Fy = Fy
        self.k0 = k0
        self.r = r
        self.alpha = alpha
        self.beta = beta
        # Response history
        self.Cstrain: float  # 上一步的应变
        self.Tstrain: float  # 当前步的应变
        self.Cstress: float  # 上一步的应力
        self.Tstress: float # 当前步的应力
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self.Cdm_pos: float  # 历史最大正应变
        self.Tdm_pos: float  # 当前最大正应变
        self.Cdm_neg: float  # 历史最大负应变
        self.Tdm_neg: float  # 当前最大负应变
        self.CFm_pos: float  # 历史最大正应力
        self.TFm_pos: float  # 当前最大正应力
        self.CFm_neg: float  # 历史最大负应力
        self.TFm_neg: float  # 当前最大负应力
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        assert self.Fy >= 0
        assert self.k0 >= 0
        assert self.r >= 0
        assert self.alpha >= 0
        assert self.beta >= 0

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
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

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain"""
        # Reset history variables to last converged state
        self.Tstrain = strain
        dStrain = strain - self.Cstrain
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
                        kr = F_flag / (u_flag - u0)  # 再加载刚度
                        self.Tstress = kr * dStrain2
                    else:
                        self.Tstress = self.Cstress + ku * dStrain
                else:
                    if u_flag > strain:
                        kr = (F_flag - self.Cstress) / (u_flag - self.Cstrain)  # 再加载刚度
                        self.Tstress = self.Cstress + kr * dStrain
                    else:
                        self.Tstress = self.Fy + (strain - self.uy) * self.r * self.k0
            else:
                u_flag = min(-self.uy, self.Cdm_neg - self.beta * (self.Cdm_neg + self.uy))
                F_flag = min(-self.Fy, self.CFm_neg - self.beta * (self.Cdm_neg + self.uy) * self.r * self.k0)
                if self.Cstress > 0:
                    ku = self.k0 * abs(self.uy / self.Cdm_neg) ** self.alpha  # 卸载刚度
                    if self.Cstress + ku * dStrain < 0:
                        dStrain1 = -self.Cstress / ku
                        dStrain2 = dStrain - dStrain1
                        u0 = self.Cstrain + dStrain1  # 滞回曲线与x轴交点横坐标
                        kr = F_flag / (u_flag - u0)  # 再加载刚度
                        self.Tstress = kr * dStrain2
                    else:
                        self.Tstress = self.Cstress + ku * dStrain
                else:
                    if u_flag < strain:
                        kr = (F_flag - self.Cstress) / (u_flag - self.Cstrain)  # 再加载刚度
                        self.Tstress = self.Cstress + kr * dStrain
                    else:
                        self.Tstress = -self.Fy + (strain + self.uy) * self.r * self.k0
            # 更新Flag点
            if dStrain > 0:
                self.Tdm_pos = max(self.Tdm_pos, self.Tstrain)
                self.TFm_pos = max(self.TFm_pos, self.Tstress)
            else:
                self.Tdm_neg = min(self.Tdm_neg, self.Tstrain)
                self.TFm_neg = min(self.TFm_neg, self.Tstress)
            # 更新切线刚度
            self.Ttangent = (self.Tstress - self.Cstress) / dStrain
        else:
            self.Tstress = self.Cstress
            self.Ttangent = self.Ctangent

    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cdm_pos = self.Tdm_pos
        self.Cdm_neg = self.Tdm_neg
        self.CFm_pos = self.TFm_pos
        self.CFm_neg = self.TFm_neg

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent
    

