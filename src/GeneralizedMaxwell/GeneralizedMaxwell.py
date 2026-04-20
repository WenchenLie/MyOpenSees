import math
import sys
import numpy as np
import matplotlib.pyplot as plt
from ..UniaxialMaterial import UniaxialMaterial


class GeneralizedMaxwell(UniaxialMaterial):
    
    def __init__(self,
        tag: int,
        k0: float,
        k1: float,
        c1: float,
        alpha1: float,
        *args,
        _iter: str | None = '-iter',
        n_iter: int = 10
    ):
        self.tag = tag
        # Materail parameters
        self.k0 = k0
        if len(args) % 3 != 0:
            raise ValueError(f'The number of arguments should be a multiple of 3 if given')
        self.n_layer = len(args) // 3 + 1
        self.k_ls = [k1] + [args[i] for i in range(0, len(args), 3)]
        self.c_ls = [c1] + [args[i+1] for i in range(0, len(args), 3)]
        self.alpha_ls = [alpha1] + [args[i+2] for i in range(0, len(args), 3)]
        self._iter = _iter
        self.n_iter = n_iter
        # Response history
        self.Cstrain: float  # 上一步的应变
        self.Tstrain: float  # 当前步的应变
        self.Cstress: float  # 上一步的应力
        self.Tstress: float # 当前步的应力
        self.Ctangent: float  # 上一步的切线刚度
        self.Ttangent: float  # 当前步的切线刚度
        self._check_paras()
        self._init_paras()

    def _check_paras(self):
        for k in self.k_ls:
            if k <= 0:
                raise ValueError(f'ki should be positive, but got {k}')
        for c in self.c_ls:
            if c <= 0:
                raise ValueError(f'ci should be positive, but got {c}')
        for alpha in self.alpha_ls:
            if alpha < 0 or alpha > 1:
                raise ValueError(f'alpha should be between 0 and 1, but got {alpha}')
        if self._iter not in ['-iter', None]:
            raise ValueError(f'_iter should be "-iter" if given, but got {self._iter}')
        if self.n_iter <= 0:
            raise ValueError(f'n_iter should be positive, but got {self.n_iter}')
        if not isinstance(self.n_iter, int):
            raise ValueError(f'n_iter should be an integer, but got {self.n_iter}')

    def _init_paras(self):
        self.Cstrain = 0
        self.Tstrain = 0
        self.Cstress = 0
        self.Tstress = 0
        self.Ctangent = 0
        self.Ttangent = 0
        
        self.Ctangent = self.k0 + sum(self.k_ls)
        self.Ttangent = self.Ctangent
        
        # 记录每个Maxwell分支 (弹簧-阻尼器串联) 的内部应力
        self.Cstress_i = [0.0] * self.n_layer
        self.Tstress_i = [0.0] * self.n_layer

    def _compute_dS(self, S_list, d_eps, dt):
        """计算内部各分支应力关于无量纲积分变量的导数"""
        dS_list = [0.0] * self.n_layer
        for i in range(self.n_layer):
            S_val = S_list[i]
            k_val = self.k_ls[i]
            c_val = self.c_ls[i]
            alpha_val = self.alpha_ls[i]
            
            if S_val > 0:
                sign_S = 1.0
            elif S_val < 0:
                sign_S = -1.0
            else:
                sign_S = 0.0
                
            abs_S = abs(S_val)
            # 计算粘性应变率项
            dashpot_vel = sign_S * math.pow(abs_S / c_val, 1.0 / alpha_val)
            dS_list[i] = k_val * d_eps - dt * k_val * dashpot_vel
        return dS_list

    def setTrialStrain(self, strain, strainRate=0):
        """传入当前步的应变值strain，使用RK4法更新内部状态"""
        # 更新总应变增量
        self.Tstrain = strain
        d_eps = self.Tstrain - self.Cstrain

        # 根据应变和应变率推导当前步的时间增量 dt
        # 如果没有提供应变率或应变率为0，则视为瞬态加载(纯弹性行为)，dt=0
        if strainRate != 0:
            dt = abs(d_eps / strainRate)
        else:
            dt = 0.0

        # 获取上一步各分支的收敛应力
        S = list(self.Cstress_i)
        
        # RK4积分的步长 (积分域为无量纲的x，从0到1)
        h = 1.0 / self.n_iter

        # 四阶龙格库塔方法 (RK4) 进行子步积分
        for _ in range(self.n_iter):
            k1 = self._compute_dS(S, d_eps, dt)
            
            S_k2 = [S[i] + 0.5 * h * k1[i] for i in range(self.n_layer)]
            k2 = self._compute_dS(S_k2, d_eps, dt)
            
            S_k3 = [S[i] + 0.5 * h * k2[i] for i in range(self.n_layer)]
            k3 = self._compute_dS(S_k3, d_eps, dt)
            
            S_k4 = [S[i] + h * k3[i] for i in range(self.n_layer)]
            k4 = self._compute_dS(S_k4, d_eps, dt)
            
            for i in range(self.n_layer):
                S[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])

        # 记录当前步的各分支试算应力和总应力
        self.Tstress_i = S
        self.Tstress = self.k0 * self.Tstrain + sum(self.Tstress_i)

        # 计算切线刚度 (为确保收敛，此处使用伪隐式后向欧拉的算法刚度公式)
        Ttangent = self.k0
        for i in range(self.n_layer):
            Si = self.Tstress_i[i]
            k_val = self.k_ls[i]
            c_val = self.c_ls[i]
            alpha_val = self.alpha_ls[i]

            if Si == 0.0:
                if alpha_val < 1.0:
                    g = 0.0
                elif alpha_val == 1.0:
                    g = k_val * dt / c_val
                else:
                    g = float('inf')
            else:
                power = 1.0 / alpha_val - 1.0
                g = (k_val * dt / (alpha_val * c_val)) * math.pow(abs(Si) / c_val, power)

            # 累加有效的分支刚度
            if not math.isinf(g):
                Ttangent += k_val / (1.0 + g)

        self.Ttangent = Ttangent

    def commitState(self):
        self.Cstrain = self.Tstrain
        self.Cstress = self.Tstress
        self.Ctangent = self.Ttangent
        self.Cstress_i = list(self.Tstress_i)

    def getStrain(self):
        return self.Tstrain

    def getStress(self):
        return self.Tstress

    def getTangent(self):
        return self.Ttangent
    

