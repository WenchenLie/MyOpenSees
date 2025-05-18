import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from ModBoucWen.ModBoucWen import ModBoucWen
from Failure.Failure import Failure
from utils.openseespy_test import generate_path, material_test


# u = generate_path([0, 
                #    5, -5, 5, -5, 5, -5,
#                    10, -10, 10, -10, 10, -10,
                #    15, -15, 15, -15, 15, -15,
#                    20, -20, 20, -20, 20, -20,
#                    25, -25, 25, -25, 25, -25,
#                    30, -30, 30, -30, 30, -30,
#                    35, -35, 35, -35, 35, -35,
#                    40, -40, 40, -40, 40, -40,
#                    45, -45, 45, -45, 45, -45,
#                    50, -50, 50, -50, 50, -50,
#                    55, -55, 55, -55, 55, -55,
#                    60, -60, 60, -60, 60, -60,
#                    65, -65, 65, -65, 65, -65,
#                    70, -70, 70, -70, 70, -70,
#                    90, -90, 100, -100, 110, -110,
                #    0], n=20, sf=2)
# u = generate_path([0,
#                    250, -250,
#                    270, -270,
#                    260, -260,
#                    0,
#                    ])
# u = generate_path([0, 100, 80, 120])
# u = generate_path([0, 60], 10)
data = np.loadtxt('data/data1.txt')
u, F_exp = data[:, 0], data[:, 1]
# u = np.loadtxt('in.txt')
# cable_length = 2000.0  # Cable length in the numerical model
# cable_length_test = 450.0  # Cable length in the experimental test
# F1 = 20e3
# k0 = cable_length_test / cable_length * 30e3
# ugap = 0.0
# F2 = 47.97e3
# k1 = cable_length_test / cable_length * 65.18e3
# k2 = cable_length_test / cable_length * 2.23e3
# beta = 0.566
# uh = cable_length / cable_length_test * 32.0564
# r1 = 0.0293
# r2 = 0.0112
# r3 = 0.8037
# uf = cable_length / cable_length_test * 51.6
# SCB_paras = [F1, k0, ugap, F2, k1, k2, beta, '-hardening', uh, r1, r2, r3, '-minmax', uf, '-configType', 1]
# paras  = [20, 50, 0, 35, 50, 2, 1, '-hardening', 50, 0.05, 0.02, 0.4, '-minmax', 80, '-configType', 2]
# paras  = [20.24, 138.58, 0, 50.67, 64.20, 2.34, 0.3, '-hardening', 41.3, 0.022, 0.01, 0.47, '-minmax', 800, '-configType', 1]
# paras1 = [20, 100, 20, 35, 50, 2, 0.3, '-hardening', 50, 0.0, 0, 0.4, '-minmax', 63]
# paras3 = [20.94, 68.24, 0, 50.22, 68.24, 2.33, 0.566, '-hardening', 27.4807, 0.0164, 0.0011, 0.7771, '-minmax', 51.6, '-configType', 1]
# PARA1 = {'F1': 20.24, 'k0': 138.58, 'ugap': 20, 'F2': 50.67, 'k1': 64.20, 'k2': 2.34, 'beta': 0.281, 'uh': 41.3, 'r1': 0.013, 'r2': 0.007, 'r3': 0.47, 'minmax': 1000}  # TSSCB-1
paras_test = [772674.8082799409, 986960.4401089358, 7.828832614554003, 1931687.0206998521, 789568.3520871487, 986.9604401089358, 0.5]
# mat1 = [20.24, 138.58, 20, 50.67, 64.20, 2.34, 0.281, '-hardening', 50.45, 0.0306, 0.0157, 0.5586]
# mat2 = [19.29, 145.04, 10, 50.61, 60.89, 2.32, 0.278, '-hardening', 40.88, 0.0289, 0.0113, 0.6707, '-minmax', 63.01]
# mat3 = [20.00, 65.18,  0,  47.97, 65.18, 2.23, 0.283, '-hardening', 27.48, 0.0148, 0,      0.9537, '-minmax', 51.48]
# mat7 = [0,     28.16,  0,  44.03, 28.16, 2.32, 0.095, '-hardening', 28.07, 0.0197, 0.0074, 0.5349, '-minmax', 48.68]
# PARA = PARA1
# paras = (
#     PARA['F1'],
#     PARA['k0'],
#     PARA['ugap'],
#     PARA['F2'],
#     PARA['k1'],
#     PARA['k2'],
#     PARA['beta'],
#     '-hardening',
#     PARA['uh'],
#     PARA['r1'],
#     PARA['r2'],
#     PARA['r3'],
#     '-minmax',
#     PARA['minmax']
# )

mat = TSSCBMaterial(1, *paras_test)  # 用python代码计算本构
F_py = []
for ui in u:
    mat.setStrain(ui)
    Fi = mat.getStress()
    F_py.append(Fi)

F_opspy, _ = np.array(material_test(u, 'TSSCB', paras_test))  # 用openseespy计算本构
# F1, _ = np.array(material_test(u, 'TSSCB', paras1))
# F2, _ = np.array(material_test(u, 'TSSCB', paras2))

plt.plot(u, F_py, label='F_py')
plt.plot(u, F_opspy, label='F_opspy')
# plt.plot(u, F1, label='F1')
# plt.plot(u, F2, label='F2')
plt.legend()
plt.show()
# np.savetxt('data/data.out', np.array([u, F]).T)
# np.savetxt('data1.out', np.array(F1))
# np.savetxt('data2.out', np.array(F2))
