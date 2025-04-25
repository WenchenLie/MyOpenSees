import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from ModBoucWen.ModBoucWen import ModBoucWen
from Failure.Failure import Failure
from utils.utils import generate_path
from utils.openseespy_test import generate_path, material_test


u = generate_path([0, 
                   5, -5, 5, -5, 5, -5,
                   10, -10, 10, -10, 10, -10,
                   15, -15, 15, -15, 15, -15,
                   20, -20, 20, -20, 20, -20,
                   25, -25, 25, -25, 25, -25,
                   30, -30, 30, -30, 30, -30,
                   35, -35, 35, -35, 35, -35,
                   40, -40, 40, -40, 40, -40,
                   45, -45, 45, -45, 45, -45,
                   50, -50, 50, -50, 50, -50,
                   55, -55, 55, -55, 55, -55,
                   60, -60, 60, -60, 60, -60,
                   65, -65, 65, -65, 65, -65,
                   70, -70, 70, -70, 70, -70,
                   90, -90, 100, -100, 110, -110,
                   0])
# u = generate_path([0,
#                    70, -70,
#                 #    70, -70,
#                    0,
#                    ])
# u = generate_path([0, 100, 80, 120])
# u = generate_path([0, 60], 10)
data = np.loadtxt('data/1完整.txt')
u, F_exp = data[:, 0], data[:, 1]
# u = np.loadtxt('in.txt')
# paras  = [20, 50, 0, 35, 50, 2, 1, '-hardening', 50, 0.05, 0.02, 0.4, '-minmax', 80, '-configType', 2]
paras  = [20.24, 138.58, 0, 50.67, 64.20, 2.34, 0.3, '-hardening', 41.3, 0.022, 0.01, 0.47, '-minmax', 800, '-configType', 1]
# paras1 = [20, 100, 20, 35, 50, 2, 0.3, '-hardening', 50, 0.0, 0, 0.4, '-minmax', 63]
# paras2 = [20, 100, 20, 35, 50, 2, 0.3, '-hardening', 50, 0.0, 0, 0.4, '-minmax', 53]
# PARA1 = {'F1': 20.24, 'k0': 138.58, 'ugap': 20, 'F2': 50.67, 'k1': 64.20, 'k2': 2.34, 'beta': 0.281, 'uh': 41.3, 'r1': 0.013, 'r2': 0.007, 'r3': 0.47, 'minmax': 1000}  # TSSCB-1
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

mat = TSSCBMaterial(1, *paras)
F = []
for ui in u:
    mat.setStrain(ui)
    Fi = mat.getStress()
    F.append(Fi)

# F, _ = np.array(material_test(u, 'TSSCB', paras))
# F1, _ = np.array(material_test(u, 'TSSCB', paras1))
# F2, _ = np.array(material_test(u, 'TSSCB', paras2))

# plt.plot(u, F_exp)
plt.plot(u, F, label='F')
# plt.plot(u, F1, label='F1')
# plt.plot(u, F2, label='F2')
plt.legend()
plt.show()

np.savetxt('data/data.out', np.array([u, F]).T)
# np.savetxt('data1.out', np.array(F1))
# np.savetxt('data2.out', np.array(F2))
