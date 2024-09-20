import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from ModBoucWen.ModBoucWen import ModBoucWen
from Failure.Failure import Failure
from utils.utils import generate_path
from utils.openseespy_test import generate_path, material_test


# u = generate_path([0, 5, -5, 10, -10, 15, -15, 20, -20, 25, -25, 30, -30, 40, -40, 50, -50, 60, -60, 70, -70, 0])
# u = generate_path([0, 60], 10)
data = np.loadtxt('data/1完整.txt')
u, F_exp = data[:, 0], data[:, 1]
PARA1 = {'F1': 20.24, 'k0': 138.58, 'ugap': 20, 'F2': 50.67, 'k1': 64.20, 'k2': 2.34, 'beta': 0.281, 'uh': 41.3, 'r1': 0.013, 'r2': 0.007, 'r3': 0.47, 'minmax': 1000}  # TSSCB-1
PARA2 = {'F1': 19.29, 'k0': 145.04, 'ugap': 10, 'F2': 50.61, 'k1': 60.89, 'k2': 2.32, 'beta': 0.278, 'uh': None, 'r1': None, 'r2': None, 'r3': None, 'minmax': None}  # TSSCB-2
PARA3 = {'F1': 20.00, 'k0': 65.18,  'ugap': 0,  'F2': 47.97, 'k1': 65.18, 'k2': 2.23, 'beta': 0.283, 'uh': None, 'r1': None, 'r2': None, 'r3': None, 'minmax': None}  # TSSCB-3
PARA6 = {'F1': 21.37, 'k0': 130.55, 'ugap': 20, 'F2': 28.89, 'k1': 43.27, 'k2': 1.33, 'beta': 0.456, 'uh': None, 'r1': None, 'r2': None, 'r3': None, 'minmax': None}  # TSSCB-6
PARA7 = {'F1': 0,     'k0': 28.16,  'ugap': 0,  'F2': 44.03, 'k1': 28.16, 'k2': 2.32, 'beta': 0.095, 'uh': None, 'r1': None, 'r2': None, 'r3': None, 'minmax': None}  # TSSCB-7

PARA = PARA1
paras = (
    PARA['F1'],
    PARA['k0'],
    PARA['ugap'],
    PARA['F2'],
    PARA['k1'],
    PARA['k2'],
    PARA['beta'],
    '-hardening',
    PARA['uh'],
    PARA['r1'],
    PARA['r2'],
    PARA['r3'],
    '-minmax',
    PARA['minmax']
)

# mat = TSSCBMaterial(1, *paras)
# F = []
# for ui in u:
#     mat.setStrain(ui)
#     Fi = mat.getStress()
#     F.append(Fi)

F, _ = material_test(u, 'TSSCB', paras)

plt.plot(u, F_exp)
plt.plot(u, F)
plt.show()
