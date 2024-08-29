import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from utils.utils import generate_path


u = generate_path([0, 5, -5, 10, -10, 20, -20, 30, -30, 0])
PARA3 = {'F1': 100, 'k0': 300, 'ugap': 10, 'F2': 300, 'k1': 120, 'k2': 10, 'beta': 0.4}

mat = TSSCBMaterial(
    1,
    PARA3['F1'],
    PARA3['k0'],
    PARA3['ugap'],
    PARA3['F2'],
    PARA3['k1'],
    PARA3['k2'],
    PARA3['beta'],
)

F = []
for ui in u:
    mat.setStrain(ui)
    Fi = mat.getStress()
    F.append(Fi)


plt.plot(u, F)
plt.show()
