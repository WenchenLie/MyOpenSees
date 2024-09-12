import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from ModBoucWen.ModBoucWen import ModBoucWen
from Failure.Failure import Failure
from utils.utils import generate_path


# u = generate_path([0, 5, -5, 10, -10, 15, -15, 20, -20, 25, -25, 30, -30, 0])
u = generate_path([0, 60], 10)
PARA3 = {'F1': 20000, 'k0': 20000, 'ugap': 20, 'F2': 40000, 'k1': 1000, 'k2': 20, 'beta': 0.2}

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


plt.plot(u, F, '-o')
plt.show()
