import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from ModBoucWen.ModBoucWen import ModBoucWen
from Failure.Failure import Failure
from utils.utils import generate_path


u = generate_path([0, 5, -5, 10, -10, 15, -15, 20, -20, 25, -25, 30, -30, 0])
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

#                tag Fy uy  alpha  n   Q    b    A  beta gamma
mat2 = ModBoucWen(2, 50, 2, 0.01,  3,  0.5, 1.005, 1, 0.5, 0.5)
mat3 = Failure(3, 2, uy=2, CPD=100)

F = []
for ui in u:
    mat3.setStrain(ui)
    Fi = mat3.getStress()
    F.append(Fi)


plt.plot(u, F)
plt.show()
