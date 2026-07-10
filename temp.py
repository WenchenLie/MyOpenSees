import numpy as np
import matplotlib.pyplot as plt

from src.ModBoucWen.ModBoucWen import ModBoucWen
from src.Degradation.Degradation import Degradation
from utils.material_test import generate_path, MaterialTest


u = generate_path([0, 
                   10, -10, 10, -10, 10, -10,
                   20, -20, 20, -20, 20, -20,
                   30, -30, 30, -30, 30, -30,
                   40, -40, 40, -40, 40, -40,
                   50, -50, 50, -50, 50, -50,
                   0], n=100, sf=1)
# u = generate_path([0, 80], n=100)
paras = [20, 1, 0.005, 1, 0, 1, 1, 0.5, 0.5]

mat0 = ModBoucWen(1, *paras)
mat = Degradation(2, 1,
                  uy = 1,
                  u_th = 30,
                  a1 = 0.01,
                  a2 = 0)

F = []
tangent = []
for i, val in enumerate(u):
    mat.setStrain(val)
    F.append(mat.getStress())
    tangent.append(mat.getTangent())
plt.plot(u, F, label='BW')
plt.legend()
plt.show()
