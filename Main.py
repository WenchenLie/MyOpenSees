import numpy as np
import matplotlib.pyplot as plt

from TSSCB.TSSCB import TSSCBMaterial
from ModBoucWen.ModBoucWen import ModBoucWen
from Failure.Failure import Failure
from utils.openseespy_test import generate_path, material_test


u = generate_path([0, 
                   10, -10, 10, -10, 10, -10,
                   20, -20, 20, -20, 20, -20,
                   30, -30, 30, -30, 30, -30,
                   40, -40, 40, -40, 40, -40,
                   50, -50, 50, -50, 50, -50,
                   0], n=100, sf=1)
paras = [30, 2, 0.01, 1, 0.5, 1.01, 1, 0.5, 0.5]
mat = ModBoucWen(1, *paras)  # hysteretic response obtained by Python script
F_py = []
for ui in u:
    mat.setStrain(ui)
    Fi = mat.getStress()
    F_py.append(Fi)
F_opspy, _ = np.array(material_test(u, 'ModBoucWen', paras))  # hysteretic responses obtained by OpenSeesPy

plt.plot(u, F_opspy, label='F_opspy')
plt.plot(u, F_py, label='F_py')
plt.legend()
plt.show()
