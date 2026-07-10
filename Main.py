import numpy as np
import matplotlib.pyplot as plt

from src.ModBoucWen.ModBoucWen import ModBoucWen
from src.TwoStage.TwoStage import TwoStage
from src.Failure.Failure import Failure
from utils.material_test import generate_path, MaterialTest


u = generate_path([0, 
                   10, -10, 10, -10, 10, -10,
                   20, -20, 20, -20, 20, -20,
                   30, -30, 30, -30, 30, -30,
                   40, -40, 40, -40, 40, -40,
                   50, -50, 50, -50, 50, -50,
                   0], n=100, sf=1)
paras = [20, 20, 0.05, 29, 20, 0.2, 30]
tester = MaterialTest(u, 'TwoStage', paras)
F_py, _ = tester.test_py(TwoStage)
F_opspy, _ = tester.test_opspy()
F_ext, _ = tester.test_ext()
F_ops, _ = tester.test_ops()

plt.plot(u, F_py, label='F_py')
plt.plot(u, F_opspy, label='F_opspy')
plt.plot(u, F_ext, label='F_ext')
plt.plot(u, F_ops, '--', label='F_ops')
plt.legend()
plt.show()
