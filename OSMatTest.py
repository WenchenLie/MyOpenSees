import numpy as np
import openseespy.opensees as ops
import matplotlib.pyplot as plt


def material_test(
        strain: list[float],
        mat_type: str,
        para: list,
    ) -> list[float]:
    ops.wipe()
    ops.uniaxialMaterial(mat_type, 1, *para)
    ops.testUniaxialMaterial(1)
    stess = []
    for val in strain:
        ops.setStrain(val)
        stess.append(ops.getStress())
    return stess


def generate_path(disp_level: list, n: int=200):
    u = []
    for i, disp in enumerate(disp_level[1:]):
        for j in range(n):
            ui = disp_level[i] + (disp_level[i + 1] - disp_level[i]) * j / n
            u.append(ui)
    else:
        u.append(disp_level[-1])
    return u




if __name__ == "__main__":

    # u = np.loadtxt('data/u.txt')
    u = np.array([0, 10, -10, 20, -20, 30, -30, 40, -40, 55, -55, 0]) / 10
    # u = np.array([0, 55]) / 10
    u = generate_path(u)  # 目标位移幅值点


    mat = 'TSSCB'
    paras = (48.085369617556054, 2526.6187266788756, 0.190315100216019, 240.42684808778026, 57.63129379564053, 2.074726576643059, 0.38270000000000004)

    F = material_test(u, mat, paras)
    F = np.array(F)
    plt.plot(u, F)
    plt.show()
