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

    # data = np.loadtxt('SMA6.txt')
    u = generate_path([0, 10, -10, 20, -20, 30, -30, 40, -40, 55, -55, 0])  # 目标位移幅值点


    mat = 'TSSCB'
    paras_a = [200, 500, 20, 400, 400, 5, 0.3]

    F = material_test(u, mat, paras_a)
    F = np.array(F)
    plt.plot(u, F)
    plt.show()
