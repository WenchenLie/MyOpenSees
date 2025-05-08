import numpy as np
import bin.opensees as ops
import matplotlib.pyplot as plt


def material_test(
        strain: list[float],
        mat_type: str,
        para: list,
    ) -> tuple[list, list]:
    ops.wipe()
    ops.uniaxialMaterial(mat_type, 1, *para)
    ops.testUniaxialMaterial(1)

    stess = []
    tangent = []
    for val in strain:
        ops.setStrain(val)
        stess.append(ops.getStress())
        tangent.append(ops.getTangent())
    return stess, tangent


def generate_path(disp_level: list, n: int=200, sf: float=1):
    u = []
    for i, disp in enumerate(disp_level[1:]):
        for j in range(n):
            ui = disp_level[i] + (disp_level[i + 1] - disp_level[i]) * j / n
            u.append(ui * sf)
    else:
        u.append(disp_level[-1] * sf)
    return u

