from typing import Type, TypeVar
from src.UniaxialMaterial import UniaxialMaterial


T = TypeVar('T', bound=UniaxialMaterial)

def test_opspy(
    strain: list[float],
    mat_type: str,
    para: list,
) -> tuple[list, list]:
    """基于openseespy的材料测试

    Args:
        strain (list[float]): 应变序列
        mat_type (str): 材料名称
        para (list): 参数

    Returns:
        tuple[list, list]: 应力、切线刚度
    """
    import bin.opensees as ops
    ops.wipe()
    ops.uniaxialMaterial(mat_type, 1, *para)
    ops.testUniaxialMaterial(1)

    stess = []
    tangent = []
    for val in strain:
        # ops.setStrain(val)
        # 只有我编译的oponseespy有setTrailStrain和commitState函数
        # 如果没有这些函数，可以直接使用setStrain
        ops.setTrialStrain(val)
        ops.commitState()
        stess.append(ops.getStress())
        tangent.append(ops.getTangent())
    return stess, tangent


def test_py(
    strain: list[float],
    mat_cls: Type[T],
    para: list,
) -> tuple[list[float], list[float]]:
    """基于Python的材料测试

    Args:
        strain (list[float]): 应变序列
        mat_cls (Type[T]): 材料类(UniaxialMaterial的子类)
        para (list): 参数

    Returns:
        tuple[list[float], list[float]]: 应力、切线刚度
    """
    mat = mat_cls(1, *para)
    stress = []
    tangent = []
    for val in strain:
        mat.setStrain(val)
        stress.append(mat.getStress())
        tangent.append(mat.getTangent())
    return stress, tangent



def generate_path(disp_level: list, n: int=200, sf: float=1):
    u = []
    for i, disp in enumerate(disp_level[1:]):
        for j in range(n):
            ui = disp_level[i] + (disp_level[i + 1] - disp_level[i]) * j / n
            u.append(ui * sf)
    else:
        u.append(disp_level[-1] * sf)
    return u

