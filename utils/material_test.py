from typing import Type, TypeVar
from src.UniaxialMaterial import UniaxialMaterial


T = TypeVar('T', bound=UniaxialMaterial)

def test_opspy(
    strain: list[float],
    mat_type: str,
    paras_args: tuple,
    paras_kwargs: dict = {},
    strainRate: list[float] = None,
) -> tuple[list, list]:
    """基于项目中openseespy.pyd的材料测试

    Args:
        strain (list[float]): 应变序列
        mat_type (str): 材料名称
        paras_args (tuple): 参数
        paras_kwargs (dict): 参数
        strainRate (list[float], optional): 应变率序列

    Returns:
        tuple[list, list]: 应力、切线刚度
    """
    import bin.opensees as ops
    ops.wipe()
    print(mat_type, 1, *paras_args, *paras_kwargs.values())
    ops.uniaxialMaterial(mat_type, 1, *paras_args, *paras_kwargs.values())
    ops.testUniaxialMaterial(1)

    stess = []
    tangent = []
    for i, val in enumerate(strain):
        # ops.setStrain(val)
        # 只有我编译的oponseespy有setTrailStrain和commitState函数
        # 如果没有这些函数，可以直接使用setStrain
        if strainRate is None:
            ops.setTrialStrain(val)
        else:
            ops.setTrialStrain(val, strainRate[i])
        ops.commitState()
        stess.append(ops.getStress())
        tangent.append(ops.getTangent())
    return stess, tangent


def test_py(
    strain: list[float],
    mat_cls: Type[T],
    paras_args: tuple,
    paras_kwargs: dict,
    strainRate: list[float] = None,
) -> tuple[list[float], list[float]]:
    """基于Python材料类的材料测试

    Args:
        strain (list[float]): 应变序列
        mat_cls (Type[T]): 材料类(UniaxialMaterial的子类)
        paras_args (list): 参数
        paras_kwargs (dict): 参数
        strainRate (list[float], optional): 应变率序列

    Returns:
        tuple[list[float], list[float]]: 应力、切线刚度
    """
    mat = mat_cls(1, *paras_args, **paras_kwargs)
    stress = []
    tangent = []
    for i, val in enumerate(strain):
        if strainRate is None:
            mat.setStrain(val)
        else:
            mat.setStrain(val, strainRate[i])
        stress.append(mat.getStress())
        tangent.append(mat.getTangent())
    return stress, tangent

def test_ext(
    strain: list[float],
    mat_type: str,
    paras_args: tuple,
    paras_kwargs: dict,
    strainRate: list[float] = None,) -> tuple[list[float], list[float]]:
    """导入python的c扩展进行材料测试

    Args:
        strain (list[float]): 应变序列
        mat_type (str): 材料名称
        paras_args (list): 参数
        paras_kwargs (dict): 参数
        strainRate (list[float], optional): 应变率序列

    Returns:
        tuple[list[float], list[float]]: 应力、切线刚度
    """
    import importlib
    module = importlib.import_module(f'ext.{mat_type}')
    mat_cls = getattr(module, mat_type)
    mat = mat_cls(1, *paras_args, **paras_kwargs)
    stress = []
    tangent = []
    for i, val in enumerate(strain):
        if strainRate is None:
            mat.setStrain(val)
        else:
            mat.setStrain(val, strainRate[i])
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

