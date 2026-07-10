import os
import subprocess
from pathlib import Path
from typing import Type, TypeVar
import numpy as np
from src.UniaxialMaterial import UniaxialMaterial


T = TypeVar('T', bound=UniaxialMaterial)

def generate_path(disp_level: list, n: int=200, sf: float=1):
    u = []
    for i, disp in enumerate(disp_level[1:]):
        for j in range(n):
            ui = disp_level[i] + (disp_level[i + 1] - disp_level[i]) * j / n
            u.append(ui * sf)
    else:
        u.append(disp_level[-1] * sf)
    return u

class MaterialTest:
    def __init__(self,
        strain: list[float],
        mat_type: str,
        paras_args: tuple,
        paras_kwargs: dict = {},
        strainRate: list[float] = None,
        mat_tag: int = 1,
        exe_path: Path = Path('bin/OpenSees.exe'),
    ):
        """材料测试类，支持四种测试类型（py、由cython编译的c扩展，openseespy、OpenSees.exe）

        Args:
            strain (list[float]): 应变序列
            mat_type (str): 材料名称
            paras_args (tuple): 参数
            paras_kwargs (dict): 参数
            strainRate (list[float], optional): 应变率序列
            mat_tag (int, optional): 材料序号
        """
        self.strain = strain
        self.mat_type = mat_type
        self.paras_args = paras_args
        self.paras_kwargs = paras_kwargs
        self.strainRate = strainRate
        self.mat_tag = mat_tag
    
    def test_opspy(self) -> tuple[list, list]:
        import bin.opensees as ops
        ops.wipe()
        ops.uniaxialMaterial(self.mat_type, self.mat_tag, *self.paras_args, *self.paras_kwargs.values())
        ops.testUniaxialMaterial(self.mat_tag)

        stess = []
        tangent = []
        for i, val in enumerate(self.strain):
            # ops.setStrain(val)
            # 只有我编译的oponseespy有setTrailStrain和commitState函数
            # 如果没有这些函数，可以直接使用setStrain
            if self.strainRate is None:
                ops.setTrialStrain(val)
            else:
                ops.setTrialStrain(val, self.strainRate[i])
            ops.commitState()
            stess.append(ops.getStress())
            tangent.append(ops.getTangent())
        return stess, tangent

    def test_py(self, mat_cls: Type[T]) -> tuple[list, list]:
        mat = mat_cls(self.mat_tag, *self.paras_args, **self.paras_kwargs)
        stress = []
        tangent = []
        for i, val in enumerate(self.strain):
            if self.strainRate is None:
                mat.setStrain(val)
            else:
                mat.setStrain(val, self.strainRate[i])
            stress.append(mat.getStress())
            tangent.append(mat.getTangent())
        return stress, tangent

    def test_ext(self) -> tuple[list, list]:
        import importlib
        module = importlib.import_module(f'ext.{self.mat_type}')
        mat_cls = getattr(module, self.mat_type)
        mat = mat_cls(self.mat_tag, *self.paras_args, **self.paras_kwargs)
        stress = []
        tangent = []
        for i, val in enumerate(self.strain):
            if self.strainRate is None:
                mat.setStrain(val)
            else:
                mat.setStrain(val, self.strainRate[i])
            stress.append(mat.getStress())
            tangent.append(mat.getTangent())
        return stress, tangent

    def test_ops(self, exe_path: Path = Path('bin/OpenSees.exe')) -> tuple[list, list]:
        folder = Path('temp')
        if not folder.exists():
            os.makedirs(folder)
        tcl_file = folder / "run_material_test.tcl"
        reaction_file = folder / "reaction.out"
        tangent_file = folder / "tangent.out"
        for f in [reaction_file, tangent_file]:
            if f.exists():
                os.remove(f)
        
        if not exe_path.exists():
            print(f"找不到指定的执行文件：{exe_path}，请检查路径是否正确。")
            return None, None

        mat_command = f"uniaxialMaterial {self.mat_type} {self.mat_tag} {' '.join(map(str, self.paras_args))}"
        if self.paras_kwargs:
            mat_command += f" {' '.join([f'{k} {v}' for k, v in self.paras_kwargs.items()])}"

        tcl_script = f"""wipe
model testUniaxial
{mat_command}
uniaxialTest {self.mat_tag}
set file_reaction [open {reaction_file.as_posix()} "w"]
set file_tangent [open {tangent_file.as_posix()} "w"]
"""

        for target_disp in self.strain:
            tcl_script += f"strainUniaxialTest {target_disp}\n"
            tcl_script += "set stress [stressUniaxialTest]\n"
            tcl_script += "set tagent [tangUniaxialTest]\n"
            tcl_script += "puts $file_reaction $stress\n"
            tcl_script += "puts $file_tangent $tagent\n"
        tcl_script += "wipe\n"

        with open(tcl_file, "w", encoding="utf-8") as f:
            f.write(tcl_script)

        try:
            result = subprocess.run([exe_path, tcl_file], 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    text=True)
            if result.returncode != 0:
                print("OpenSees 运行出错：")
                print(result.stderr)
                return None, None
        except FileNotFoundError:
            print(f"找不到指定的执行文件：{exe_path}，请检查路径是否正确。")
            return None, None

        try:
            stress = np.loadtxt(reaction_file)
            tangent = np.loadtxt(tangent_file)
            return stress, tangent
        except Exception as e:
            print(f"读取输出文件时发生错误: {e}")
            return None, None
