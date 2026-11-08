# 双阶自复位模型（TSSCB）
双阶自复位模型（TSSCB）是一种在传统旗帜型自复位模型的基础上改进的模型，用于模拟双阶自复位支撑的力学行为。
在小位移时模型处于第一阶段，并表现出与理性弹塑性模型一致的滞回响应(通过参数`F1`和`k0`控制)，
在位移超出阈值(`ugap`)时，进入第二阶段，并表现出旗帜型滞回环(通过参数`F2`, `k1`, `k2`, 和`beta`控制)。
该模型可考虑SMA材料的硬化特性，通过参数`uh`, `r1`, `r2`, `r3`控制，从而模拟马氏体硬化导致的强度进一步
增大和残余变形的增加。
此外，参数`uf`为极限位移，当位移绝对值超过该值时认为SMA材料的失效。
SMA材料失效后，由于耗能元件还可继续变形并耗能，因为模型的失效后行为也应纳入考虑。
通过参数`congifType`控制两种不同构造的支撑的失效后行为。
```tcl
TSSCB $tag $F1 $k0 $ugap $F2 $k1 $k2 $beta <-hardening $uh $r1 $r2 $r3> <-minmax $uf> <-configType $configType> <-up $up>
```
* tag: 材料编号  
* F1: 摩擦滑移力  
* k0: 第一阶段初始刚度  
* ugap: 间隙大小  
* F2: 自复位模型的伪屈服力  
* k1: 自复位模型的初始刚度  
* k2: 自复位模型的屈服后刚度  
* beta: 自复位模型的耗能系数(0 <= beta <= 1)  
* uh: 硬化开始位移  
* r1: 控制第二阶段的起始阶段的承载力退化  
* r2: 控制第二阶段的后续承载力退化  
* r3: 控制硬化位移后承载力硬化的程度  
* uf: SMA线缆断裂位移
* configType: 构造类型参数(默认为1)。1: 第二阶段中仅一半数量的耗能元件参与工作；2: 任意时刻中所有耗能元件均参与工作。
* up: SMA线缆开始断裂至完全断裂间的可滑动距离，避免断裂后承载力的瞬间下降

SMA 完全断裂后的摩擦力上限：`configType=1` 为 `F1/2`，`configType=2` 为 `F1`。
当 `up>0` 时，断裂过渡段也以相应摩擦力为终点；构造 1 断裂后的摩擦屈服位移为 `F1/(2*k0)`。

2026-11-08 同步修复 Python、C++、Cython 的 gap 阶段衔接：从边界实际内力连续进入下一阶段，
边界退化力在 gap 处计算，并返回一致的分段切线。`TSSCB.c` 是 `.pyx` 重新生成的 C 源码，
不是已编译扩展；更新源码后仍须自行重建 `.pyd`/OpenSees 才能用于原生接口。

纯 Python 回归检查位于 `checks/tsscb_gap_regression.py`；材料级滞回检查图由
`utils/tsscb_hysteresis_review.py` 生成，默认保存至 `data/tsscb_review_20261108/`，
包括 PDF/SVG/PNG、完整加载与响应 CSV、参数和源码哈希。这些曲线使用 Python 源码，
不是旧 `.pyd` 的结果，也不是地震时程响应或实验观测。

注：当`ugap`取为0时（即无间隙），建议将`F1`取为装置的实际摩擦力大小（影响SMA线缆断裂后行为），`k0`取为与`k1`相同值。  

## 参考文献：  

[1] [Lie W, Shi F, Zhou Y, Yuan W, Erbolat A, Fang D. Development of a two-stage self-centering brace with dissipating-centering sequence for potential multi-objective seismic design. Eng Struct. 2025;336:120464.](https://doi.org/10.1016/j.engstruct.2025.120464)
