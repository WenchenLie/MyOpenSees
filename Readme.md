# OpenSees个人二次开发
## 1 改进BoucWen模型（ModBoucWen）
在传统BoucWen模型的基础上增加了循环硬化特性，通过参数`Q`和`b`进行控制，
其中`Q`控制最大硬化的程度`(Q >= 0)`，`b`控制硬化的速率`(b >= 1)`。
```tcl
ModBoucWen $tag $Fy $uy $alpha $n $Q $b $A $beta $gamma <$iter>
```
* tag: 材料编号  
* Fy: 屈服力  
* uy: 屈服位移  
* alpha: 屈服后刚度比  
* n: 弹塑性过渡段平滑系数  
* Q: 循环硬化系数(控制最大硬化程度, Q >= 0)  
* b: 循环硬化系数(控制硬化速度, b >= 1)  
* A: 滞回环形状参数(建议取1)  
* beta: 滞回环形状参数(建议取0.5)  
* gamma: 滞回环形状参数(建议取0.5)  
* iter: 迭代次数(默认10)  
## 2 双阶自复位模型（TSSCB）
双阶自复位模型（TSSCB）是一种在传统旗帜型自复位模型的基础上改进的模型，用于模拟双阶自复位支撑的力学行为。
在小位移时模型处于第一阶段，并表现出与理性弹塑性模型一致的滞回响应(通过参数`F1`和`k0`控制)，
在位移超出阈值(`ugap`)时，进入第二阶段，并表现出旗帜型滞回环(通过参数`F2`, `k1`, `k2`, 和`beta`控制)。
该模型可考虑SMA材料的硬化特性，通过参数`uh`, `r1`, `r2`, `r3`控制，从而模拟马氏体硬化导致的强度进一步
增大和残余变形的增加。
此外，参数`uf`为极限位移，当位移绝对值超过该值时认为SMA材料的失效。
SMA材料失效后，由于耗能元件还可继续变形并耗能，因为模型的失效后行为也应纳入考虑。
通过参数`congifType`控制两种不同构造的支撑的失效后行为。
```tcl
TSSCB $tag $F1 $k0 $ugap $F2 $k1 $k2 $beta <-hardening $uh $r1 $r2 $r3> <-minmax $uf> <-configType $configType>
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
## 3 失效准则（Failure）
一种wrapper材料，对于被引用的材料，提供三种可选的失效准则，分别为: (1) 基于位移的失效，通过参数`minStrain`或
`maxStrain`控制；(2) 基于力的失效，通过参数`minStress`或`maxStress`控制；(3) 基于累积塑性位移的失效，
通过参数`maxCPD`控制(需指定屈服位移`uy`才能计算累积塑性位移)。
```tcl
Failure $tag $otherTag <-minStrain $minStrain> <-maxStrain $maxStrain> <-minStress $minStress> <-maxStress $maxStress> <-maxCPD $uy $maxCPD>
```
* tag: 材料编号  
* otherTag: 引用的材料编号  
* minStrain: 最小位移  
* maxStrain: 最大位移  
* minStress: 最小力  
* maxStress: 最大力  
* uy: 屈服位移  
* maxCPD: 最大累积塑性位移
