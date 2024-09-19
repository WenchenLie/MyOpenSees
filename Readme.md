# OpenSees个人二次开发
## 1 改进BoucWen模型
```tcl
ModBoucWen $tag $Fy $uy $alpha $alpha $n $Q $b $A $beta $gamma <$iter>
```
*tag: 材料编号  
*Fy: 屈服力  
*uy: 屈服位移  
*alpha: 屈服后刚度比  
*n: 弹塑性过渡段平滑系数  
*Q: 循环硬化系数(控制最大硬化程度, Q >= 0)  
*b: 循环硬化系数(控制硬化速度, b >= 1)  
*A: 滞回环形状参数(建议取1)  
*beta: 滞回环形状参数(建议取0.5)  
*gamma: 滞回环形状参数(建议取0.5)  
*iter: 迭代次数(默认10)  
## 2 双阶自复位模型（TSSCB）
```tcl
TSSCB $tag $F1 $k0 $ugap $F2 $k1 $k2 $beta <-hardening $uh $r1 $r2 $r3> <-minmax $uf>
```
*tag: 材料编号  
*F1: 摩擦滑移力  
*k0: 第一阶段初始刚度  
*ugap: 间隙大小  
*F2: 自复位模型(第1组)的伪屈服力  
*k1: 自复位模型(第1组)的初始刚度  
*k2: 自复位模型(第1组)的屈服后刚度  
*beta: 自复位模型(第1组)的耗能系数  
*uh: 硬化开始位移
*r1: 控制第二阶段的起始阶段的承载力退化
*r2: 控制第二阶段的后续承载力退化
*r3: 控制硬化位移后承载力硬化的程度
*uf: SMA线缆断裂位移
