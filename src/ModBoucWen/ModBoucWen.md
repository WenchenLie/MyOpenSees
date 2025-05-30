# 改进BoucWen模型（ModBoucWen）
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

## 参考文献
[1] [Lie W, Wu C, Luo W, et al. Cyclic behaviour of a novel torsional steel-tube damper. J Constr Steel Res. 2022;188:107010.](https://doi.org/10.1016/j.jcsr.2021.107010)
