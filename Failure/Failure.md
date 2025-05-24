# 失效准则（Failure）
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
