# 广义非线性Maxwell模型（GeneralizedMaxwell）
广义非线性Maxwell模型包括一个弹性弹簧和任意数量（至少1组）的非线性Maxwell模型，采用四阶龙格库塔计算力响应，材料定义格式如下：
```tcl
GeneralizedMaxwell $tag $k0 $k1 $c1 $alpha1 <k2 $c2 $alpha2 ...> <-iter $n_iter>
```
* tag: 材料编号  
* k0: 弹性弹簧刚度  
* ki: 第i个Maxwell模型的弹簧刚度  
* ci: 第i个Maxwell模型的黏壶阻尼系数  
* alphai: 第i个Maxwell模型的黏壶阻尼指数  
* n_iter: 每个增量步的迭代次数，默认10  

注：i ≥ 1，至少有1组Maxwell模型。