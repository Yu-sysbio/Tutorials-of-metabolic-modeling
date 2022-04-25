# Flux estimation

代谢网络模型可以用来计算胞内代谢流分布，有助于定量理解细胞代谢状态。可用于比较同一菌株在不同条件下的代谢差异，也可以比较不同菌株之间的差异。

这里将结合[_Biotechnol Bioeng_ 2020](https://onlinelibrary.wiley.com/doi/10.1002/bit.27488)的实例，介绍如何利用基于约束的方法计算胞内代谢流量。

## 概念介绍
基于约束的方法就是将代谢网络模型转换为线性方程组，利用线性规划的方法，寻找满足给定约束的解，这里的解即为胞内代谢流量。值得注意的是，由于模型中包含的约束有限，通常无法确定唯一解。为解决这个问题，一般可以使用两类方法：
1. biased方法

代表性的是FBA（Flux Balance Analysis）。该方法是在给定的约束下，还需要人为确定一个目标函数，一方面既能找到目标函数对应的最值，另一方面也能给出在满足约束且达到最值时胞内的流量分布。然而此时给出的结果也是从众多可能的流量分布中挑选一组，因此存在不确定性。

2. unbiased方法

例如FVA（Flux Variability Analysis）和Sampling。这些方法不会给定一组明确的流量分布，而是给出在当前约束下所有的可能性，因此有助于消除不确定性，可与FBA方法互补。

建议阅读以下参考文献：

基于约束的方法[_Nat Rev Microbiol_ 2012](https://www.nature.com/articles/nrmicro2737)，FBA方法[_Nat Biotechnol_ 2010](https://www.nature.com/articles/nbt.1614)，pFBA方法[_Mol Syst Biol_ 2010](https://www.embopress.org/doi/full/10.1038/msb.2010.47)，FVA方法[_Metab Eng_ 2003](https://www.sciencedirect.com/science/article/pii/S1096717603000582)，Sampling方法[_J Biol Chem_ 2009](https://doi.org/10.1074/jbc.R800048200)