# Tutorials-of-metabolic-modeling

## 关于我
长期从事代谢网络模型相关工作，以一作身份在国际期刊如PNAS（4篇）和MSB等发文10+篇（见[ResearchGate](https://www.researchgate.net/profile/Yu-Chen-104/publications)），并受邀审稿（见[Publons](https://publons.com/researcher/3561775/yu-chen/peer-review/)）。


## 关于你
具有微生物学和生物化学的知识储备，对代谢网络模型感兴趣，打算在未来研究中使用模型。建议提前阅读模型相关的经典文献以掌握一定基础知识，例如模型构建（[_Nat Protoc_ 2010](https://www.nature.com/articles/nprot.2009.203)）、流量平衡分析（[_Nat Biotechnol_ 2010](https://www.nature.com/articles/nbt.1614)）等。


## 关于本教程
将介绍细胞代谢网络模型的构建和使用的实际操作。

### 教程内容
将涉及基础、应用和进阶三部分：
* 基础部分（Basic modeling）主要是模型基本信息介绍以及模型构建；
* 应用部分（Model applications）包括细胞代谢流计算，微生物细胞工厂分析和设计等；
* 进阶部分（Advanced modeling）涉及新型模型，特别是蛋白约束模型。

### 需要用的软件和工具
请自行下载并安装以下软件和工具包：
* MATLAB（本教程的所有操作都是在MATLAB中实现的）
* GitHub（非必须，但如果你将主要从事模型研究强烈建议下载安装并注册账号）
* [COBRA](https://github.com/opencobra/cobratoolbox)（模型构建和分析的工具箱，点击[链接](https://opencobra.github.io/cobratoolbox/stable/installation.html)有详细的安装说明，Solver最好选择Gurobi）

**注意** 安装完毕后，每次使用本教程的时候都需要运行startme文件以加载COBRA工具包，或者自行输入initCobraToolbox进行加载。

## 其它
如果在使用过程中有任何问题，可以在[Issues](https://github.com/Yu-sysbio/Tutorials-of-metabolic-modeling/issues)提出。

如对教程内容有任何建议，请与我联系（yuchen.scholar@gmail.com）。 

* 本教程将长期但不定期更新，建议点击右上角Watch可订阅每次更新。
* 最后更新时间：2022-04-26