%% introduction
% 预计用时：< 1 min
% 以酿酒酵母最新模型Yeast8为例，介绍利用COBRA进行模型读取、调整、求解等操作。

% Author: Yu Chen (yuchen.scholar@gmail.com)

%% 模型读取

% 首先需要下载模型，一般可以在原始文献中找到下载地址 (本例子的模型来自于https://github.com/SysBioChalmers/yeast-GEM)
% 模型文件一般有多种格式，常见的是".xml" ".mat" ".txt" ".xlsx"等
% ".xml"和".mat"是machine-readable的文件
% ".xml"是最为通用的格式，基本适配于大部分模型分析工具包；而".mat"是MATLAB文件，仅适用于基于MATLAB的工具包。
% ".txt"和".xlsx"是human-readable的文件，能够直观了解模型信息，但一般无法作为COBRA的读取文件，特定格式的".xlsx"除外

% 针对".xml"格式的模型，可用下面这一行代码读取
model = readCbModel('yeast-GEM.xml'); % 读取模型并命名为model

% 针对".mat"格式的模型，可用MATLAB自带的function进行读取
% load('yeast-GEM.mat');


%% 模型信息

% 双击工作区Workspace的model可了解模型信息，在introduction.jpg中有具体介绍


%% 模型调整

% 1. 反应速率的调整（以葡萄糖交换反应为例）
%    查看模型中葡萄糖交换反应默认的速率：
printRxnFormula(model,'rxnAbbrList',{'r_1714'},'metNameFlag',1,'printBounds',1); % 模型中葡萄糖交换反应的ID为"r_1714"

%    修改葡萄糖交换反应的速率：
model = changeRxnBounds(model,'r_1714',-10,'b');
    % 以上是在模型中设置葡萄糖吸收速率为10 mmol/gCDW/h
    % changeRxnBounds可用来改变模型中反应的速率
    % "-10"代表的是吸收速率为10 mmol/gCDW/h（相应的正数代表分泌）
    % "b"代表同时改变该反应的上限和下限（若为"u"则改变上限，若为"l"则改变下限）

%    查看模型中葡萄糖交换反应的速率：
printRxnFormula(model,'rxnAbbrList',{'r_1714'},'metNameFlag',1,'printBounds',1); % 模型中葡萄糖交换反应的ID为"r_1714"

%    显示模型中包含的约束
printConstraints(model,-100,100); 

% 2. 目标函数的调整
%    在命令窗口显示模型中的目标函数
printObjective(model); % 可看到目标函数对应的反应ID为r_2111
%    以下显示该反应
printRxnFormula(model,'rxnAbbrList',{'r_2111'},'metNameFlag',1,'printBounds',1);

%    修改模型中的目标函数（以最大化乙醇的生成为例）：
model = changeObjective(model,'r_1761',1);
% changeObjective用来设置模型中的目标函数，"1"代表最大化（"-1"代表最小化），"r_1761"是模型中乙醇交换反应的ID

%    显示模型中的目标函数
printObjective(model); % 可看到目标函数对应的反应已经修改为r_1761
%    以下显示该反应
printRxnFormula(model,'rxnAbbrList',{'r_1761'},'metNameFlag',1,'printBounds',1);


%% 模型求解

% 利用FBA对已经设置好的模型求解
FBAsolution = optimizeCbModel(model); % 将计算结果命名为FBAsolution

% 双击工作区的FBAsolution可了解模型信息
% 其中最为关键的是 "f"（代表目标函数的值）和 "x"（代表每个反应的速率，与模型中的反应ID一一对应）

