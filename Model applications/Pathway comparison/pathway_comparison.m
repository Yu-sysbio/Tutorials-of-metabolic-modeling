%% pathway comparison
% 在酿酒酵母模型中添加外源反应使其能够更高效地合成Spermidine，并比较外源途径和原生途径的得率差异。

%% 读取并设置野生型模型

model_wt = readCbModel('yeast_7.6_COBRA.xml'); % 读取野生型模型并命名为model_wt
% "yeast_7.6_COBRA.xml"是已经公开的酿酒酵母模型，虽然不是最新版，但不影响作为例子进行实操

% 野生型酿酒酵母可以合成Spermidine，为计算基于葡萄糖的最大得率，需要有以下几个准备工作：
% 1. 给定1个单位的葡萄糖吸收速率
model_wt = changeRxnBounds(model_wt,'r_1714',-1,'b');
% % changeRxnBounds可改变模型中反应的速率，"r_1714"是模型中葡萄糖交换反应的ID，"-1"代表的是葡萄糖吸收速率为1 mmol/CDW/h（相应的正数代表分泌），"b"代表同时改变该反应的上限和下限

% 2. 设定培养基组分，即所需要的营养物质都要给足，除了葡萄糖外的其它物质比如氮源、氧气等都要保证能够无限量供应
% 由于模型默认使用最小培养基，即必需的营养物质的交换反应的速率下限都为-1000（意味着吸收速率最大能到1000），所以此步骤可以跳过

% 3. 关闭生长的反应（传统计算最大得率的时候是不考虑生长的）
model_wt = changeRxnBounds(model_wt,'r_2111',0,'b');
% 类似的，"r_2111"是模型中生长反应的ID
% 另外，对于别的模型，也需要关闭模型中默认的非生长偶联的维持能的消耗，yeast_7.6模型中没有该反应因此无需修改

% 4. 最大化Spermidine的分泌速率
model_wt = changeObjective(model_wt,'r_2051',1);
% changeObjective用来设置模型中的目标函数，"1"代表最大化（"-1"代表最小化），"r_2051"是模型中Spermidine交换反应的ID
% 此外要保证模型分泌Spermidine的反应是通畅的
model_wt = changeRxnBounds(model_wt,'r_2051',0,'l'); % 下限是0
model_wt = changeRxnBounds(model_wt,'r_2051',1000,'u'); % 上限是1000


%% 在野生型模型中添加外源反应
% 这里通过引入两套外源途径以构建两个模型（model_1和model_2）：
% 1. 在野生型模型中添加Arginine decarboxylase(EC 4.1.1.19)和Agmatinase(EC3.5.3.11)，并敲除Ornithine decarboxylase
model_1 = model_wt; % 是在野生型模型基础上进行修改，因此可以将野生型模型赋值给新的模型model_1
model_1 = changeRxnBounds(model_1,'r_0817',0,'b'); % 敲除Ornithine decarboxylase，可通过关闭该酶对应的反应（"r_0817"）实现
model_1 = addReaction(model_1,'M1_1','reactionFormula','s_0965[c_03] + s_0794[c_03] -> s_0456[c_03] + agmatine');
% addReaction用来在模型中添加新反应，"M1_1"是赋予的新反应的ID，后面是化学方程式，其中s开头的是模型中原本包含的代谢物，"agmatine"是不含有的所以没有对应的代谢物ID，需要自己添加
model_1 = addReaction(model_1,'M1_2','reactionFormula','agmatine + s_0803[c_03] -> s_1389[c_03] + s_1552[c_03]');
printRxnFormula(model_1,'rxnAbbrList',model_1.rxns(end-1:end),'metNameFlag',1); % 显示所添加的反应

% 2. 在野生型模型中添加Carboxyspermidine dehydrogenase (EC 1.5.1.43)和Carboxyspermidine decarboxylase(EC 4.1.1.-)，并敲除Spermidine synthase
model_2 = model_wt; 
model_2 = changeRxnBounds(model_2,'r_1001',0,'b'); % 敲除Spermidine synthase，可通过关闭该酶对应的反应（"r_1001"）实现
model_2 = addReaction(model_2,'M2_1','reactionFormula','s_1389[c_03] + s_0978[c_03] + s_0794[c_03] + s_1212[c_03] -> Carboxyspermidine + s_1207[c_03] + s_0803[c_03]');
model_2 = addReaction(model_2,'M2_2','reactionFormula','Carboxyspermidine + s_0794[c_03] -> s_1439[c_03] + s_0456[c_03]');
printRxnFormula(model_2,'rxnAbbrList',model_2.rxns(end-1:end),'metNameFlag',1);

% 注意：由于model_1和model_2是在野生型模型基础上进行修改的，所以model_wt中的所有信息都被继承到了新模型中，比如培养基组分、目标函数等，不再需要重复设置

%% 计算各途径基于葡萄糖的最大得率
% 计算基于葡萄糖的最大得率分为两步：
% 1. 最大化产品生成
sol_wt = optimizeCbModel(model_wt,'max'); % 求解已经设置好的野生型模型
% optimizeCbModel用来求解模型，其中"max"代表对模型中的目标函数进行最大化
sol_1 = optimizeCbModel(model_1,'max'); % 求解已经设置好的model_1
sol_2 = optimizeCbModel(model_2,'max'); % 求解已经设置好的model_2

% 2. 找到Spermidine分泌速率和葡萄糖吸收速率，计算得率
% Spermidine分泌速率实际上是求解得到的目标函数的值，即
q_spermidine_wt = sol_wt.f; % sol_wt中包含所有求解的信息，其中f是目标函数的值
q_spermidine_1 = sol_1.f;
q_spermidine_2 = sol_2.f;
% 而葡萄糖吸收速率其实是已经设置好的，即1 mmol/CDW/h
% 因此，三种途径的基于葡萄糖的最大得率如下
Y_wt = q_spermidine_wt/1; % 单位是mol/mol
Y_1 = q_spermidine_1/1;
Y_2 = q_spermidine_2/1;


