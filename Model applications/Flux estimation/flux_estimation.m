%% flux estimation
% 比较凝结芽孢杆菌（Bacillus coagulans）在不同pH条件下的胞内代谢流量差异。
% 背景介绍：同一株凝结芽孢杆菌分别在pH为6和pH为5.5的条件下进行恒化培养，稀释速率为0.2 /h

% 注：为简化教程，模型以及模拟方法与原始文献中略有差异，因此结果也有所差别。

%% 读取模型
load('iBag597.mat'); % 得到的模型文件为model
% ".mat"文件可直接使用MATLAB自带的function读取。

%% 设置模型
% 因为要比较两个条件，所以需要构建两个模型，分别是pH6和pH5.5的模型。
% 这两个模型由于都是针对同一菌株，所以在结构上没有差别。差别在于测量到的交换流量，因此需要用到条件特异性的数据来约束模型。

% 实验数据可以自行写入到Excel等类型的文件中，再导入到MATLAB。
[~, ~, data_bc] = xlsread('data_bc_chemostat.xlsx'); % 导入已经在Excel中提前备好的实验数据
exchange_rxns = data_bc(2:end,1); % 读取需要修改流量的反应ID
lb_ph6 = cell2mat(data_bc(2:end,ismember(data_bc(1,:),'LB_ph6'))); % 读取pH6条件下对应反应的流量的下限
ub_ph6 = cell2mat(data_bc(2:end,ismember(data_bc(1,:),'UB_ph6'))); % 读取pH6条件下对应反应的流量的上限
lb_ph55 = cell2mat(data_bc(2:end,ismember(data_bc(1,:),'LB_ph55'))); % 读取pH5.5条件下对应反应的流量的下限
ub_ph55 = cell2mat(data_bc(2:end,ismember(data_bc(1,:),'UB_ph55'))); % 读取pH5.5条件下对应反应的流量的上限

% pH6的模型，命名为model_ph6
model_ph6 = model; % 将出发模型赋值给新的模型model_high
model_ph6 = changeRxnBounds(model_ph6,exchange_rxns,lb_ph6,'l'); % 批量设置模型中对应反应的实验测得的流量
model_ph6 = changeRxnBounds(model_ph6,exchange_rxns,ub_ph6,'u');

% pH5.5的模型，命名为model_ph55
model_ph55 = model; % 将出发模型赋值给新的模型model_high
model_ph55 = changeRxnBounds(model_ph55,exchange_rxns,lb_ph55,'l'); % 批量设置模型中对应反应的实验测得的流量
model_ph55 = changeRxnBounds(model_ph55,exchange_rxns,ub_ph55,'u');

%% FBA计算代谢流量
% FBA需要认为确定一个目标函数，这里定为最大化非生长偶联的ATP消耗（NGAM），也可以理解为最大化ATP生成。
model_ph6 = changeObjective(model_ph6,'NGAM');
model_ph55 = changeObjective(model_ph55,'NGAM');
% 求解
sol_ph6 = optimizeCbModel(model_ph6,'max','one');
sol_ph55 = optimizeCbModel(model_ph55,'max','one');
% optimizeCbModel用来求解模型，其中"max"代表对模型中的目标函数进行最大化
% 通常而言传统FBA只能从众多解中可能毫无根据地选择一组解
% 而使用pFBA可以在所有解中选择总流量最小的一组解，意味着细胞会尽可能减少酶的使用以维持给定的状态，因而更具生物学意义
% optimizeCbModel中的"one"即代表在所有解中选择总流量最小的，是非常常用的一种方式
fba_ph6 = sol_ph6.x;
fba_ph55 = sol_ph55.x;

%% FVA计算每个反应的最大和最小值
% FVA能给出在满足当前约束时每个反应的最大值和最小值
[minFlux_ph6, maxFlux_ph6, ~] = fastFVA(model_ph6);
[minFlux_ph55, maxFlux_ph55, ~] = fastFVA(model_ph55);

%% Sampling获得10000组流量
% Sampling即是在满足当前约束下，穷举所有可能的流量结果。取样次数越多，越能看出每个反应的速率可能性分布

% 将FBA计算得到的目标函数的最值作为约束添加到模型中
modelSampling_ph6 = changeRxnBounds(model_ph6,'NGAM',sol_ph6.f*0.99,'b');
modelSampling_ph55 = changeRxnBounds(model_ph55,'NGAM',sol_ph55.f*0.99,'b');
% 有时需要将计算得到的最值稍微调小，比如乘以0.99，使得模型有解

% 设置Sampling的相关参数
options.nStepsPerPoint = 200;
options.nPointsReturned = 10000;

% Sampling计算10000组流量
[~,samples_ph6] = sampleCbModel(modelSampling_ph6, [], [], options);
[~,samples_ph55] = sampleCbModel(modelSampling_ph55, [], [], options);

%% 可视化
% 选取两个反应的FBA，FVA和Sampling结果进行比较

% 定位所关注的三个反应在模型中的位置
xpk_idx = ismember(model.rxns,'R0032'); % phosphoketolase
fba_idx = ismember(model.rxns,'R0008'); % fructose‐bisphosphate aldolase

clr_6 = [5,113,176]/255; % 给pH6的结果特定颜色
clr_55 = [202,0,32]/255; % 给pH5.5的结果特定颜色

figure();
subplot(2,2,1);
hold on;
errorbar(1,fba_ph55(xpk_idx),abs(fba_ph55(xpk_idx)-minFlux_ph55(xpk_idx)),abs(fba_ph55(xpk_idx)-maxFlux_ph55(xpk_idx)),'o','Color',clr_55,'linestyle','none','MarkerSize',10,'LineWidth',1,'CapSize',20);
errorbar(2,fba_ph6(xpk_idx),abs(fba_ph6(xpk_idx)-minFlux_ph6(xpk_idx)),abs(fba_ph6(xpk_idx)-maxFlux_ph6(xpk_idx)),'o','Color',clr_6,'linestyle','none','MarkerSize',10,'LineWidth',1,'CapSize',20);
title('phosphoketolase');
xlim([0 3]);ylim([0 12]);ylabel('Flux (mmol/gCDW/h)');xticks([1 2]);xticklabels({'pH5.5','pH6'});

subplot(2,2,2);
hold on;
errorbar(1,fba_ph55(fba_idx),abs(fba_ph55(fba_idx)-minFlux_ph55(fba_idx)),abs(fba_ph55(fba_idx)-maxFlux_ph55(fba_idx)),'o','Color',clr_55,'linestyle','none','MarkerSize',10,'LineWidth',1,'CapSize',20);
errorbar(2,fba_ph6(fba_idx),abs(fba_ph6(fba_idx)-minFlux_ph6(fba_idx)),abs(fba_ph6(fba_idx)-maxFlux_ph6(fba_idx)),'o','Color',clr_6,'linestyle','none','MarkerSize',10,'LineWidth',1,'CapSize',20);
title('fructose‐bisphosphate aldolase');
xlim([0 3]);ylim([0 12]);ylabel('Flux (mmol/gCDW/h)');xticks([1 2]);xticklabels({'pH5.5','pH6'});

subplot(2,2,3);
hold on;
[y1, x1] = hist(samples_ph55(xpk_idx, :), 100);
[y2, x2] = hist(samples_ph6(xpk_idx, :), 100);
plot(x1, y1, 'Color', clr_55);
plot(x2, y2, 'Color', clr_6);
f1 = fill([x1,fliplr(x1)],[y1,zeros(1,length(x1))],clr_55);
set(f1,'edgealpha',0,'facealpha',0.3);
f2 = fill([x2,fliplr(x2)],[y2,zeros(1,length(x2))],clr_6);
set(f2,'edgealpha',0,'facealpha',0.3);
title('phosphoketolase');
xlim([0 12]);xlabel('Flux (mmol/gCDW/h)');ylabel('Number of Samples');

subplot(2,2,4);
hold on;
[y1, x1] = hist(samples_ph55(fba_idx, :), 100);
[y2, x2] = hist(samples_ph6(fba_idx, :), 100);
plot(x1, y1, 'Color', clr_55);
plot(x2, y2, 'Color', clr_6);
f1 = fill([x1,fliplr(x1)],[y1,zeros(1,length(x1))],clr_55);
set(f1,'edgealpha',0,'facealpha',0.3);
f2 = fill([x2,fliplr(x2)],[y2,zeros(1,length(x2))],clr_6);
set(f2,'edgealpha',0,'facealpha',0.3);
title('fructose‐bisphosphate aldolase');
xlim([0 12]);xlabel('Flux (mmol/gCDW/h)');ylabel('Number of Samples');


