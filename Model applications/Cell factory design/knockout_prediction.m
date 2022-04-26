%% knockout prediction
% 预计用时：< 1 min
% 利用OptKnock算法预测提高酿酒酵母生产2,3-丁二醇（2,3-butanediol）的反应敲除靶点
% OptKnock是预测需要敲除的反应而非基因，因此在实验中需要把反应对应的所有isozyme都敲除

% 注：为简化教程，模型以及模拟方法与原始文献中略有差异，因此结果也有所差别。

% 更多案例见COBRA教程(https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialOptKnock.html)

% Author: Yu Chen (yuchen.scholar@gmail.com)

%% 读取并设置模型
model = readCbModel('iMM904.xml'); % 读取酵母模型（不是最新版模型）

% 找到模型中生长的反应和产品的交换反应
biomass = 'biomass_SC5_notrace'; % 模型中生长的反应ID
product = 'EX_btd_RR_e_'; % 模型中2,3-丁二醇的交换反应的ID

% 确定目标函数
model = changeObjective(model,biomass,1); % 设置目标函数为最大化生长

% 在模型中添加合理的约束
model = changeRxnBounds(model,'EX_glc_e_',-10,'l'); % 设置葡萄糖吸收速率最多不超过10 mmol/gCDW/h
model = changeRxnBounds(model,'EX_o2_e_',-2,'l'); % 设置氧气吸收速率最多不超过2 mmol/gCDW/h
model = changeRxnBounds(model,'ATPM',1,'b'); % 设置非生长偶联的ATP消耗为1 mmol/gCDW/h
% 此外还需要确保其它必需的营养物质可自由摄取，iMM904默认是可以的，若使用其它生物的模型需要自行确认

%% 模拟用野生型酿酒酵母生产2,3-丁二醇
solWT = optimizeCbModel(model,'max'); % 求解模型，并将求解结果命名为solWT
productFluxWT = solWT.x(strcmp(model.rxns,product)); % 在得到的解中找到产品生成速率
% ".x"代表求得的代谢流量，需要用反应的ID定位到所关注的反应
growthRateWT = solWT.f; % 在得到的解中找到生长速率
% ".f"代表目标函数的值，而模型的目标函数设置的恰好是生长

fprintf('The production rate of 2,3-butanediol before optimization is %.3f mmol/gCDW/h \n', productFluxWT);
fprintf('The growth rate before optimization is %.3f /h \n', growthRateWT);
% 上面两行代码是将结果显示在命令窗口

%% 设置OptKnock算法所需要的参数

% 最多给出多少组改造策略
threshold = 5; % 最多给出5组改造策略

% 给定每组策略中同时敲除多少个反应
numDelRxns = 3; % 同时敲除3个反应（可从1开始依次增加，这里仅以3为例）

% 给定需要对哪些反应进行敲除预测
selectedRxnList = {'ALCD2ir';'ALCD2irm';'ALCD2x';'ASPTA';'CYTK1';'GLUDyi';'GTPCI';'MDH';'NADK';'NDPK3';'PDHm';'TMDPP'};
% 通常需要对模型中所有的酶催化反应进行搜索，可用下面这一行代码，但搜索的反应越多所需要的时间越长，因此本教程只搜索如上若干反应作为例子。
% selectedRxnList = model.rxns(~ismember(model.rules,''));

options = struct('targetRxn', product, 'numDel', numDelRxns);
constrOpt = struct('rxnList', {{biomass}},'values', 0.5*solWT.f, 'sense', 'G'); %规定生长至少要达到野生型的50%

% 以下即是利用OptKnock算法寻找需要敲除的反应
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 1;
while nIter < threshold
    fprintf('...Performing optKnock analysis...\n')
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt);
    else
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt, previousSolutions, 1);
    end
    productFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, product));
    growthRateM1 = optKnockSol.fluxes(strcmp(model.rxns, biomass));
    setM1 = optKnockSol.rxnList;
    if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions = contPreviousSolutions + 1;
        fprintf('optKnock found a optKnock set of large %d composed by ', length(setM1));
        for j = 1:length(setM1)
            if j == 1
                fprintf('%s', setM1{j});
            elseif j == length(setM1)
                fprintf(' and %s', setM1{j});
            else
                fprintf(', %s', setM1{j});
            end
        end
        fprintf('\n');
        fprintf('The production of 2,3-butanediol after optimization is %.3f mmol/gCDW/h \n', productFluxM1);
        fprintf('The growth rate after optimization is %.3f /h \n', growthRateM1);
        fprintf('...Performing coupling analysis...\n');
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, setM1, product);
        fprintf('The solution is of type: %s\n', type);
        fprintf('The maximun growth rate given the optKnock set is %.3f /h \n', maxGrowth);
        fprintf(['The maximun and minimun production of 2,3-butanediol given the optKnock set is ' ...
                 '%.3f and %.3f mmol/gCDW/h, respectively \n\n'], minProd, maxProd);
        if strcmp(type, 'growth coupled')
            singleProductionEnvelope(model, setM1, product, biomass, 'savePlot', 0, 'showPlot', 1);
        end
    else
        if nIter == 1
            fprintf('optKnock was not able to found an optKnock set\n');
        else
            fprintf('optKnock was not able to found additional optKnock sets\n');
        end
        break;
    end
    nIter = nIter + 1;
end

