%% startme
% 加载COBRA和RAVEN工具包，这样才能正常使用这两个工具包，即调用这两个工具包的所有函数

cd ../cobratoolbox;                     % 进入COBRA文件夹（如果路径不对需要手动修改至相应路径）
initCobraToolbox;                       % 初始化COBRA
cd ../Tutorials-of-metabolic-modeling/; % 从COBRA文件夹回到本教程的文件夹
addpath(genpath('../RAVEN/'));          % 添加RAVEN所在路径（如果路径不对需要手动修改至相应路径）
