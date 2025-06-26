% =========================================================================
% MATLAB 程序：蠕变数据多模型拟合分析
% =========================================================================
% 功能：
% 1. 弹窗让用户选择包含蠕变数据的 CSV 文件。
% 2. 自动识别文件中的不同工况（例如，不同温度）。
% 3. 对每个工况的数据，使用多种唯象模型进行非线性最小二乘拟合。
%    - Burger's Model (4-Element)
%    - Findley Power Law Model
%    - Andrade's Model
% 4. 在命令窗口中输出每个模型的拟合参数和R-squared值。
% 5. 为每个工況生成图表，对比原始数据和所有模型的拟合曲线。
%
% 要求:
% - MATLAB R2016b 或更高版本。
% - Curve Fitting Toolbox™。
%
% CSV文件格式要求:
% - 必须包含列标题。
% - 必须包含 'Time' 和 'Displacement' 列。
% - 可选包含 'Temperature' 或其他分组变量列。
% =========================================================================

%% 1. 初始化和用户输入
clear; 
clc; 
close all;

fprintf('蠕变数据多模型拟合程序\n');
fprintf('---------------------------------\n');

% --- 弹窗选择文件 ---
[file, path] = uigetfile('*.csv', '请选择包含蠕变数据的CSV文件');
if isequal(file, 0)
    disp('用户取消了操作');
    return;
end
fullFilePath = fullfile(path, file);
fprintf('已选择文件: %s\n', fullFilePath);

% --- 读取数据 ---
try
    data = readtable(fullFilePath);
catch ME
    error('无法读取CSV文件。请检查文件格式是否正确。错误信息: %s', ME.message);
end

% --- 获取列名 ---
colNames = data.Properties.VariableNames;
fprintf('\n文件包含的列: %s\n', strjoin(colNames, ', '));

% 确定时间、位移和分组变量列
timeCol = 'Time';
dispCol = 'Displacement';
groupCol = ''; % 分组变量，例如 'Temperature' 或 'Stress'

% 自动检测分组变量（如果存在）
if ismember('Temperature', colNames)
    groupCol = 'Temperature';
elseif ismember('Stress', colNames)
    groupCol = 'Stress';
end

if ~ismember(timeCol, colNames) || ~ismember(dispCol, colNames)
    error("CSV文件必须包含 '%s' 和 '%s' 列。", timeCol, dispCol);
end

if ~isempty(groupCol)
    fprintf("将使用 '%s' 列对数据进行分组拟合。\n", groupCol);
    uniqueGroups = unique(data.(groupCol));
else
    fprintf("未找到分组列，将所有数据作为一个整体进行拟合。\n");
    uniqueGroups = 1; % 虚拟分组
    data.Group = ones(height(data), 1); % 创建虚拟分组列
    groupCol = 'Group';
end

fprintf('---------------------------------\n\n');

%% 2. 定义蠕变模型

% --- 模型定义 (fittype格式) ---
% a) Burger's Model (4-element)
% y(t) = p1 + p2*(1-exp(-t/p3)) + p4*t
% p1: 瞬时弹性应变 (或位移)
% p2: 粘弹性应变 (或位移)
% p3: 粘弹性延迟时间 (Kelvin-Voigt模型)
% p4: 稳态蠕变速率 (Maxwell模型)
burgerModel = fittype('p1 + p2 * (1 - exp(-t/p3)) + p4 * t', ...
    'independent', 't', 'coefficients', {'p1', 'p2', 'p3', 'p4'});

% b) Findley Power Law Model
% y(t) = y0 + A * t^n
% y0: 瞬时应变 (或位移)
% A: 幅值系数
% n: 时间指数 (通常 < 1)
findleyModel = fittype('y0 + A * t.^n', ...
    'independent', 't', 'coefficients', {'y0', 'A', 'n'});

% c) Andrade's Creep Law
% y(t) = A * t^(1/3) + B*t
% A*t^(1/3): 瞬态蠕变项
% B*t: 稳态蠕变项
andradeModel = fittype('A * t.^(1/3) + B * t',...
    'independent', 't', 'coefficients', {'A', 'B'});


models = {burgerModel, findleyModel, andradeModel};
modelNames = {'Burger (4-Element)', 'Findley Power Law', 'Andrade'};

%% 3. 循环处理每个数据组并进行拟合

% 循环遍历每个独特的分组（例如，每个温度）
for i = 1:length(uniqueGroups)
    
    currentGroupVal = uniqueGroups(i);
    
    % --- 提取当前组的数据 ---
    groupData = data(data.(groupCol) == currentGroupVal, :);
    t = groupData.(timeCol);
    y = groupData.(dispCol);
    
    % 确保数据从(0,0)或接近0开始，以获得更好的拟合效果
    y = y - y(1);
    t = t - t(1);
    
    fprintf('=== 开始处理工况: %s = %.2f ===\n', groupCol, currentGroupVal);
    
    % --- 创建新的图表 ---
    figure('Name', sprintf('蠕变拟合: %s = %.2f', groupCol, currentGroupVal), 'NumberTitle', 'off');
    hold on;
    grid on;
    
    % 绘制原始数据点
    plot(t, y, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', '原始数据');
    
    % --- 存储结果 ---
    fitResults = cell(length(models), 1);
    gofResults = cell(length(models), 1);
    
    % --- 对当前数据组应用所有模型 ---
    for j = 1:length(models)
        
        fprintf('  正在使用模型: %s\n', modelNames{j});
        
        % --- 设置拟合选项 (初始值和边界) ---
        % 合理的初始值是拟合成功的关键！
        y_max = max(y);
        t_max = max(t);
        
        opts = fitoptions(models{j});
        
        % 根据不同模型设置不同的初始值和边界
        switch modelNames{j}
            case 'Burger (4-Element)'
                opts.StartPoint = [0, y_max/2, t_max/4, y_max/t_max];
                opts.Lower = [0, 0, 0, 0];
                opts.Upper = [y_max/2, y_max*2, t_max*2, (y_max/t_max)*5];
            case 'Findley Power Law'
                opts.StartPoint = [0, y_max/2, 0.3];
                opts.Lower = [0, 0, 0];
                opts.Upper = [y_max/2, y_max*5, 1];
            case 'Andrade'
                 opts.StartPoint = [y_max/2, y_max/t_max];
                 opts.Lower = [0, 0];
                 opts.Upper = [y_max*2, (y_max/t_max)*5];
        end

        % --- 执行拟合 ---
        try
            [fit_obj, gof] = fit(t, y, models{j}, opts);
            fitResults{j} = fit_obj;
            gofResults{j} = gof;
            
            % 在图上绘制拟合曲线
            plot(fit_obj, 'predobs'); % 'predobs' 会同时绘制置信区间
            
        catch ME
            fprintf('    ** 模型 "%s" 拟合失败: %s **\n', modelNames{j}, ME.message);
            fitResults{j} = [];
            gofResults{j} = [];
        end
    end
    
    % --- 美化图表并输出结果 ---
    title(sprintf('蠕变位移 vs. 时间 (%s = %.2f)', groupCol, currentGroupVal));
    xlabel('时间 (s)');
    ylabel('位移 (mm)');
    legend('Location', 'NorthWest');
    hold off;
    
    % 在命令窗口打印该工况的拟合参数总结
    fprintf('\n  --- 拟合结果总结: %s = %.2f ---\n', groupCol, currentGroupVal);
    for j = 1:length(models)
        if ~isempty(fitResults{j})
            fprintf('  模型: %s\n', modelNames{j});
            disp(fitResults{j});
            fprintf('    R-squared: %.4f\n\n', gofResults{j}.rsquare);
        end
    end
    fprintf('===============================================\n\n');
    
end

disp('所有工况处理完毕！');