% 清除工作区和命令行
clear; clc; close all;

%% 1. 加载 T2 Solid Echo 数据 (CSV文件)
% 用户可以通过文件选择对话框选择CSV文件。
% 假设CSV文件包含两列：第一列是时间 (例如，微秒 us)，第二列是信号强度。
% 请确保你的时间单位与后续的T2_bins的单位保持一致！

fprintf('请选择你的 T2 Solid Echo 数据 CSV 文件...\n');
[filename, pathname] = uigetfile({'*.csv', 'CSV Files (*.csv)'}, '选择 T2 Solid Echo 数据文件');

if isequal(filename, 0)
    disp('用户取消了文件选择。');
    return; % 如果用户取消，则退出脚本
else
    filepath = fullfile(pathname, filename);
    fprintf('正在加载文件: %s\n', filepath);
    
    try
        % 尝试读取CSV文件。假设第一列是时间，第二列是信号强度。
        % 如果你的CSV文件有表头，或者列的顺序不同，你需要调整 csvread 或 readtable 的参数。
        data = csvread(filepath); % 适用于纯数字CSV，无表头
        
        t = data(:, 1);    % 时间向量，第一列
        M_data = data(:, 2); % 信号强度向量，第二列
        
        % 确保信号非负 (物理约束)
        M_data(M_data < 0) = 0; 

        fprintf('数据加载完成。时间点数量: %d\n', length(t));
        
    catch ME
        warning('文件读取或数据处理出错: %s\n', ME.message);
        disp('请检查CSV文件格式是否正确，并确保只包含时间和信号强度两列数字数据。');
        return;
    end
end

% 请注意：以下模拟数据生成部分已被注释掉，因为我们将使用加载的真实数据
% %% 1. 模拟 T2 Solid Echo 数据 (如果已有真实数据，请跳过此步骤并加载数据)
% num_points = 100; 
% max_time = 1000; % us
% t = linspace(10, max_time, num_points)'; 
% T2_component1 = 20; Amplitude1 = 0.6; 
% T2_component2 = 200; Amplitude2 = 0.4; 
% M_ideal = Amplitude1 * exp(-t / T2_component1) + ...
%           Amplitude2 * exp(-t / T2_component2);
% noise_level = 0.02; 
% M_data = M_ideal + noise_level * randn(size(M_ideal));
% M_data(M_data < 0) = 0; 
% fprintf('模拟数据生成完成。\n');


%% 2. 定义 T2 值范围和离散化
% 这是NNLS反演的关键参数。T2值的范围应涵盖样品中所有可能的T2值。
% 减少 num_T2_bins 的值会使每个T2 bin更宽，从而导致输出的峰更宽。
% 增加此值会使曲线更平滑，但峰可能显得更尖锐。

T2_min_log = log10(1);   % 最小T2值的对数 (us)
T2_max_log = log10(10000); % T2值的最大对数值更改为 log10(10000) = 4
num_T2_bins = 200;       % 为了更好地观察正则化效果，我们将其设回200。
                         % 正则化本身会带来平滑效果，所以不用刻意减少这里的值。

% 在对数尺度上均匀分布T2值，然后转换回线性尺度
T2_bins_log = linspace(T2_min_log, T2_max_log, num_T2_bins);
T2_bins = 10.^T2_bins_log; % T2值向量

fprintf('T2 离散化范围: [%.1f us, %.1f us], 离散点数: %d\n', min(T2_bins), max(T2_bins), num_T2_bins);


%% 3. 构建核矩阵 A
% 核矩阵 A 的大小是 (num_points x num_T2_bins)
% A(j, i) = exp(-t(j) / T2_bins(i))

A = zeros(length(t), num_T2_bins); % 使用实际加载的数据长度
for i = 1:num_T2_bins
    A(:, i) = exp(-t / T2_bins(i));
end

fprintf('核矩阵 A 构建完成，大小: %d x %d\n', size(A,1), size(A,2));


%% 4. 使用 NNLS 算法进行反演 (加入 Tikhonov 正则化)

% 定义正则化参数 (lambda)
% 这个值非常重要，需要根据数据和期望的平滑程度进行调整。
% 尝试从一个小值开始（例如 1e-3, 1e-2）并逐渐增大，观察T2分布谱的变化。
% 对于不同的数据集，最佳的lambda值可能不同。
lambda = 1; % 调整此值以控制平滑程度

fprintf('正则化参数 lambda = %.4f\n', lambda);

% 构建二阶导数矩阵 L (用于Tikhonov正则化)
% L是一个 (num_T2_bins - 2) x num_T2_bins 的矩阵，它近似计算二阶导数
L = zeros(num_T2_bins - 2, num_T2_bins);
for i = 1:(num_T2_bins - 2)
    L(i, i)     = 1;
    L(i, i + 1) = -2;
    L(i, i + 2) = 1;
end

% 构建增广矩阵 A_prime 和 增广向量 M_prime
A_prime = [A; lambda * L];
M_prime = [M_data; zeros(size(L, 1), 1)]; % 正则化项对应的目标值为0

fprintf('开始带正则化的 NNLS 反演...\n');
% 使用 lsqnonneg 求解增广问题，它仍然会确保解是非负的
F_T2_distribution = lsqnonneg(A_prime, M_prime);
fprintf('带正则化的 NNLS 反演完成。\n');


% 归一化 T2 分布谱
% 这里我们按照常用方法，将其总和归一化为1。
area_T2_distribution = sum(F_T2_distribution);
F_T2_distribution_normalized = F_T2_distribution / area_T2_distribution;
fprintf('T2 分布谱已归一化，总和为: %.4f\n', sum(F_T2_distribution_normalized));


%% 5. 绘制结果
figure; % 创建一个新图窗

% 子图 1: 原始数据与拟合曲线 (时间轴使用对数坐标，因为它更适合衰减曲线)
subplot(2,1,1);
plot(t, M_data, 'o', 'DisplayName', '原始数据');
hold on;
M_fit = A * F_T2_distribution; % 注意：这里仍然使用原始的 A 来计算拟合曲线，而不是 A_prime
plot(t, M_fit, '-', 'LineWidth', 1.5, 'DisplayName', 'NNLS拟合曲线 (正则化)');
title('T2 Solid Echo 衰减曲线与 NNLS 拟合');
xlabel('时间 (us)');
ylabel('信号强度');
legend('show');
grid on;
set(gca, 'XScale', 'log'); % 衰减曲线通常在对数时间轴上显示更好
xlim([min(t) max(t)]);

% 子图 2: 连续 T2 分布谱 (贡献度归一化，T2弛豫时间采用对数坐标)
subplot(2,1,2);
% 使用 plot 函数绘制实线曲线，使用归一化后的贡献度，并设置对数X轴。
plot(T2_bins, F_T2_distribution_normalized, 'b-', 'LineWidth', 1.5, 'DisplayName', 'T2 分布 (归一化, 正则化)'); 

set(gca, 'XScale', 'log'); % 将T2弛豫时间轴重新设置为对数坐标
title(['归一化 T2 分布谱 (对数 T2 坐标, \lambda = ', num2str(lambda, '%.1e'), ')']); % 标题中显示lambda值
xlabel('T2 弛豫时间 (us)');
ylabel('归一化相对贡献度'); % 标签改为“归一化相对贡献度”
xlim([min(T2_bins) max(T2_bins)]); % 根据实际T2范围调整X轴显示
grid on;
ylim([0 inf]); % 确保Y轴从0开始，自动调整上限


% 调整布局
sgtitle('T2 Solid Echo 数据反演分析'); % 总标题


%% 6. 输出拟合分布谱数据到CSV文件

% 定义输出文件名
% 可以在原始CSV文件名的基础上添加后缀，或者定义一个新名称
[~, name, ~] = fileparts(filename); % 获取原始文件名（不包含扩展名）
output_filename = [name '_T2_distribution_regularized.csv']; % 输出文件添加 regularized 后缀

fprintf('\n正在将 T2 分布谱数据输出到文件: %s\n', output_filename);

% 将 T2_bins 和 F_T2_distribution_normalized 合并为一个矩阵
% 使用 (:) 操作符确保它们是列向量
output_data = [T2_bins(:), F_T2_distribution_normalized(:)];

% 推荐使用 writetable，因为它支持写入列名并且更灵活
T = array2table(output_data, 'VariableNames', {'T2_Relaxation_Time_us', 'Normalized_Contribution'});
writetable(T, output_filename);

fprintf('T2 分布谱数据已成功输出。\n');