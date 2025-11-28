% MIT License
% 
% Copyright (c) 2025 Jecricho Knox
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clc,clear
% 高斯-勒让德龙格库塔法求解开普勒问题（变步长版本）
% 使用局部外推法进行步长控制

% 物理参数
mu = 1.0;      % 引力参数 (GM)
e = 0.8;       % 偏心率
a = 20;       % 半长轴
c0 = e*a;

% 初始条件 (椭圆轨道)
q0 = [a*(1-e); 0];                    % 初始位置
p0 = [0; sqrt((mu/a)*(1+e)/(1-e))];   % 初始动量
y0 = [q0; p0];                        % 初始状态向量 [q; p]

% 时间参数
T = 2*pi*sqrt(a^3/mu); % 轨道周期
tspan = [0, 500*T];    % 模拟500个周期

% 变步长控制参数
h0 = T/200;           % 初始步长
h_min = T/1e6;        % 最小步长
h_max = T/10;         % 最大步长
tol = 1e-14;          % 误差容限
safe_fac = 0.7;       % 安全系数
order = 12;           % GL6方法是12阶方法

% 高斯-勒让德 RK6 系数
[c,A,b] = Guass_legendre_butcher_6();

% 初始化存储变量
t_list = tspan(1);
y_list = y0;
h_list = h0;

% 存储能量和角动量误差
E_error_list = 0;
L_error_list = 0;

% 精确的能量和角动量值
E_exact = 0.5*dot(p0, p0) - mu/norm(q0);
L_exact = q0(1)*p0(2) - q0(2)*p0(1);

% 主循环 - 变步长积分
tic
t_current = tspan(1);
y_current = y0;
h_current = h0;

iteration_count = 0;
step_rejected = 0;

while t_current < tspan(2)
    iteration_count = iteration_count + 1;
    
    % 限制步长不超过终点
    if t_current + h_current > tspan(2)
        h_current = tspan(2) - t_current;
    end
    
    % 使用局部外推法：计算全步长和两个半步长
    [y1_full, success] = glrk_step_improved(t_current, y_current, h_current, c, A, b, mu);
    
    if ~success
        % 如果全步长失败，直接拒绝并缩小步长
        h_current = max(h_min, h_current * 0.5);
        step_rejected = step_rejected + 1;
        continue;
    end
    
    % 计算两个半步长
    [y_half, success1] = glrk_step_improved(t_current, y_current, h_current/2, c, A, b, mu);
    if ~success1
        h_current = max(h_min, h_current * 0.5);
        step_rejected = step_rejected + 1;
        continue;
    end
    
    [y2_full, success2] = glrk_step_improved(t_current + h_current/2, y_half, h_current/2, c, A, b, mu);
    if ~success2
        h_current = max(h_min, h_current * 0.5);
        step_rejected = step_rejected + 1;
        continue;
    end
    
    % 计算误差估计（使用两个结果的差异）
    error_est = norm(y1_full - y2_full);
    
    % 计算缩放后的误差（考虑相对容差）
    scale = max(norm(y_current), 1.0);
    scaled_error = error_est / (scale * tol);
    
    if scaled_error <= 1.0  % 接受这一步
        % 使用更精确的两个半步长结果（局部外推）
        y_next = y2_full;
        t_next = t_current + h_current;
        
        % 存储结果
        t_list(end+1) = t_next;
        y_list(:, end+1) = y_next;
        h_list(end+1) = h_current;
        
        % 计算守恒量误差
        q = y_next(1:2);
        p = y_next(3:4);
        E_current = 0.5*dot(p, p) - mu/norm(q);
        L_current = q(1)*p(2) - q(2)*p(1);
        
        E_error_list(end+1) = abs(E_current - E_exact);
        L_error_list(end+1) = abs(L_current - L_exact);
        
        % 调整下一步步长
        if error_est > 0
            h_new = h_current * safe_fac * (1.0 / scaled_error)^(1/(order+1));
        else
            h_new = h_current * 2.0;  % 如果误差为0，可以大胆增加步长
        end
        
        h_current = min(max(h_new, h_min), h_max);
        t_current = t_next;
        y_current = y_next;
        
    else  % 拒绝这一步，缩小步长重试
        h_current = max(h_min, h_current * safe_fac * (1.0 / scaled_error)^(1/(order+1)));
        step_rejected = step_rejected + 1;
    end
    
    % 防止无限循环
    if iteration_count > 100000
        warning('达到最大迭代次数，可能步长过小');
        break;
    end
end
toc

fprintf('总步数: %d\n', length(t_list));
fprintf('拒绝步数: %d\n', step_rejected);
fprintf('接受率: %.2f%%\n', (length(t_list)-1)/(length(t_list)-1+step_rejected)*100);

% 转换为固定格式的输出
t = t_list;
y = y_list;
h_used = h_list;

% 提取位置坐标
qx = y(1, :);
qy = y(2, :);

% 绘图
figure(Color='w');

% 子图1: 轨迹
subplot(2,2,1);
scatter(qx, qy, 10, 'filled', 'MarkerEdgeColor','none');
hold on;
plot(qx(1), qy(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x');
ylabel('y');
title(sprintf('开普勒问题(e=%.1f)轨迹（变步长）',e));
axis equal;
grid on;
legend('轨迹', '起点', '中心天体', 'Location', 'best');

% 子图2: 使用的步长
subplot(2,2,2);
semilogy(t(1:end-1), h_used(2:end), 'g-', 'LineWidth', 1);
xlabel('时间（t）');
ylabel('步长');
title('步长变化');
grid on;

% 子图3: 能量和角动量误差（对数坐标）
subplot(2,2,3);
semilogy(t, E_error_list, 'r-', 'LineWidth', 1.5, 'DisplayName', '总能量');
hold on;
semilogy(t, L_error_list, 'b-', 'LineWidth', 1.5, 'DisplayName', '角动量');
xlabel('时间（t）');
ylabel('相对误差y对数');
legend('show', 'Location','southeast');
grid on;

% 子图4: 轨迹的椭圆偏差
subplot(2,2,4);
semilogy(t, abs((qx+c0).^2 + qy.^2./(1-e^2) - a^2), 'm-', 'LineWidth', 1);
xlabel('时间（t）');
ylabel('$\frac{x^2}{a^2} + \frac{y^2}{b^2} - 1$','Interpreter','latex','FontSize',12);
title('轨道的椭圆偏差');
grid on;

% 显示统计信息
fprintf('平均能量误差: %.4e\n', mean(E_error_list));
fprintf('平均角动量误差: %.4e\n', mean(L_error_list));
fprintf('最大步长: %.4e, 最小步长: %.4e\n', max(h_used), min(h_used));

function [y_next, success] = glrk_step_improved(t, y, h, c, A, b, mu)
    % 改进的高斯勒让德方法单步积分 - 使用高斯-塞德尔风格迭代
    % 返回: y_next - 下一步状态, success - 是否成功
    
    max_iter = 25;  % 可以减少迭代次数，因为收敛更快
    tol_iter = 5e-15;
    max_tol = 5e1;
    s = length(c);
    n = length(y);
    
    % 初始化阶段值（使用当前y作为初始猜测）
    Y = repmat(y, 1, s);
    F = zeros(n, s);
    
    success = true;
    
    % 先计算初始的F值
    for i = 1:s
        F(:, i) = kepler_rhs(t + c(i)*h, Y(:, i), mu);
    end
    
    % 高斯-塞德尔风格迭代
    for iter = 1:max_iter
        Y_old = Y;
        
        % 顺序更新每个阶段值，立即使用最新信息
        for i = 1:s
            % 计算当前阶段的更新
            sum_AF = zeros(n, 1);
            for j = 1:s
                sum_AF = sum_AF + A(i,j) * F(:, j);
            end
            Y(:, i) = y + h * sum_AF;
            
            % 立即更新当前阶段的导数
            F(:, i) = kepler_rhs(t + c(i)*h, Y(:, i), mu);
        end
        
        % 检查收敛性
        max_change = 0;
        for i = 1:s
            max_change = max(max_change, norm(Y(:, i) - Y_old(:, i))/norm(Y_old(:, i)));
        end
        
        if max_change < tol_iter || max_change > max_tol
            break;
        end
    end
    
    if max_change > max_tol
        warning('在t=%.8f处迭代不收敛，\n当前步长：%.2e，最大变化: %.2e', t, h, max_change);
        success = false;
    end

    if iter == max_iter && max_change > tol_iter
        warning('在t=%.8f处迭代未完全收敛，\n当前步长为：%.2e，最大变化: %.2e', t, h, max_change);
        success = false;
    end
    
    % 计算下一步
    sum_bF = zeros(n, 1);
    for i = 1:s
        sum_bF = sum_bF + b(i) * F(:, i);
    end
    y_next = y + h * sum_bF;
end

function dydt = kepler_rhs(t, y, mu)
    % 开普勒问题的右手边函数
    q = y(1:2);
    p = y(3:4);
    
    r = norm(q);
    a = -mu * q / r^3;
    
    dydt = [p; a];
end

