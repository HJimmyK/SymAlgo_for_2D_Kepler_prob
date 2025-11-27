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

function [t, y, vy] = rkn76t_adaptive(f, tspan, y0, vy0, tol)
% RKN76T_ADAPTIVE 自适应步长 Runge-Kutta-Nyström 7(6)T 方法
%   用于求解二阶常微分方程 y'' = f(t, y)
%
% 输入:
%   f    - 函数句柄，表示 y'' = f(t, y)
%   tspan - 时间区间 [t0, tf]
%   y0    - 初始位置
%   vy0   - 初始速度
%   tol   - 容差参数
%
% 输出:
%   t     - 时间点
%   y     - 位置解
%   vy    - 速度解

    % 初始化Butcher表系数
    [c, A, b, b_hat, b_dot_hat] = rkn76t_coefficients();
    
    % 初始化变量
    t0 = tspan(1);
    tf = tspan(2);
    h = (tf - t0) * tol;  % 初始步长
    safety = .97;          % 安全系数
    min_h = (tf - t0) * tol * 1e-1; % 最小步长
    
    % 预分配存储（可根据需要调整大小）
    chunk_size = 1000;
    t = zeros(1, chunk_size);
    y = zeros(length(y0), chunk_size);
    vy = zeros(length(vy0), chunk_size);
    
    % 设置初始条件
    t(1) = t0;
    y(:, 1) = y0(:);
    vy(:, 1) = vy0(:);
    
    % 主循环
    n = 1;
    while t(n) < tf
        % 确保不超出最终时间
        if t(n) + h > tf
            h = tf - t(n);
        end
        
        % 计算RKN阶段
        [y_new, vy_new, err] = rkn76t_step(f, t(n), y(:, n), vy(:, n), h, c, A, b, b_hat, b_dot_hat);
        
        % 检查误差是否可接受
        if err <= tol
            % 接受步长
            n = n + 1;
            
            % 如果需要，扩展存储
            if n > size(t, 2)
                t = [t, zeros(1, chunk_size)];
                y = [y, zeros(length(y0), chunk_size)];
                vy = [vy, zeros(length(vy0), chunk_size)];
            end
            
            % 存储结果（使用高阶估计）
            t(n) = t(n-1) + h;
            y(:, n) = y_new;
            vy(:, n) = vy_new;
        end
        
        h_new = h * safety * (tol / err)^(1/7);  % 7阶方法
        
        % 限制步长变化
        h_new = min(max(h_new, min_h), 2*h);
        h = h_new;
    end
    
    % 裁剪输出数组
    t = t(1:n);
    y = y(:, 1:n);
    vy = vy(:, 1:n);
end

function [c, A, b, b_hat, b_dot_hat] = rkn76t_coefficients()
    % RKN76T_COEFFICIENTS 返回RKN7(6)T方法的Butcher表系数
    
    % 阶段数
    s = 9;
    
    % 时间节点
    c = [0; 1/10; 1/5; 3/8; 1/2; (7-sqrt(21))/14; (7+sqrt(21))/14; 1; 1];
    
    % 矩阵A（下三角）
    A = zeros(s, s-1);
    A(2,1) = 1/200;
    A(3,1) = 1/150;       A(3,2) = 1/75;
    A(4,1) = 171/8192;    A(4,2) = 45/4096;     A(4,3) = 315/8192;
    A(5,1) = 5/288;       A(5,2) = 25/528;      A(5,3) = 25/672;      A(5,4) = 16/693;
    
    % Stage 6
    A(6,1) = (1003 - 205*sqrt(21))/12348;
    A(6,2) = -25*(751 - 173*sqrt(21))/90552;
    A(6,3) = 25*(624 - 137*sqrt(21))/43218;
    A(6,4) = -128*(361 - 79*sqrt(21))/237699;
    A(6,5) = (3411 - 745*sqrt(21))/24696;
    
    % Stage 7
    A(7,1) = (793 + 187*sqrt(21))/12348;
    A(7,2) = -25*(331 + 113*sqrt(21))/90552;
    A(7,3) = 25*(1044 + 247*sqrt(21))/43218;
    A(7,4) = -128*(14885 + 3779*sqrt(21))/9745659;
    A(7,5) = (3327 + 797*sqrt(21))/24696;
    A(7,6) = -(581 + 127*sqrt(21))/1722;
    
    % Stage 8
    A(8,1) = -(157 - 3*sqrt(21))/378;
    A(8,2) = 25*(143 - 10*sqrt(21))/2772;
    A(8,3) = -25*(876 + 55*sqrt(21))/3969;
    A(8,4) = 1280*(913 + 18*sqrt(21))/596673;
    A(8,5) = -(1353 + 26*sqrt(21))/2268;
    A(8,6) = 7*(1777 + 377*sqrt(21))/4428;
    A(8,7) = 7*(5 - sqrt(21))/36;
    
    % Stage 9
    A(9,1) = 1/20;
    A(9,5) = 8/45;
    A(9,6) = 7*(7+sqrt(21))/360;
    A(9,7) = 7*(7-sqrt(21))/360;
    
    % 6阶权重
    b = [1/20; 0; 0; 0; 8/45; 7*(7+sqrt(21))/360; 7*(7-sqrt(21))/360; 0 ;0];
    
    % 7阶权重
    b_hat = [1/20; 0; 0; 0; 8/45; 7*(7+sqrt(21))/360; 7*(7-sqrt(21))/360; -1/20; 1/20];

    b_dot_hat = [1/20; 0; 0; 0; 16/45; 49/180; 49/180; 1/20; 0];
end

function [y_new, vy_new, err] = rkn76t_step(f, t, y, vy, h, c, A, b, b_hat, b_dot_hat)
    % RKN76T_STEP 执行单步RKN7(6)T计算
    
    s = length(c);  % 阶段数
    n = length(y);  % 维度
    
    % 初始化阶段变量
    k = zeros(n, s);
    
    % 计算各个阶段
    for i = 1:s
        % 计算当前阶段的位置和速度
        y_stage = y + h * c(i) * vy;
        for j = 1:i-1
            y_stage = y_stage + h^2 * A(i,j) * k(:,j);
        end
        
        % 计算当前阶段的加速度
        k(:,i) = f(t + c(i)*h, y_stage);
    end
    
    % 计算7阶解（高阶估计）
    y_new = y + h * vy;
    vy_new = vy;
    for i = 1:s
        y_new = y_new + h^2 * b(i) * k(:,i);
        vy_new = vy_new + h * b_dot_hat(i) * k(:,i);
    end
    
    % 计算误差估计（基于位置和速度的差异）
    err = h^2*norm(k(:,s-1)-k(:,s))/20;
end